"""Automatic style check for the repository (standard library only).

Usage: python tests/style/check_style.py [ROOT]

Prints one `path:line: rule: message` line per violation and exits 1 when any exist.
"""

import ast
import re
import sys
from pathlib import Path

SCANNED_DIRECTORIES = ("src", "include", "cpp", "tests", "python/src", "python/tests", "python/benchmarks")
IGNORED_DIRECTORIES = ("build", "python/.venv")
CPP_SUFFIXES = {".c", ".h", ".cc", ".hh", ".cpp", ".hpp", ".cxx"}
MODULE_SUFFIXES = CPP_SUFFIXES | {".py"}
MAX_MODULES_PER_FOLDER = 8
MAX_FUNCTION_LINES = 60
MAX_FUNCTION_PARAMS = 5
MAX_NESTING = 3

NON_CODE = re.compile(
    r"""R"(?P<delim>[^()\s]{0,16})\(.*?\)(?P=delim)"
    |//[^\n]*
    |/\*.*?\*/
    |"(?:\\.|[^"\\\n])*"
    |'(?:\\.|[^'\\\n])*'""",
    re.DOTALL | re.VERBOSE,
)
ELSE_TOKEN = re.compile(r"\belse\b")
ELSE_KEYWORD_LINE = re.compile(r"else\s*:|elif\b")


def ast_nodes(*names):
    return tuple(getattr(ast, name) for name in names if hasattr(ast, name))


BLOCK_NODES = ast_nodes("If", "For", "AsyncFor", "While", "With", "AsyncWith", "Try", "TryStar", "Match")
FUNCTION_NODES = (ast.FunctionDef, ast.AsyncFunctionDef)
ELSE_BEARING_NODES = ast_nodes("If", "For", "AsyncFor", "While", "Try", "TryStar")


class Violation:
    def __init__(self, path, line, rule, message):
        self.path = path
        self.line = line
        self.rule = rule
        self.message = message

    def __str__(self):
        return f"{self.path}:{self.line}: {self.rule}: {self.message}"


def blank_out(match):
    return re.sub(r"[^\n]", " ", match.group(0))


def strip_non_code(text):
    return NON_CODE.sub(blank_out, text)


def check_cpp_else(path, text):
    violations = []
    for number, line in enumerate(strip_non_code(text).splitlines(), start=1):
        if line.lstrip().startswith("#"):
            continue
        if ELSE_TOKEN.search(line):
            violations.append(Violation(path, number, "else", "`else` is forbidden; use a guard clause or early return"))
    return violations


def else_line(lines, node):
    first = node.orelse[0]
    for number in range(first.lineno, 0, -1):
        stripped = lines[number - 1].lstrip()
        if ELSE_KEYWORD_LINE.match(stripped):
            return number
    return first.lineno


def check_python_else(path, tree, lines):
    violations = []
    for node in ast.walk(tree):
        if isinstance(node, ELSE_BEARING_NODES) and node.orelse:
            message = "`else`/`elif` is forbidden; use a guard clause, early return or a table"
            violations.append(Violation(path, else_line(lines, node), "else", message))
    return violations


def is_simple_operand(node):
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return isinstance(node.operand, ast.Constant)
    if isinstance(node, ast.Attribute):
        return is_attribute_chain(node.value)
    return isinstance(node, (ast.Name, ast.Constant))


def is_attribute_chain(node):
    if isinstance(node, ast.Attribute):
        return is_attribute_chain(node.value)
    return isinstance(node, ast.Name)


def ternary_problem(node):
    if node.end_lineno != node.lineno:
        return "spans more than one line"
    nested = [child for child in ast.walk(node) if isinstance(child, ast.IfExp) and child is not node]
    if nested:
        return "is nested"
    if not (is_simple_operand(node.body) and is_simple_operand(node.orelse)):
        return "has a branch that is not a name, literal or attribute"
    return None


def check_python_ternaries(path, tree):
    violations = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.IfExp):
            continue
        problem = ternary_problem(node)
        if problem:
            message = f"ternary {problem}; use `if` with an early return"
            violations.append(Violation(path, node.lineno, "ternary", message))
    return violations


def function_line_count(node, lines):
    skipped = set()
    first_statement = node.body[0]
    is_docstring = isinstance(first_statement, ast.Expr) and isinstance(first_statement.value, ast.Constant)
    if is_docstring and isinstance(first_statement.value.value, str):
        skipped = set(range(first_statement.lineno, first_statement.end_lineno + 1))
    count = 0
    for number in range(node.lineno, node.end_lineno + 1):
        stripped = lines[number - 1].strip()
        if stripped and not stripped.startswith("#") and number not in skipped:
            count += 1
    return count


def parameter_count(node):
    arguments = node.args
    named = arguments.posonlyargs + arguments.args
    if named and named[0].arg in ("self", "cls"):
        named = named[1:]
    extra = [arguments.vararg, arguments.kwarg]
    return len(named) + len(arguments.kwonlyargs) + sum(1 for item in extra if item is not None)


def deepest_block(node, depth):
    deepest = (depth, node)
    for child in ast.iter_child_nodes(node):
        if isinstance(child, FUNCTION_NODES + (ast.ClassDef, ast.Lambda)):
            continue
        child_depth = depth
        if isinstance(child, BLOCK_NODES):
            child_depth += 1
        candidate = deepest_block(child, child_depth)
        if candidate[0] > deepest[0]:
            deepest = candidate
    return deepest


def check_function(path, node, lines):
    violations = []
    length = function_line_count(node, lines)
    if length > MAX_FUNCTION_LINES:
        message = f"`{node.name}` has {length} lines (max {MAX_FUNCTION_LINES})"
        violations.append(Violation(path, node.lineno, "function-length", message))
    params = parameter_count(node)
    if params > MAX_FUNCTION_PARAMS:
        message = f"`{node.name}` has {params} parameters (max {MAX_FUNCTION_PARAMS})"
        violations.append(Violation(path, node.lineno, "function-params", message))
    depth, block = deepest_block(node, 0)
    if depth > MAX_NESTING:
        message = f"`{node.name}` nests {depth} blocks (max {MAX_NESTING})"
        violations.append(Violation(path, block.lineno, "nesting", message))
    return violations


def check_python(path, text):
    try:
        tree = ast.parse(text)
    except SyntaxError as error:
        return [Violation(path, error.lineno or 1, "syntax", error.msg)]
    lines = text.splitlines()
    violations = check_python_else(path, tree, lines) + check_python_ternaries(path, tree)
    for node in ast.walk(tree):
        if isinstance(node, FUNCTION_NODES):
            violations.extend(check_function(path, node, lines))
    return violations


def check_file(path, relative):
    text = path.read_text(encoding="utf-8", errors="replace")
    if path.suffix == ".py":
        return check_python(relative, text)
    return check_cpp_else(relative, text)


def is_ignored(relative):
    return any(relative == ignored or relative.startswith(ignored + "/") for ignored in IGNORED_DIRECTORIES)


def scanned_files(root):
    for directory in SCANNED_DIRECTORIES:
        for path in sorted((root / directory).rglob("*")):
            relative = path.relative_to(root).as_posix()
            if path.is_file() and path.suffix in MODULE_SUFFIXES and not is_ignored(relative):
                yield path, relative


def check_folder_modules(files):
    modules = {}
    for path, relative in files:
        if path.name == "__init__.py":
            continue
        modules.setdefault(Path(relative).parent.as_posix(), set()).add(path.stem)
    violations = []
    for folder, stems in sorted(modules.items()):
        if len(stems) > MAX_MODULES_PER_FOLDER:
            message = f"{len(stems)} modules (max {MAX_MODULES_PER_FOLDER}); split into context subfolders"
            violations.append(Violation(folder, 1, "folder-modules", message))
    return violations


def main(root=None):
    if root is None:
        root = Path(__file__).resolve().parents[2]
    root = Path(root)
    files = list(scanned_files(root))
    violations = check_folder_modules(files)
    for path, relative in files:
        violations.extend(check_file(path, relative))
    for violation in violations:
        print(violation)
    return 1 if violations else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:2]))
