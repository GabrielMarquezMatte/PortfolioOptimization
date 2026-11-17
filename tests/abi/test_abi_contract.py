"""Unit tests of scripts/abi_contract.py. Run: python -m unittest discover -s tests/abi"""

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import abi_contract  # noqa: E402


class ContractLines(unittest.TestCase):
    def test_multiline_declaration_is_one_line(self) -> None:
        header = "PO_API po_status po_f(const size_t* a,\n    size_t n);"
        self.assertEqual(abi_contract.contract_lines(header), ["function po_status po_f(const size_t*, size_t)"])

    def test_define_drops_trailing_comment(self) -> None:
        self.assertEqual(abi_contract.contract_lines("#define PO_X 1 /* note */"), ["define PO_X 1"])

    def test_pointer_to_const_pointer(self) -> None:
        header = "PO_API po_status po_g(const po_spec* const* specs, size_t n_specs);"
        self.assertEqual(
            abi_contract.contract_lines(header), ["function po_status po_g(const po_spec* const*, size_t)"]
        )

    def test_void_parameter_and_pointer_return(self) -> None:
        self.assertEqual(
            abi_contract.contract_lines("PO_API const char* po_last_error(void);"),
            ["function const char* po_last_error(void)"],
        )

    def test_typedefs_and_platform_macros(self) -> None:
        header = "#define PO_API __declspec(dllexport)\ntypedef int32_t po_status;\ntypedef struct po_context po_context;"
        self.assertEqual(
            abi_contract.contract_lines(header), ["typedef int32_t po_status", "typedef struct po_context po_context"]
        )

    def test_parameter_names_do_not_matter(self) -> None:
        a = abi_contract.contract_lines("PO_API void po_h(double* out, size_t count);")
        b = abi_contract.contract_lines("PO_API void po_h(double* values, size_t n);")
        self.assertEqual(a, b)


class Unparsed(unittest.TestCase):
    """A declaration the parser cannot read must be an error, never a line missing from the contract."""

    def test_function_pointer_parameter_is_refused(self) -> None:
        with self.assertRaisesRegex(ValueError, "po_on_progress"):
            abi_contract.contract_lines("PO_API void po_on_progress(void (*callback)(int), void* data);")

    def test_trailing_attribute_is_refused(self) -> None:
        with self.assertRaisesRegex(ValueError, "po_old"):
            abi_contract.contract_lines("PO_API void po_old(void) PO_DEPRECATED;")

    def test_name_without_prefix_is_refused(self) -> None:
        with self.assertRaisesRegex(ValueError, "other_name"):
            abi_contract.contract_lines("PO_API void other_name(void);")

    def test_anonymous_enum_and_struct_are_refused(self) -> None:
        with self.assertRaisesRegex(ValueError, "enum"):
            abi_contract.contract_lines("enum { A = 1 };")
        with self.assertRaisesRegex(ValueError, "struct"):
            abi_contract.contract_lines("typedef struct { int a; } po_pair;")

    def test_extern_c_braces_are_allowed(self) -> None:
        text = 'extern "C"\n{\nPO_API void po_f(void);\n}'
        self.assertEqual(abi_contract.contract_lines(text), ["function void po_f(void)"])


class Comparisons(unittest.TestCase):
    def test_removed_line_fails_base(self) -> None:
        self.assertEqual(abi_contract.removed_lines(["define A 1", "define B 2"], ["define A 1"]), ["define B 2"])

    def test_changed_line_fails_base(self) -> None:
        self.assertEqual(abi_contract.removed_lines(["define A 1"], ["define A 2"]), ["define A 1"])

    def test_added_line_passes_base(self) -> None:
        self.assertEqual(abi_contract.removed_lines(["define A 1"], ["define A 1", "define B 2"]), [])


class ExportList(unittest.TestCase):
    CONTRACT = ["define PO_X 1", "function int32_t po_a(void)", "function void po_b(po_x*)"]

    def test_def_lists_functions_only(self) -> None:
        text = abi_contract.def_file(self.CONTRACT, "x.dll")
        self.assertEqual(text.split(), ["LIBRARY", "x.dll", "EXPORTS", "po_a", "po_b"])

    def test_dll_exports_must_match(self) -> None:
        # The export table of `objdump -p`: every exported name, ours or not.
        table = (
            "[Ordinal/Name Pointer] Table -- Ordinal Base 1\n"
            "\t          Ordinal   Hint Name\n"
            "\t[   0] +base[   1]  0000 po_a\n"
            "\t[   1] +base[   2]  0001 po_b\n"
        )
        self.assertEqual(abi_contract.export_differences(self.CONTRACT, table), [])
        extra = table + "\t[   2] +base[   3]  0002 _ZNSt9exceptionD2Ev\n"
        self.assertEqual(
            abi_contract.export_differences(self.CONTRACT, extra),
            ["exported but not in the contract: _ZNSt9exceptionD2Ev"],
        )
        missing = "\t[   0] +base[   1]  0000 po_a\n"
        self.assertEqual(
            abi_contract.export_differences(self.CONTRACT, missing), ["in the contract but not exported: po_b"]
        )
        self.assertEqual(
            abi_contract.export_differences(self.CONTRACT, ""),
            ["in the contract but not exported: po_a", "in the contract but not exported: po_b"],
        )


class CommittedContract(unittest.TestCase):
    def test_header_matches_committed_contract(self) -> None:
        header = (ROOT / "include" / "portopt" / "portopt.h").read_text(encoding="utf-8")
        committed = (ROOT / "include" / "portopt" / "portopt.abi").read_text(encoding="utf-8")
        self.assertEqual(abi_contract.contract_lines(header), committed.splitlines())

    def test_stable_tags_start_at_one(self) -> None:
        tags = ["v1.0.0-rc.1", "v0.9.0", "v1.0.0", "v1.2.0", "v2.0.0-rc.1"]
        self.assertEqual(abi_contract.latest_stable(tags), "v1.2.0")
        self.assertIsNone(abi_contract.latest_stable(["v0.9.0", "v1.0.0-rc.1"]))

    def test_check_command_reports_a_difference(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            contract = Path(folder) / "portopt.abi"
            contract.write_text("define PO_ABI_VERSION 0\n", encoding="utf-8")
            header = Path(folder) / "portopt.h"
            header.write_text("#define PO_ABI_VERSION 1\n", encoding="utf-8")
            self.assertEqual(abi_contract.main(["check", "--header", str(header), "--contract", str(contract)]), 1)
            self.assertEqual(abi_contract.main(["write", "--header", str(header), "--contract", str(contract)]), 0)
            self.assertEqual(abi_contract.main(["check", "--header", str(header), "--contract", str(contract)]), 0)


if __name__ == "__main__":
    unittest.main()
