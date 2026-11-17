"""The shared library exports the C interface and nothing else (python/scripts/check_native_deps.py)."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python" / "scripts"))

import check_native_deps  # noqa: E402

READELF = """Symbol table '.dynsym' contains 5 entries:
   Num:    Value          Size Type    Bind   Vis      Ndx Name
     0: 0000000000000000     0 NOTYPE  LOCAL  DEFAULT  UND
     1: 0000000000000000     0 FUNC    GLOBAL DEFAULT  UND memcpy@GLIBC_2.14 (3)
     2: 0000000000001000    40 FUNC    GLOBAL DEFAULT   12 po_solve
     3: 0000000000002000    16 FUNC    GLOBAL DEFAULT   12 _ZNSt13runtime_errorC1EPKc
     4: 0000000000003000     8 OBJECT  WEAK   DEFAULT   20 _ZTVSt9exception@@GLIBCXX_3.4
""".splitlines()


class ForeignExports(unittest.TestCase):
    def test_linux_lists_defined_symbols_that_are_not_the_interface(self) -> None:
        self.assertEqual(
            check_native_deps.foreign_exports("Linux", READELF),
            ["_ZNSt13runtime_errorC1EPKc", "_ZTVSt9exception"],
        )

    def test_macos_names_carry_an_underscore(self) -> None:
        nm = ["0000000000001000 T _po_solve", "0000000000002000 T __ZNSt13runtime_errorC1EPKc"]
        self.assertEqual(check_native_deps.foreign_exports("Darwin", nm), ["_ZNSt13runtime_errorC1EPKc"])

    def test_linker_boilerplate_is_not_an_export_of_ours(self) -> None:
        # musl and aarch64 toolchains leave the C runtime's section symbols in the dynamic table.
        readelf = READELF[:3] + [
            "     7: 0000000000000600     0 FUNC    GLOBAL DEFAULT   10 _init",
            "     8: 0000000000000700     0 FUNC    GLOBAL DEFAULT   11 _fini",
            "     9: 0000000000000800     0 FUNC    GLOBAL DEFAULT   10 .init",
            "    10: 0000000000000900     0 FUNC    GLOBAL DEFAULT   11 .fini",
        ]
        self.assertEqual(check_native_deps.foreign_exports("Linux", readelf), [])

    def test_section_symbols_are_not_exports(self) -> None:
        # aarch64 lists a section symbol such as .data as GLOBAL; it names no function or object.
        readelf = READELF[:3] + ["    11: 0000000000004000     0 SECTION GLOBAL DEFAULT   22 .data"]
        self.assertEqual(check_native_deps.foreign_exports("Linux", readelf), [])

    def test_only_the_interface_passes(self) -> None:
        only_interface = [line for line in READELF if "po_solve" in line or "memcpy" in line or "Num:" in line]
        self.assertEqual(check_native_deps.foreign_exports("Linux", only_interface), [])


if __name__ == "__main__":
    unittest.main()
