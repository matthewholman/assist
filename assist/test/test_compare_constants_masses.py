from __future__ import annotations

import os
import struct
import ctypes
import unittest

import assist


class TestCompareConstantsAndMasses(unittest.TestCase):
    """
    Compare extracted constants and masses between:
    - ASCII-derived binary (.440/.441) planets file
    - SPK planets kernel (.bsp) comments parser

    This validates that any remaining behavioral differences are not due to
    constants/masses extraction drift.
    """

    @staticmethod
    def _parse_ascii_constants(path: str) -> dict[str, float]:
        # Mirror `assist_ascii_init` header parsing enough to locate the constants table.
        JPL_N = 15
        JPL_NUT = 11
        JPL_TDB = 14

        with open(path, "rb") as f:
            f.seek(0x0A5C)
            beg, end, inc = struct.unpack("<ddd", f.read(24))
            num = struct.unpack("<i", f.read(4))[0]
            cau = struct.unpack("<d", f.read(8))[0]
            cem = struct.unpack("<d", f.read(8))[0]

            off = []
            ncf = []
            niv = []
            for _ in range(12):
                off.append(struct.unpack("<i", f.read(4))[0])
                ncf.append(struct.unpack("<i", f.read(4))[0])
                niv.append(struct.unpack("<i", f.read(4))[0])

            _ver = struct.unpack("<i", f.read(4))[0]
            off.append(struct.unpack("<i", f.read(4))[0])
            ncf.append(struct.unpack("<i", f.read(4))[0])
            niv.append(struct.unpack("<i", f.read(4))[0])

            f.seek(0x00FC)
            names = [f.read(6).decode("ascii") for _ in range(400)]
            f.seek(0x0B28)
            names.extend(f.read(6).decode("ascii") for _ in range(400, num))

            # columns 14/15, read as in C (we don't use values, but advances stream)
            for _ in range(2):
                off.append(struct.unpack("<i", f.read(4))[0])
                ncf.append(struct.unpack("<i", f.read(4))[0])
                niv.append(struct.unpack("<i", f.read(4))[0])

            ncm = [3] * JPL_N
            ncm[JPL_NUT] = 2
            ncm[JPL_TDB] = 1
            off = [x - 1 for x in off]

            rec = 16  # 2 doubles
            for p in range(JPL_N):
                rec += 8 * ncf[p] * niv[p] * ncm[p]

            f.seek(rec)
            con = list(struct.unpack("<" + "d" * num, f.read(8 * num)))

        m = {n.strip(): float(v) for n, v in zip(names, con)}
        # Keep header copies too (useful for debugging)
        m["cau"] = float(cau)
        m["cem"] = float(cem)
        _ = beg, end, inc  # unused, but validated by parse
        return m

    def test_de440_constants_and_masses_match_between_440_and_bsp(self) -> None:
        root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        ascii_bin_path = os.path.join(root, "data", "linux_p1550p2650.440")
        spk_path = os.path.join(root, "data", "de440.bsp")
        if not (os.path.exists(ascii_bin_path) and os.path.exists(spk_path)):
            self.skipTest("Required data files not present in data/.")

        ascii_bin = self._parse_ascii_constants(ascii_bin_path)

        # Call into libassist for SPK constant+mass extraction.
        clib = assist.clibassist

        class MassData(ctypes.Structure):
            _fields_ = [
                ("names", ctypes.POINTER(ctypes.c_char_p)),
                ("values", ctypes.POINTER(ctypes.c_double)),
                ("count", ctypes.c_size_t),
            ]

        class SpkConstantsAndMasses(ctypes.Structure):
            _fields_ = [
                ("AU", ctypes.c_double),
                ("EMRAT", ctypes.c_double),
                ("J2E", ctypes.c_double),
                ("J3E", ctypes.c_double),
                ("J4E", ctypes.c_double),
                ("J2SUN", ctypes.c_double),
                ("RE", ctypes.c_double),
                ("CLIGHT", ctypes.c_double),
                ("ASUN", ctypes.c_double),
                ("masses", MassData),
            ]

        clib.assist_load_spk_constants_and_masses.argtypes = [ctypes.c_char_p]
        clib.assist_load_spk_constants_and_masses.restype = SpkConstantsAndMasses
        clib.assist_free_spk_constants_and_masses.argtypes = [ctypes.POINTER(SpkConstantsAndMasses)]
        clib.assist_free_spk_constants_and_masses.restype = None

        res = clib.assist_load_spk_constants_and_masses(spk_path.encode("ascii"))
        try:
            # Core constants should match.
            self.assertEqual(ascii_bin["AU"], float(res.AU))
            self.assertEqual(ascii_bin["EMRAT"], float(res.EMRAT))
            self.assertEqual(ascii_bin["J2E"], float(res.J2E))
            self.assertEqual(ascii_bin["J3E"], float(res.J3E))
            self.assertEqual(ascii_bin["J4E"], float(res.J4E))
            self.assertEqual(ascii_bin["J2SUN"], float(res.J2SUN))
            self.assertEqual(ascii_bin["RE"], float(res.RE))
            self.assertEqual(ascii_bin["CLIGHT"], float(res.CLIGHT))
            self.assertEqual(ascii_bin["ASUN"], float(res.ASUN))

            masses = {
                res.masses.names[i].decode("ascii"): float(res.masses.values[i])
                for i in range(res.masses.count)
            }

            # Planet GMs should match.
            for k in ["GMS", "GM1", "GM2", "GMB", "GM4", "GM5", "GM6", "GM7", "GM8", "GM9"]:
                self.assertIn(k, ascii_bin)
                self.assertIn(k, masses)
                self.assertEqual(ascii_bin[k], masses[k])

            # Asteroid masses (MAxxxx) should match, too.
            ascii_bin_ma = {k: v for k, v in ascii_bin.items() if k.startswith("MA")}
            spk_ma = {k: v for k, v in masses.items() if k.startswith("MA")}
            self.assertEqual(set(ascii_bin_ma.keys()), set(spk_ma.keys()))
            for k in ascii_bin_ma.keys():
                self.assertEqual(ascii_bin_ma[k], spk_ma[k])
        finally:
            clib.assist_free_spk_constants_and_masses(ctypes.byref(res))


if __name__ == "__main__":  # pragma: no cover
    unittest.main()


