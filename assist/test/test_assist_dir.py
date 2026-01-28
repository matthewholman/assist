import os
import unittest
import tempfile

import assist


class TestAssistDir(unittest.TestCase):
    """
    Tests for ASSIST_DIR-based automatic planet ephemeris discovery.
    These exercise the C-side discovery logic in assist_ephem_init.
    """

    @staticmethod
    def _project_root() -> str:
        # This file lives at <root>/assist/test/test_assist_dir.py
        return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    def _with_assist_dir(self, path: str):
        """
        Context manager to temporarily set ASSIST_DIR.
        """
        class _Ctx:
            def __init__(self, new_path: str):
                self.new_path = new_path
                self.old = None

            def __enter__(self):
                self.old = os.environ.get("ASSIST_DIR")
                os.environ["ASSIST_DIR"] = self.new_path

            def __exit__(self, exc_type, exc, tb):
                if self.old is None:
                    os.environ.pop("ASSIST_DIR", None)
                else:
                    os.environ["ASSIST_DIR"] = self.old

        return _Ctx(path)

    def test_assist_dir_discovers_spk_de440(self):
        """
        With ASSIST_DIR pointing at the project root (which contains data/de440.bsp),
        calling Ephem() with no explicit path should auto-discover the SPK kernel
        and use the SPK-based planets provider.
        """
        root = self._project_root()
        de440_path = os.path.join(root, "data", "de440.bsp")
        if not os.path.exists(de440_path):
            self.skipTest("SPK de440.bsp not available in data/; skipping ASSIST_DIR SPK test.")

        with self._with_assist_dir(root):
            ephem = assist.Ephem()
            # planets_source == 0 corresponds to FILE_FORMAT_VALID_BSP
            self.assertEqual(ephem.planets_source, 0)
            # Sanity check: reuse basic Sun position from existing SPK tests
            p = ephem.get_particle(0, 0.0)
            self.assertAlmostEqual(p.x, -0.0071371791616079054)

    def test_assist_dir_discovers_ascii_bin_440_when_only_that_is_present(self):
        """
        Build a temporary ASSIST_DIR tree that only contains an ASCII-derived
        binary (.440/.441) planets file (and optional asteroid SPK). Ephem()
        should then auto-discover the ASCII-derived ephemeris and use the
        matching provider.
        """
        root = self._project_root()
        src_440 = os.path.join(root, "data", "linux_p1550p2650.440")
        if not os.path.exists(src_440):
            self.skipTest("linux_p1550p2650.440 not available in data/; skipping ASSIST_DIR .440 test.")

        src_ast = os.path.join(root, "data", "sb441-n16.bsp")

        with tempfile.TemporaryDirectory() as tmpdir:
            data_dir = os.path.join(tmpdir, "data")
            os.mkdir(data_dir)

            # Prefer symlinks to avoid copying large ephemeris files.
            def _safe_symlink(src: str, dst: str) -> bool:
                try:
                    os.symlink(src, dst)
                    return True
                except (OSError, NotImplementedError):
                    return False

            if not _safe_symlink(src_440, os.path.join(data_dir, "linux_p1550p2650.440")):
                self.skipTest("Symlinks not available; skipping ASSIST_DIR .440 test to avoid copying large data files.")

            if os.path.exists(src_ast):
                _safe_symlink(src_ast, os.path.join(data_dir, "sb441-n16.bsp"))

            with self._with_assist_dir(tmpdir):
                ephem = assist.Ephem()
                # planets_source == 2 corresponds to FILE_FORMAT_ASCII_BIN
                self.assertEqual(ephem.planets_source, 2)
                # Sanity check: Sun position matches existing JPL .440 tests
                p = ephem.get_particle(0, 0.0)
                self.assertAlmostEqual(p.x, -0.007137179161607906)


if __name__ == '__main__':
    unittest.main()







