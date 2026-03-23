import os
import shutil
import subprocess
import sys
import sysconfig
import tempfile
import unittest


class TestRelocatableLinking(unittest.TestCase):
    """
    Verify that the ASSIST C library can still find REBOUND after the
    install tree is moved to a different location on disk.

    This guards against regressions in the runtime link settings
    (e.g. rpath / install_name) used when building libassist.
    """

    def test_assist_and_rebound_work_after_relocation(self) -> None:
        try:
            import assist  # type: ignore
        except Exception as e:  # pragma: no cover - environment issue, not a test failure
            self.skipTest(f"assist import failed in test environment: {e!r}")
        try:
            import rebound  # type: ignore
        except Exception as e:  # pragma: no cover - environment issue, not a test failure
            self.skipTest(f"rebound import failed in test environment: {e!r}")

        # Locate the current assist package and its sibling libassist shared library.
        suffix = sysconfig.get_config_var("EXT_SUFFIX") or ".so"
        assist_pkg_dir = os.path.dirname(os.path.abspath(assist.__file__))
        assist_parent = os.path.dirname(assist_pkg_dir)
        libassist_src = os.path.join(assist_parent, "libassist" + suffix)
        if not os.path.exists(libassist_src):
            self.skipTest(
                f"Relocatable linking test requires a binary build with libassist{suffix} "
                f"adjacent to the assist package (looked in {assist_parent!r})."
            )

        # Locate the REBOUND Python package directory (which contains its shared libraries).
        rebound_pkg_dir = os.path.dirname(os.path.abspath(rebound.__file__))
        if not os.path.isdir(rebound_pkg_dir):
            self.skipTest(f"Could not locate rebound package directory at {rebound_pkg_dir!r}.")

        with tempfile.TemporaryDirectory() as tmpdir:
            # Construct a minimal "site-packages" style tree in a new location.
            site_dir = os.path.join(tmpdir, "site")
            os.mkdir(site_dir)

            # Copy the Python assist package, the libassist shared object, and the
            # entire rebound package directory (including its shared libraries)
            # into the new location.
            shutil.copytree(assist_pkg_dir, os.path.join(site_dir, "assist"))
            shutil.copy2(libassist_src, os.path.join(site_dir, os.path.basename(libassist_src)))
            shutil.copytree(rebound_pkg_dir, os.path.join(site_dir, "rebound"))

            # Also copy the librebound shared object which is usually in site-packages root
            # (sibling to rebound dir)
            rebound_parent = os.path.dirname(rebound_pkg_dir)
            for f in os.listdir(rebound_parent):
                if f.startswith("librebound") and f.endswith(suffix):
                    shutil.copy2(os.path.join(rebound_parent, f), site_dir)

            # Now spawn a *clean* Python process that only sees this relocated tree.
            code = r"""
import os
import sys

import assist

from assist import Ephem

e = Ephem()
p = e.get_particle(0, 0.0)  # Sun at t=0
print("OK", p.x)
"""

            env = os.environ.copy()
            env["PYTHONPATH"] = site_dir
            # Ensure we test the relocated install tree (and not a forced libassist
            # override from the parent process, which would defeat the purpose of
            # this test).
            env.pop("ASSIST_LIBASSIST_PATH", None)
            # Ensure ASSIST_DIR points to the original source so we can find data files
            # (or we could copy data/ to site_dir, but this is simpler)
            if "ASSIST_DIR" not in env:
                 # Fallback: assume the test is running from the repo root or near it
                 repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
                 env["ASSIST_DIR"] = repo_root

            # Avoid picking up any user/system site-packages that might hide
            # relocation problems.
            env["PYTHONNOUSERSITE"] = "1"

            # -S prevents automatic import of site, which further reduces the
            # chance that the original install tree masks relocation issues.
            cmd = [sys.executable, "-S", "-c", code]
            proc = subprocess.run(
                cmd,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            if proc.returncode != 0:
                self.fail(
                    "Relocated import or linking failed.\n"
                    f"STDOUT:\n{proc.stdout}\n"
                    f"STDERR:\n{proc.stderr}\n"
                )

            self.assertIn("OK", proc.stdout)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()


