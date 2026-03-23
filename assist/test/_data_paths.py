from __future__ import annotations

import os


def project_root() -> str:
    # This file lives at <root>/assist/test/_data_paths.py
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def assist_dir() -> str:
    """
    Return the ASSIST_DIR root for tests.

    CI sets ASSIST_DIR explicitly. Locally, fall back to the repo root so tests
    can be run from a source checkout without extra environment setup.
    """
    return os.environ.get("ASSIST_DIR") or project_root()


def data_path(name: str) -> str:
    return os.path.join(assist_dir(), "data", name)

