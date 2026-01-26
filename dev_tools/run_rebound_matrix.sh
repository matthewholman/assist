#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

VERSIONS=(
  "4.4.11"
  "4.5.0"
)

REBOUND_SRC_BASE="${ROOT_DIR}/.rebound-src/github"

echo "ASSIST repo: ${ROOT_DIR}"
echo ""

python_rebound_lib_path() {
  # Prints the absolute path to the librebound shared library bundled with the
  # installed `rebound` Python package.
  "$1" - <<'PY'
import pathlib
import sysconfig
import rebound

site = pathlib.Path(rebound.__file__).resolve().parent.parent
candidates = []
suffix = sysconfig.get_config_var("EXT_SUFFIX") or ".so"
name = "librebound" + suffix
p = site / name
if p.exists():
    print(str(p))
else:
    # Fallback (should not be needed, but helps if naming changes)
    for pat in ("librebound*.so", "librebound*.dylib", "librebound*.dll"):
        candidates += sorted(site.glob(pat))
    print(str(candidates[0]) if candidates else "")
PY
}

python_ext_suffix() {
  "$1" - <<'PY'
import sysconfig
print(sysconfig.get_config_var("EXT_SUFFIX") or ".so")
PY
}

verify_python_env() {
  "$1" - <<'PY'
import rebound
import assist

print("python:", __import__("sys").executable)
print("rebound:", rebound.__version__)
print("assist:", assist.__version__, assist.__build__, assist.__githash__)
print("libassist:", assist.clibassist._name)
PY
}

run_python_tests_for_version() {
  local ver="$1"
  local py="${ROOT_DIR}/.venv-rebound-${ver}/bin/python"

  if [[ ! -x "${py}" ]]; then
    echo "ERROR: missing venv at ${py}"
    echo "Create it with: python3 -m venv .venv-rebound-${ver} && .venv-rebound-${ver}/bin/python -m pip install -U pip setuptools wheel && .venv-rebound-${ver}/bin/python -m pip install rebound==${ver} numpy"
    exit 1
  fi

  echo "=== Python tests (rebound==${ver}) ==="

  local suffix
  suffix="$(python_ext_suffix "${py}")"

  local librebound
  librebound="$(python_rebound_lib_path "${py}")"
  if [[ -z "${librebound}" ]]; then
    echo "ERROR: could not locate librebound shared library in venv for rebound==${ver}"
    exit 1
  fi

  # In this repo, `assist/__init__.py` loads `../libassist*.so` from the repo root,
  # and that binary uses `@rpath/librebound*.so`. For local testing, the simplest
  # way to ensure we're using the intended REBOUND is to place the matching
  # `librebound*` next to `libassist*` in the repo root.
  cp -f "${librebound}" "${ROOT_DIR}/$(basename "${librebound}")"
  cp -f "${librebound}" "${ROOT_DIR}/librebound.so"

  # Force a clean rebuild of the in-place Python extension for this venv.
  rm -rf "${ROOT_DIR}/build"
  rm -f "${ROOT_DIR}/libassist${suffix}"
  "${py}" "${ROOT_DIR}/setup.py" build_ext --inplace --force >/dev/null

  verify_python_env "${py}"
  echo ""

  # Existing runner
  make -C "${ROOT_DIR}" test_python PYTHON="${py}"
  echo ""
}

run_c_tests_for_version() {
  local ver="$1"
  local reb_dir="${ROOT_DIR}/.rebound-src/github/rebound-${ver}"

  echo "=== C tests (REBOUND source ${ver}) ==="
  if [[ ! -f "${reb_dir}/src/Makefile.defs" ]]; then
    echo "ERROR: missing REBOUND source build files at ${reb_dir}"
    echo "Expected to find ${reb_dir}/src/Makefile.defs"
    echo ""
    echo "Tip: download the tagged source from GitHub:"
    echo "  curl -L -o ${ROOT_DIR}/.rebound-src/rebound-${ver}.tar.gz https://github.com/hannorein/rebound/archive/refs/tags/${ver}.tar.gz"
    exit 1
  fi

  if [[ -f "${reb_dir}/version.txt" ]]; then
    echo -n "REB_DIR version.txt: "
    cat "${reb_dir}/version.txt"
  else
    echo "REB_DIR: ${reb_dir} (no version.txt found)"
  fi

  echo "REB_DIR: ${reb_dir}"

  # Ensure we don't accidentally reuse stale symlinks/executables when switching versions.
  make -C "${ROOT_DIR}/src" clean REB_DIR="${reb_dir}" >/dev/null 2>&1 || true
  rm -f "${ROOT_DIR}"/unit_tests/*/rebound
  rm -f "${ROOT_DIR}"/unit_tests/*/libassist.so
  rm -f "${ROOT_DIR}"/unit_tests/*/librebound.so
  rm -f "${ROOT_DIR}"/unit_tests/*/librebound*.so

  # Existing runner
  make -C "${ROOT_DIR}" test REB_DIR="${reb_dir}"
  echo ""
}

for ver in "${VERSIONS[@]}"; do
  run_python_tests_for_version "${ver}"
done

for ver in "${VERSIONS[@]}"; do
  run_c_tests_for_version "${ver}"
done

echo "All matrix runs completed successfully."


