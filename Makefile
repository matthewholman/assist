libassist: 
	$(MAKE) -C src 
	@ln -f -s src/libassist.so .
	@if [ "$(MAKELEVEL)" -eq "0" ]; then echo "To run tests, run 'make test' or go to a subdirectory of examples/ and execute make there."; fi
	
all: libassist

clean:
	$(MAKE) -C src clean
	@if [ -d doc ] && [ -f doc/Makefile ]; then $(MAKE) -C doc clean; fi
	@-rm -f *.so
	@pip uninstall -y assist 2>/dev/null || true
	@python setup.py clean --all 2>/dev/null || true
	@rm -rf assist.*

.PHONY: doc
doc: 
	@if [ -d doc/doxygen ]; then cd doc/doxygen && doxygen; fi
	@if [ -d doc ] && [ -f doc/Makefile ]; then $(MAKE) -C doc html; fi
		

# Iterate through each directory in the unit_tests directory and compile the test files
# then run them
test:
	@echo "Running tests ..."
	@set -e; \
	for dir in $(wildcard unit_tests/*); do \
		if [ -d $$dir ]; then \
			echo "Entering directory $$dir"; \
			if $(MAKE) -C $$dir; then \
				if [ -f "$$dir/rebound" ]; then \
					(cd $$dir && LD_LIBRARY_PATH=$(REB_DIR):. ./rebound); \
				else \
					echo "No rebound executable found in $$dir"; \
					exit 1; \
				fi \
			else \
				echo "Make failed in $$dir"; \
				exit 1; \
			fi \
		fi \
	done

# Python unit tests (requires a Python env with rebound + numpy installed).
#
# By default, prefer the repo venv if present; otherwise use python3 on PATH.
PYTHON ?= .venv/bin/python
PYTHON_ABS := $(abspath $(PYTHON))

test_python:
	@set -e; \
	if [ -x "$(PYTHON)" ]; then PY="$(PYTHON)"; else PY=python3; fi; \
	echo "Running Python tests with $$PY ..."; \
	$$PY -m unittest discover -s assist/test -p 'test_*.py'

# Same-libassist runner: build `src/libassist.so` against the selected Python's
# `rebound` wheel and force Python tests to load that exact binary via
# ASSIST_LIBASSIST_PATH. This guarantees Python + C tests exercise the same
# compiled ASSIST library.
NO_REB_DIR ?= /__no_rebound_source__

test_python_same_libassist:
	@set -e; \
	if [ -x "$(PYTHON)" ]; then PY="$(PYTHON_ABS)"; else PY=python3; fi; \
	echo "Building src/libassist.so against $$PY's rebound (python-wheel mode) ..."; \
	$(MAKE) -C src REB_DIR=$(NO_REB_DIR) PYTHON=$$PY; \
	echo "Running Python tests with $$PY (ASSIST_LIBASSIST_PATH=src/libassist.so) ..."; \
	ASSIST_LIBASSIST_PATH=src/libassist.so $$PY -m unittest discover -s assist/test -p 'test_*.py'

test_same_libassist: test_python_same_libassist
	@set -e; \
	if [ -x "$(PYTHON)" ]; then PY="$(PYTHON_ABS)"; else PY=python3; fi; \
	echo "Running C tests (make test) against $$PY's rebound (python-wheel mode) ..."; \
	$(MAKE) test REB_DIR=$(NO_REB_DIR) PYTHON=$$PY

# Focused position/trajectory tolerance tests (SPK + ASCII-derived .440/.441 paths).
test_positions_python:
	@set -e; \
	if [ -x "$(PYTHON)" ]; then PY="$(PYTHON)"; else PY=python3; fi; \
	echo "Running Python SPK position-tolerance tests with $$PY ..."; \
	$$PY -m unittest discover -s assist/test -p 'test_*_spk.py'; \
	echo "Running Python ASCII(.440/.441) position-tolerance tests with $$PY ..."; \
	$$PY -m unittest discover -s assist/test -p 'test_*_ascii.py'; \
	echo "Running Python SPK vs ASCII(.440/.441) drift test with $$PY ..."; \
	$$PY -m unittest assist.test.test_apophis_drift_spk_vs_ascii

test_positions_c:
	@echo "Running C position-tolerance tests (SPK + ASCII(.440/.441)) ..."
	@set -e; \
	for dir in \
		unit_tests/holman_spk unit_tests/holman_ascii \
		unit_tests/holman_reverse_spk unit_tests/holman_reverse_ascii \
		unit_tests/apophis_spk unit_tests/apophis_ascii \
		unit_tests/apophis_drift \
		unit_tests/5303_Ceres_spk unit_tests/5303_Ceres_ascii \
		unit_tests/interpolation_spk unit_tests/interpolation_ascii \
		unit_tests/convert_simulation_spk unit_tests/convert_simulation_ascii \
		unit_tests/variational_spk unit_tests/variational_ascii \
		unit_tests/roundtrip_spk unit_tests/roundtrip_ascii \
		unit_tests/roundtrip_adaptive_spk unit_tests/roundtrip_adaptive_ascii \
	; do \
		echo "Entering directory $$dir"; \
		$(MAKE) -C $$dir; \
		(cd $$dir && ./rebound); \
	done

test_positions: test_positions_python test_positions_c