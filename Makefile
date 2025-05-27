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