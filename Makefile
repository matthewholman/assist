libassist: 
	$(MAKE) -C src 
	@ln -f -s src/libassist.so .
	@if [ "$(MAKELEVEL)" -eq "0" ]; then echo "To compile the example problems, go to a subdirectory of examples/ and execute make there."; fi
	
all: libassist

clean:
	$(MAKE) -C src clean
	$(MAKE) -C doc clean
	@-rm -f *.so
	@pip uninstall -y assist
	@python setup.py clean --all
	@rm -rf assist.*

.PHONY: doc
doc: 
	cd doc/doxygen && doxygen
	$(MAKE) -C doc html
		

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