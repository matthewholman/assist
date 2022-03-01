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
		
