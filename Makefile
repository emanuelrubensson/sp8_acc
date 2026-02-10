all:
	cd sp8cc && $(MAKE)

install:
	pip install .

uninstall:
	pip uninstall sp8py

check:
	(cd test && $(MAKE) check)
	(cd testpy && $(MAKE) check)

clean:
	$(MAKE) -C test clean
	$(MAKE) -C sp8cc clean
	rm -rf sp8py/__pycache__/