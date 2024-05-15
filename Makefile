all:
	cd sp8cc && $(MAKE)

install:
	pip install .

uninstall:
	pip uninstall sp8py

check:
	(cd test && $(MAKE) check)

clean:
	rm sp8cc/sp8cc.so
