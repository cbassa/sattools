# Makefile

prefix = /usr/local
bindir = $(prefix)/bin

all:
	mkdir -p ./bin
	cp scripts/tleupdate bin/
	$(MAKE) -C src

clean:
	rm -f *.o
	rm -f *~
	rm -f src/*.o
	rm -f src/*~
	rm -rf bin
