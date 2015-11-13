# Copyright 2010 Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
#
# Copying and distribution of this file, with or without
# modification, are permitted in any medium without royalty provided
# the copyright notice and this notice are preserved.  This file is
# offered as-is, without any warranty.

# source code, C language
CSRC	= io_png.c example/readpng.c example/mmms.c example/axpb.c 
# source code, all languages (only C here)
SRC	= $(CSRC)
# object files (partial compilation)
OBJ	= $(CSRC:.c=.o)
# binary executable programs
BIN	= $(filter example/%, $(CSRC:.c=))

# standard C compiler optimization options
COPT	= -O2 -funroll-loops -fomit-frame-pointer
# complete C compiler options
CFLAGS	= -ansi -pedantic -Wall -Wextra -Werror -pipe $(COPT)
# linker options
LDFLAGS	= -lpng -lm
# library build dependencies (none)
LIBDEPS =

# use local embedded libraries
ifdef LOCAL_LIBS
# library location
LIBDIR = ./libs/build/lib
INCDIR = ./libs/build/include
# libpng is required
LIBDEPS += libpng
# compile options to use the local libpng header
CFLAGS 	+= -I$(INCDIR) -DIO_PNG_LOCAL_LIBPNG
# link options to use the local libraries
LDFLAGS = $(LIBDIR)/libpng.a $(LIBDIR)/libz.a -lm
endif

# default target: the example programs
default: $(BIN)

# build the png library
.PHONY	: libpng
libpng	:
	$(MAKE) -C libs libpng

# partial C compilation xxx.c -> xxx.o
%.o	: %.c $(LIBDEPS)
	$(CC) $< -c $(CFLAGS) -I. -o $@

# final link of an example program
example/%	: example/%.o io_png.o $(LIBDEPS)
	$(CC) $< io_png.o $(LDFLAGS) -o $@

# cleanup
.PHONY	: clean distclean scrub
clean	:
	$(RM) $(OBJ)
	$(RM) *.flag
	$(MAKE) -C libs $@
distclean	: clean
	$(RM) $(BIN)
	$(RM) -r srcdoc
	$(MAKE) -C libs $@
scrub	: distclean
	$(MAKE) -C libs $@

################################################
# extra tasks

PROJECT	= io_png
DATE	= $(shell date -u +%Y%m%d)
RELEASE_TAG   = 0.$(DATE)

.PHONY	: srcdoc lint beautify test release
# source documentation
srcdoc	: $(SRC)
	doxygen doc/doxygen.conf
# code cleanup
beautify	: $(CSRC)
	for FILE in $^; do \
		expand $$FILE | sed 's/[ \t]*$$//' > $$FILE.$$$$ \
		&& indent -kr -i4 -l78 -nut -nce -sob -sc \
			$$FILE.$$$$ -o $$FILE \
		&& rm $$FILE.$$$$; \
	done
# static code analysis
lint	: $(CSRC)
	for FILE in $^; do \
		clang --analyze -ansi -I. $$FILE || exit 1; done;
	for FILE in $^; do \
		splint -ansi-lib -weak -I. $$FILE || exit 1; done;
	$(RM) *.plist
# code tests
test	: $(CSRC)
	sh -e test/run.sh && echo SUCCESS || ( echo ERROR; return 1)
# release tarball
release	: beautify lint test distclean
	git archive --format=tar --prefix=$(PROJECT)-$(RELEASE_TAG)/ HEAD \
	        | gzip > ../$(PROJECT)-$(RELEASE_TAG).tar.gz
