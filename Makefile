# C source code
CSRC	= io_png/io_png.c \
			mt19937ar.c
# C++ source code
CXXSRC	= DCTdenoising.cpp \
		DCT2D16x16.cpp \
		DCT2D.cpp \
		demo_DCTdenoising.cpp 

CXXSRC2 = img_diff_ipol.cpp img_mse_ipol.cpp

# all source code
SRC	= $(CSRC) $(CXXSRC) $(CXXSRC2) 

# C objects
COBJ	= $(CSRC:.c=.o)
# C++ objects
CXXOBJ	= $(CXXSRC:.cpp=.o)
CXXOBJ2	= $(CXXSRC2:.cpp=.o)
# all objects
OBJ	= $(COBJ) $(CXXOBJ)
OBJ2    = $(CXXOBJ2)
# binary target
BIN	= demo_DCTdenoising 
BIN2    = img_diff_ipol img_mse_ipol

default	: $(BIN) $(BIN2)

# C optimization flags
COPT	= -O3 -ftree-vectorize -funroll-loops

# C++ optimization flags
CXXOPT	= $(COPT)

# C compilation flags
CFLAGS	= $(COPT)  \
	-Wno-write-strings -ansi
# C++ compilation flags
CXXFLAGS	= $(CXXOPT)  \
	-Wno-write-strings -Wno-deprecated -ansi
# link flags
LDFLAGS	= -lpng -lm

# use local embedded libraries with `make LOCAL_LIBS=1`
ifdef LOCAL_LIBS
# library location
LIBDIR = io_png/libs/build/lib
INCDIR = io_png/libs/build/include
# libpng is required
LIBDEPS += libpng
# compile options to use the local libpng header
CFLAGS 	+= -I$(INCDIR) -D_LOCAL_LIBS
# link options to use the local libraries
LDFLAGS = $(LIBDIR)/libpng.a $(LIBDIR)/libz.a -lm
# io_png.o needs png.h
io_png/io_png.o	:  $(LIBDEPS)
endif

# use openMP with `make OMP=1`
ifdef OMP
CFLAGS	+= -fopenmp
CXXFLAGS	+= -fopenmp
LDFLAGS += -lgomp
else
CFLAGS	+= -Wno-unknown-pragmas
CXXFLAGS  += -Wno-unknown-pragmas
endif

# build the local png library
.PHONY	: libpng
libpng	:
	$(MAKE) -C io_png/libs libpng

# partial compilation of C source code
%.o: %.c %.h
	$(CC) -c -o $@  $< $(CFLAGS)
# partial compilation of C++ source code
%.o: %.cpp %.h
	$(CXX) -c -o $@  $< $(CXXFLAGS)

# link all the object code
$(BIN): $(OBJ) $(LIBDEPS)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

# link all the object code
$(BIN2) : % : %.o  io_png/io_png.o 
	$(CXX)  -o  $@ $^ $(LDFLAGS)


# housekeeping
.PHONY	: clean distclean
clean	:
	$(RM) $(OBJ)
	$(RM) $(OBJ2)
	$(MAKE) -C ./io_png/libs $@
distclean	: clean
	$(RM) $(BIN)
	$(RM) $(BIN2)
	$(MAKE) -C ./io_png/libs $@

