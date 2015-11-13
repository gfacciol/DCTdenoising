% io_png: simplified front-end to libpng

* overview
* license
* requirements
* usage
 * read
 * write
 * example
* compilation
 * local libraries
* todo
* copyright

# OVERVIEW

io_png.c contains high-level routines to read and write PNG images
using libpng. It only handles common use cases, and provides a
simplified interface.

# LICENSE

io_png.c is distributed under a GPL3+ or BSD license, at your
option. See the included copyright notice, conditions and disclaimer
for details.

# REQUIREMENTS

libpng is required, version >= 1.2.2. The source code and binaries
can be found at http://www.libpng.org/pub/png/libpng.html.

Note that libpng requires zlib for compression. The source code and
binaries can be found at http://www.zlib.net/.

io_png.c is ANSI C, and should compile on any system with any ANSI C
compiler.

# USAGE

Compile io_png.c with your program, and include io_png.h to get the
function declarations. You can use io_png.c with C or C++ code.

## READ

A PNG image is read into a single array. For multiple channel images,
the output array is de-interlaced and successively contains each
channel. For example, a color image with 30 rows and 40 columns is
read into a single array of 3600 cells, with:

* the first 1200 cells (30 x 40) containing the red channel
* the next 1200 cells containing the green channel
* the last 1200 cells containing the blue channel

In each channel, the image is stored row after row.

No image structure is needed, and the image size information is
collected via pointer parameters.

The two main front-end functions are, depending of the data type you want
to use:

* io_png_read_u8():
  read a PNG image as an unsigned char array
  - 16bit images are converted to 8bit with precision loss
  - 1, 2 and 4bit images are converted to 8bit without precision loss
* io_png_read_f32():
  read a PNG image as a float array
  - 16bit images are first converted to 8bit with precision loss
  - integer values are then converted to float

These functions have the same syntax:

    io_png_read_xxx(fname, &nx, &ny, &nc)
    - fname: file name; the standard input stream is used if fname is "-"
    - nx, ny, nc: variables to fill with the image size

Four secondary read functions can be used to force a color model:

* io_png_read_u8_rgb():
  convert gray images to RGB and strip the alpha channel
* io_png_read_u8_gray():
  convert RGB images to gray and strip the alpha channel
* io_png_read_f32_rgb():
  convert gray images to RGB and strip the alpha channel
* io_png_read_f32_gray():
  convert RGB images to gray and strip the alpha channel

These functions have the same syntax as the previous ones, except that
they don't need the &nc parameter.

## WRITE

A PNG image is written from a single array, with the same layout as
the one received from the read functions.

2 front-end functions are available; they all have the same syntax,
io_png_write_xxx(fname, data, nx, ny, nc):
    - fname: file name, the standard output stream is used if fname is "-"
    - data: image array
    - nx, ny, nc: image size

* io_png_write_u8(): write a PNG image from an unsigned char array
* io_png_write_f32(): write a PNG image from a float array
  - the array values are first rounded,
  - then limited to [0..255], with values lower than 0 set to 0 and
    values higher than 255 set to 255

## EXAMPLE

see example/readpng.c and example/mmms.c

# COMPILATION

You can compile the example codes located in the example folder using
the provided makefile, with the `make` command.

## LOCAL LIBRARIES

If libpng is not installed on your system, of if you prefer a local
static build, a mechanism is available to automatically download,
build and include libpng in your program:

1. run `make libpng`;
   this uses the makefiles from the `libs` folder to download and
   compile libpng and zlib, and builds the libraries into `libs/build`;
2. use the "-DIO_PNG_LOCAL_LIBPNG -I./libs/build/include" options to compile
   io_png.c;
3. add ./libs/build/lib/libpng.a ./libs/build/lib/libz.a to the list of
  files being linked into your program

This is automatically handled in the provided makefile for the example
code example/readpng.c; simply use the `make LOCAL_LIBS=1` command
instead of `make`.

# TODO

* handle 16bit data
* cmake support
* C++ wrappers (vector output, merged functions)
* implement proper gamma and RGBY conversion
* handle data as float before re-quantization


# COPYRIGHT

Copyright 2011 Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>

Copying and distribution of this README file, with or without
modification, are permitted in any medium without royalty provided
the copyright notice and this notice are preserved.  This file is
offered as-is, without any warranty.
