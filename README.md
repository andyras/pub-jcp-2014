pub-jcp-2014
========

This repository contains data and figures used in "Theory of Ultrafast Photoinduced Electron Transfer from a Bulk Semiconductor to a Quantum Dot" as well as the source code used to generate the data and figures.

Pre-compiled data and figures are in the directory

    pre-compiled/

and the source code and Makefile are in the directory

    you-compile/

## Software requirements

To compile the dynamics code and create the data, you need:

0. A compatible OS.  This code has been tested on Linux and Mac OS X 10.9.2.
1. A C++ compiler, <a href="http://software.intel.com/en-us/intel-education-offerings">Intel</a> or <a href="http://gcc.gnu.org/">GCC</a> recommended
2. An implementation of the <a href="http://www.netlib.org/lapack/">LAPACK</a> linear algebra library, either Intel <a href="http://software.intel.com/en-us/intel-mkl">MKL</a> (included with Intel compiler suites) or another distribution
3. The <a href="http://computation.llnl.gov/casc/sundials/download/download.html">SUNDIALS</a> differential ordinary equation solver package
4. <a href="http://gnuplot.sourceforge.net/">gnuplot</a> (version 4.6.0+ recommended)

### Optional requirements

1.  To convert images to PDF from the default EPS, you need one of the `ps2pdf` (via <a href="http://www.ghostscript.com/">Ghostscript</a>) or `convert` (via <a href="http://www.imagemagick.org/">Imagemagick<a/>) utilities.
2.  To convert images to PNG from the default EPS, you need the <a href="http://www.imagemagick.org/">`convert`</a> utility.

Alternatively, you can modify the gnuplot `.plt` files to use a different terminal.

## Instructions

0.  Install the requirements as listed above.
1.  In the Makevars file, set the variables appropriately.
2.  To make everything, type `make` at the command line.

    To make just the code, type `make code`.

    To make just the data, type `make data`.

    To make just the plots, type `make plots`.

Note: figures which do not use data are not generated by this repository.
