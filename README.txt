        MTPSD - Multi-Taper Power Spectral Density estimator

Copyright (c) 2010 Antonio Sanchez <antonio@eigenspectrum.com>

Licence: GPLv3

Project web page: http://sourceforge.net/projects/mtpsd/

INSTALLATION
============

Dependencies:  FFTW3, LAPACK, (optional: GNU Octave)

Windows:  Adjust variables in the Makefile to point to your compiler/Octave installation.
          e.g.  CC=mingw32-g++-4.4.0-dw2.exe
                MINGW_PATH=/d/local/octave/mingw32
                LIB_PATH=/d/local/octave/lib

*nix:     May need to set compiler in Makefile
          e.g.  CC=g++

Call
    make

Copy 'lib/libmtpsd.a' and/or 'lib/libdpss.a' to your system's lib directory
Copy 'include/*.h' to your system's include directory
Copy 'bin/mtpsd.oct' and 'bin/dpss.oct' somewhere where Octave can find them
Copy 'bin/dpss' to your system's bin directory

See the mtpsd documentation for more options/details.