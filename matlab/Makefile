include bin/spladder_config.sh

#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2014 Gunnar Raetsch, Regina Bohnert, Philipp Drewe, Andre Kahles
# Copyright (C) 2009-2014 Max Planck Society, Sloan-Kettering Institute
#

all:	mexfiles

mexfiles:
ifeq ($(SPLADDER_INTERPRETER),octave)
	echo Entering ./mex
	cd mex ; make octave
else
	echo Entering ./mex
	cd mex ; make matlab
endif

clean:	
	echo Entering ./mex
	cd mex ; make clean

example:
	echo Entering ./examples
	cd examples ; make example

