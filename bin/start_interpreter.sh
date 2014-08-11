#/bin/bash

#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2011 Regina Bohnert, Gunnar Raetsch
# Copyright (C) 2009-2011 Max Planck Society
#

set -e

. `dirname $0`/spladder_config.sh

export MATLAB_RETURN_FILE=`mktemp -t spladder.XXXXXXXXXX.tmp` 


if [ "$SPLADDER_INTERPRETER" == 'octave' ];
then
	echo exit | ${SPLADDER_OCTAVE_BIN_PATH} -q --eval "global SHELL_INTERPRETER_INVOKE; SHELL_INTERPRETER_INVOKE=1; warning('off', 'Octave:shadowed-function'); warning('off', 'Octave:deprecated-function') ; addpath $SPLADDER_SRC_PATH;  $1('$2'); exit;" || (echo starting Octave failed; rm -f $MATLAB_RETURN_FILE; exit -1) ;
elif [ "$SPLADDER_INTERPRETER" == 'matlab' ];
then
	echo exit | ${SPLADDER_MATLAB_BIN_PATH} -nodisplay -r "global SHELL_INTERPRETER_INVOKE; SHELL_INTERPRETER_INVOKE=1; addpath $SPLADDER_SRC_PATH;  $1('$2'); exit;" || (echo starting Matlab failed; rm -f $MATLAB_RETURN_FILE; exit -1) ;
	${SPLADDER_MATLAB_BIN_PATH} -nodisplay -r "addpath $SPLADDER_SRC_PATH;  $1('$2'); exit;" || (echo starting Matlab failed; rm -f $MATLAB_RETURN_FILE; exit -1) ;
fi

test -f $MATLAB_RETURN_FILE || exit 0
ret=`cat $MATLAB_RETURN_FILE` ;
rm -f $MATLAB_RETURN_FILE
exit $ret


