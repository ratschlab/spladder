#/bin/bash

#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2013 Gunnar Raetsch, Regina Bohnert. Philipp Drewe, Andre Kahles
# Copyright (C) 2009-2013 Max Planck Society, Sloan-Kettering Institute
#

set -e 
echo "To get more information please run: configure -h" 
if [ -z "$1" -o "$1" == "-h" ];
then
		  echo Usage: "$0 -h   (this help)"
		  echo Usage: "$0 -i   (interactive configuration)"
		  echo Usage: "$0 [-a] (automatic configuration)"
  		  exit -1
fi

. ./bin/spladder_config.sh

echo "====================================== "
echo " SplAdder configure (version $SPLADDER_VERSION) " 
echo "====================================== "
echo

if [ "$1" == "-i" ]
then
  echo SplAdder base directory \(currently set to \"$SPLADDER_PATH\", suggest to set to \"`pwd`\", used if left empty\)
  read SPLADDER_PATH       
else
  echo "Running automatic configuration"
  SPLADDER_PATH=""
fi
if [ -z "$SPLADDER_PATH" ]
then
	SPLADDER_PATH=`pwd`
fi
echo '=>' Setting SplAdder base directory to \"$SPLADDER_PATH\"
echo

if [ "$1" == "-i" ];
then
  echo SAMTools directory \(currently set to \"$SPLADDER_SAMTOOLS_BIN_DIR\", system version used if left empty\)
  read SPLADDER_SAMTOOLS_BIN_DIR
else   
  echo Checking for samtools
  SPLADDER_SAMTOOLS_BIN_DIR=""
fi
if [ "$SPLADDER_SAMTOOLS_BIN_DIR" == "" ];
then
	if [ "$(which samtools)" != "" ] ;
	then
		SPLADDER_SAMTOOLS_BIN_DIR=$(which samtools)
	    if [ -f $(dirname $(which samtools))/../include/bam/sam.h ]
	    then
			SPLADDER_SAMTOOLS_INCLUDE_DIR=$(dirname $(which samtools))/../include/bam/
			echo "Include found: $SPLADDER_SAMTOOLS_INCLUDE_DIR"
	    elif [ -f $(dirname $(which samtools))/sam.h ]
	    then 
	               SPLADDER_SAMTOOLS_INCLUDE_DIR=$(dirname $(which samtools))/
                       echo "Include found: $SPLADDER_SAMTOOLS_INCLUDE_DIR"
	    else
			echo "ERROR: Include sam.h include not found"
			exit -1 ;
	    fi
	    if [ -f $(dirname $(which samtools))/../lib/libbam.a ]
	    then
			SPLADDER_SAMTOOLS_LIB_DIR=$(dirname $(which samtools))/../lib/
			echo "Library found: $SPLADDER_SAMTOOLS_LIB_DIR"
	    elif [ -f $(dirname $(which samtools))/libbam.a ]
	    then
		       SPLADDER_SAMTOOLS_LIB_DIR=$(dirname $(which samtools))/
                       echo "Library found: $SPLADDER_SAMTOOLS_LIB_DIR"
	    else
			echo "ERROR: Library libbam.a not found"
			exit -1 ;
	    fi
	else
	    echo SAMTools libraries not found
	    echo please run interactive mode: ./configure -i
	    exit -1 ;
	fi
else
	if [ ! -f $SPLADDER_SAMTOOLS_BIN_DIR ];
	then
		echo "ERROR: Binary $SPLADDER_SAMTOOLS_BIN_DIR not found"
                echo please run interactive mode: ./configure -i
		exit -1 ;
	fi

	echo SAMTools Include directory \(currently set to \"$SPLADDER_SAMTOOLS_INCLUDE_DIR\"\)
	read SPLADDER_SAMTOOLS_INCLUDE_DIR
	if [ ! -f $SPLADDER_SAMTOOLS_INCLUDE_DIR/sam.h ]
	then
		echo "ERROR: Include $SPLADDER_SAMTOOLS_INCLUDE_DIR/sam.h include not found"
                echo please run interactive mode: ./configure -i
		exit -1 ;
	fi
	
	echo SAMTools library directory \(currently set to \"$SPLADDER_SAMTOOLS_LIB_DIR\"\)
	read SPLADDER_SAMTOOLS_LIB_DIR
	if [ ! -f $SPLADDER_SAMTOOLS_LIB_DIR/libbam.a ]
	then
		echo "ERROR: Library $SPLADDER_SAMTOOLS_LIB_DIR/libbam.a include not found"
                echo please run interactive mode: ./configure -i
		exit -1 ;
	fi
fi
echo '=>' Setting samtools directory to \"$SPLADDER_SAMTOOLS_BIN_DIR\"
echo

if [ "$1" == "-i" ];
then
  echo Path to the python binary \(currently set to \"$SPLADDER_PYTHON_PATH\", system version used, if left empty\)
  read SPLADDER_PYTHON_PATH
else
  echo Checking for python and Scipy
  SPLADDER_PYTHON_PATH="" 
fi
if [ "$SPLADDER_PYTHON_PATH" == "" ];
then
	python_found=
    scipy_found=
	for i in python python2.7 python2.6 python2.5 python2.4;
    do	
      SPLADDER_PYTHON_PATH=`which $i`
  	  python_found=$i
  	  if [ "$SPLADDER_PYTHON_PATH" != "" ];
 	  then
	    scipy=`echo import scipy | $SPLADDER_PYTHON_PATH 2>&1 | grep  -e ImportError|wc -l`
		if [ $scipy == "0" ];
        then
		  scipy_found=$SPLADDER_PYTHON_PATH
          break
		fi
	  fi
    done 
    if [ "$python_found" == "" ];
    then 
      echo "EROR: Python not found"
      echo please run interactive mode: ./configure -i
      exit -1 
    fi
    if [ "$scipy_found" == "" ];
    then 
      echo "EROR: Scipy not found (for $python_found)"
      echo please run interactive mode: ./configure -i
      exit -1 
    fi
fi
echo '=>' Setting Python path to \"$SPLADDER_PYTHON_PATH\"
echo

if [ "$1" == "-i" ];
then
	echo "Please enter prefered interpreter ( octave or matlab )"
	read SPLADDER_INTERPRETER
elif [ "$1" == "-o" ] ; then
     SPLADDER_INTERPRETER="octave"
elif [ "$1" == "-m" ] ; then
     SPLADDER_INTERPRETER="matlab"
else
	SPLADDER_INTERPRETER="octave"
fi

if [ "$SPLADDER_INTERPRETER" == 'octave' ];
then
  if [ "$1" == "-i" ];
  then
	echo Please enter the full path to octave \(currently set to \"$SPLADDER_OCTAVE_BIN_PATH\", system version used, if left empty\)
	read SPLADDER_OCTAVE_BIN_PATH
  else
    SPLADDER_OCTAVE_BIN_PATH="" 
	echo checking for octave 
  fi
	if [ "$SPLADDER_OCTAVE_BIN_PATH" == "" ];
	then
	    SPLADDER_OCTAVE_BIN_PATH=`which octave` 
		if [ "$SPLADDER_OCTAVE_BIN_PATH" == "" ];
		then
			echo octave not found
                        echo please run interactive mode: ./configure -i
			exit -1
		fi
	fi
	echo '=>' Setting octave\'s path to \"$SPLADDER_OCTAVE_BIN_PATH\"
  if [ "$1" == "-i" ];
  then
	echo Please enter the full path to mkoctfile \(currently set to \"$SPLADDER_OCTAVE_MKOCT\", system version used, if left empty\)
	read SPLADDER_OCTAVE_MKOCT
  else 
    SPLADDER_OCTAVE_MKOCT="" 
  fi
	if [ "$SPLADDER_OCTAVE_MKOCT" == "" ];
	then
	    SPLADDER_OCTAVE_MKOCT=`which mkoctfile` 
		if [ "$SPLADDER_OCTAVE_MKOCT" == "" ];
		then
			SPLADDER_OCTAVE_MKOCT=$(dirname $SPLADDER_OCTAVE_BIN_PATH)/mkoctfile
			if [ ! -f SPLADDER_OCTAVE_MKOCT ];
			then
				echo mkoctfile not found
				echo please run interactive mode: ./configure -i
				exit -1
			fi
		fi
	fi
	echo '=>' Setting octave\'s path to \"$SPLADDER_OCTAVE_MKOCT\"
	echo
fi

if [ "$SPLADDER_INTERPRETER" == 'matlab' ];
then
	echo Please enter the full path to matlab \(currently set to \"$SPLADDER_MATLAB_BIN_PATH\", system version used, if left empty\)
	read SPLADDER_MATLAB_BIN_PATH
	if [ "$SPLADDER_MATLAB_BIN_PATH" == "" ];
	then
		SPLADDER_MATLAB_BIN_PATH=`which matlab`
		if [ "$SPLADDER_MATLAB_BIN_PATH" == "" ];
		then
			echo matlab not found
			echo please run interactive mode: ./configure -i
			exit -1
		fi
	fi
	if [ ! -f $SPLADDER_MATLAB_BIN_PATH ];
	then
		echo matlab not found
		echo please run interactive mode: ./configure -i
		exit -1
	fi
	echo '=>' Setting matlab\'s path to \"$SPLADDER_MATLAB_BIN_PATH\"
	echo

	echo Please enter the full path to mex binary \(currently set to \"$SPLADDER_MATLAB_MEX_PATH\", system version used if left empty\)
	read SPLADDER_MATLAB_MEX_PATH
	if [ "$SPLADDER_MATLAB_MEX_PATH" == "" ];
	then
		SPLADDER_MATLAB_MEX_PATH=`which mex`
		if [ "$SPLADDER_MATLAB_MEX_PATH" == "" ];
		then
			echo mex not found
			echo please run interactive mode: ./configure -i
			exit -1
		fi
	fi
	if [ ! -f "$SPLADDER_MATLAB_MEX_PATH" ];
	then
		echo mex not found
		echo please run interactive mode: ./configure -i
		exit -1
	fi
	echo '=>' Setting mex\' path to \"$SPLADDER_MATLAB_MEX_PATH\"
	echo

	echo Please enter the full path to the matlab include directory \(currently set to \"$SPLADDER_MATLAB_INCLUDE_DIR\", system version used, if left empty\)
	read SPLADDER_MATLAB_INCLUDE_DIR
	if [ "$SPLADDER_MATLAB_INCLUDE_DIR" == "" ];
	then
		SPLADDER_MATLAB_INCLUDE_DIR=$(dirname $SPLADDER_MATLAB_BIN_PATH)/../extern/include
	fi
	if [ ! -d "$SPLADDER_MATLAB_INCLUDE_DIR" ];
	then
		echo matlab include dir not found
		echo please run interactive mode: ./configure -i
		exit -1
	fi
	echo '=>' Setting matlab\'s include directory to \"$SPLADDER_MATLAB_INCLUDE_DIR\"
	echo

	SPLADDER_OCTAVE_BIN_PATH=
fi

cp -p bin/rdiff_config.sh bin/rdiff_config.sh.bak

grep -v -e SPLADDER_OCTAVE_BIN_PATH -e SPLADDER_OCTAVE_MKOCT -e SPLADDER_MATLAB_BIN_PATH -e SPLADDER_MATLAB_MEX_PATH -e SPLADDER_MATLAB_INCLUDE_DIR \
    -e SPLADDER_PATH -e SPLADDER_SRC_PATH -e SPLADDER_BIN_PATH \
    -e SPLADDER_INTERPRETER bin/rdiff_config.sh.bak \
    -e SPLADDER_SAMTOOLS_BIN_DIR -e SPLADDER_SAMTOOLS_LIB_DIR -e SPLADDER_SAMTOOLS_INCLUDE_DIR -e SPLADDER_PYTHON_PATH -e SCIPY_PATH -e SPLADDER_VERSION > bin/rdiff_config.sh

echo Generating config file ... 

# appending the relevant lines to rdiff_config.sh
echo export SPLADDER_VERSION=$SPLADDER_VERSION >> bin/rdiff_config.sh
echo export SPLADDER_PATH=$SPLADDER_PATH >> bin/rdiff_config.sh
echo export SPLADDER_SRC_PATH=${SPLADDER_PATH}/src >> bin/rdiff_config.sh
echo export SPLADDER_BIN_PATH=${SPLADDER_PATH}/bin >> bin/rdiff_config.sh
echo export SPLADDER_INTERPRETER=$SPLADDER_INTERPRETER >> bin/rdiff_config.sh
echo export SPLADDER_MATLAB_BIN_PATH=$SPLADDER_MATLAB_BIN_PATH >> bin/rdiff_config.sh
echo export SPLADDER_MATLAB_MEX_PATH=$SPLADDER_MATLAB_MEX_PATH >> bin/rdiff_config.sh
echo export SPLADDER_MATLAB_INCLUDE_DIR=$SPLADDER_MATLAB_INCLUDE_DIR >> bin/rdiff_config.sh
echo export SPLADDER_OCTAVE_BIN_PATH=$SPLADDER_OCTAVE_BIN_PATH >> bin/rdiff_config.sh
echo export SPLADDER_OCTAVE_MKOCT=$SPLADDER_OCTAVE_MKOCT >> bin/rdiff_config.sh
echo export SPLADDER_SAMTOOLS_BIN_DIR=$SPLADDER_SAMTOOLS_BIN_DIR >> bin/rdiff_config.sh  
echo export SPLADDER_SAMTOOLS_LIB_DIR=$SPLADDER_SAMTOOLS_LIB_DIR >> bin/rdiff_config.sh  
echo export SPLADDER_SAMTOOLS_INCLUDE_DIR=$SPLADDER_SAMTOOLS_INCLUDE_DIR >> bin/rdiff_config.sh  
echo export SPLADDER_PYTHON_PATH=$SPLADDER_PYTHON_PATH >> bin/rdiff_config.sh

echo Done.
echo 

echo Please use \'make\' to compile the mex files before using SplAdder.
echo To test SplAdder use \'make example\' or \'make threeexamples\'.
echo
