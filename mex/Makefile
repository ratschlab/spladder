include ../bin/spladder_config.sh

all:	octave
octave:	get_reads.mex get_header.mex
matlab:	get_reads.mexa64 get_header.mexa64 

get_reads.mexa64:	get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp read.h
	${SPLADDER_MATLAB_MEX_PATH} -g get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp -I${SPLADDER_SAMTOOLS_INCLUDE_DIR} -L${SPLADDER_SAMTOOLS_LIB_DIR} -lbam -lz -lcurses $(SPLADDER_MATLAB_INCLUDE_DIR)

get_header.mexa64:	get_header.cpp mex_input.cpp
	${SPLADDER_MATLAB_MEX_PATH} -g get_header.cpp mex_input.cpp -I${SPLADDER_SAMTOOLS_INCLUE_DIR} -L${SPLADDER_SAMTOOLS_LIB_DIR} -lbam -lz -lcurses $(SPLADDER_MATLAB_INCLUDE_DIR)

get_reads.mex:	get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp
	${SPLADDER_OCTAVE_MKOCT} --mex get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp -I${SPLADDER_SAMTOOLS_INCLUDE_DIR} -L${SPLADDER_SAMTOOLS_LIB_DIR} -lbam -lz -lcurses

get_header.mex:	get_header.cpp mex_input.cpp
	${SPLADDER_OCTAVE_MKOCT} --mex get_header.cpp mex_input.cpp -I${SPLADDER_SAMTOOLS_INCLUDE_DIR} -L${SPLADDER_SAMTOOLS_LIB_DIR} -lbam -lz -lcurses

clean:	
	rm -f *.o *.mexa64 *.mex
