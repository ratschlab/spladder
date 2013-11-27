#include <stdio.h>
#include <mex.h>

#ifndef __MEX_INPUT_h__
#define __MEX_INPUT_h__
	char *get_string(const mxArray *prhs);
	bool get_bool(const mxArray *prhs);
	int get_int(const mxArray *prhs);
#endif
