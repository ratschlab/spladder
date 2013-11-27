#include <mex.h>
#include <sam.h>
#include "mex_input.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    char *fname = get_string(prhs[0]);

    // open alignment file
    samfile_t* infile = samopen(fname, "rb" , 0);
    bam_header_t* header = infile->header ;

    plhs[0] = mxCreateCellMatrix((mwSize)header->n_targets,2);
    size_t cnt = 0;
    for (size_t i = 0; i < header->n_targets; i++) {
        mxSetCell(plhs[0], cnt++, mxCreateString(header->target_name[i]));
    }
    for (size_t i = 0; i < header->n_targets; i++) {
        //mxArray *tmp = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
        //double *tmp2 = mxGetPr(tmp);
        //tmp2[0] = header->target_len[i];
        //mxSetCell(plhs[0], cnt++, tmp);
        mxSetCell(plhs[0], cnt++, mxCreateDoubleScalar((double) header->target_len[i]));
    }

    // close alignment file
    samclose(infile);

    return;
}
	
