/* written by Jonas Behr, Regina Bohnert and Gunnar Raetsch, FML Tuebingen, Germany, 2010 */

#include <stdio.h>
#include <string.h>
#include <signal.h>
#include <mex.h>
#include <algorithm>
#include <vector>
	using std::vector;
#include "get_reads_direct.h"
#include "mex_input.h"
#include "read.h"

#define MAXLINE 10000

/*
 * input: 
 * 1 bam file
 * 2 chromosome
 * 3 region start (1-based index)
 * 4 region end (1-based index)
 * 5 strand (either '+' or '-' or '0')
 * [6] collapse flag: if true the reads are collapsed to a coverage track
 * [7] subsample percentage: percentage of reads to be subsampled (in per mill)
 * [8] intron length filter
 * [9] exon length filter
 * [10] mismatch filter
 * [11] bool: use mapped reads for coverage
 * [12] bool: use spliced reads for coverage
 * [13] return maxminlen
 * [14] return pair coverage and pair index list
 * [15] only_clipped 
 * [15] switch of pair filter
 * [16] pair flag filter (0 no flag info, 1 right flag info, 2 left flag info, 3 both flag info)
 *
 * output: 
 * 1 coverage
 * [2] intron cell array
 * [3] pair coverage
 * [4] pair list
 *
 * example call: 
 * [cov introns] = get_reads('polyA_left_I+_el15_mm1_spliced.bam', 'I', 10000, 12000, '-', 1, 30);
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
	if (nrhs<5 || nrhs>17 || (nlhs<1 || nlhs>4)) {
		fprintf(stderr, "usage: [x [introns] [pair]] = get_reads(fname, chr, start, end, strand, [collapse], [subsample], [max intron length], [min exonlength], [max mismatches], [mapped], [spliced], [maxminlen], [pair], [only clipped], [all pairs], [pair flag filter]);\n");
		return; 
	}
	
	/* obligatory arguments
	 * **********************/
	char *fname = get_string(prhs[0]);
	//fprintf(stdout, "arg1: %s\n", fname);
	char *chr = get_string(prhs[1]);
	//fprintf(stdout, "arg2: %s\n", chr);
	int from_pos = get_int(prhs[2]);
	//fprintf(stdout, "arg3: %d\n", from_pos);
	int to_pos = get_int(prhs[3]);
	//fprintf(stdout, "arg4: %d\n", to_pos);
	char *strand = get_string(prhs[4]);
	//fprintf(stdout, "arg5: %s\n", strand);

	if (from_pos>to_pos)
		 mexErrMsgTxt("Start (arg 3) must be <= end (arg 4)\n");

	if (strand[0]!='+' && strand[0]!='-' && strand[0]!='0') 
		mexErrMsgTxt("Unknown strand (arg 5): either + or - or 0");

	/* optional arguments
	 * ******************/	
	int collapse = 0;
	if (nrhs>=6)
		collapse = get_int(prhs[5]);
	
	int subsample = 1000;	
	if (nrhs>=7)
		subsample = get_int(prhs[6]);
		
	int intron_len_filter = 1e9;
	if (nrhs>=8)
		intron_len_filter = get_int(prhs[7]);

	int exon_len_filter = -1;
	if (nrhs>=9)
		exon_len_filter = get_int(prhs[8]);

	int filter_mismatch = 1e9;
	if (nrhs>=10)
		filter_mismatch = get_int(prhs[9]);

	int mapped = 1;
	if (nrhs>=11)
		mapped = get_int(prhs[10]);

	int spliced = 1;
	if (nrhs>=12)
		spliced = get_int(prhs[11]);

	int maxminlen = 0;
	if (nrhs>=13)
		maxminlen = get_int(prhs[12]);

	int pair_cov = 0;
	if (nrhs>=14)
		pair_cov = get_int(prhs[13]);

	int only_clipped = 0;
	if (nrhs>=15)
		only_clipped = get_int(prhs[14]);

	int no_pair_filter = 0;
	if (nrhs>=16)
		no_pair_filter = get_int(prhs[15]);

    int pair_flag_filter = -1; 
    if (nrhs>=17)
        pair_flag_filter = get_int(prhs[16]);


	/* call function to get reads
	 * **************************/	
	char region[MAXLINE];
	sprintf(region, "%s:%i-%i", chr, from_pos, to_pos);

	vector<CRead*> all_reads;
	
	get_reads_from_bam(fname, region, &all_reads, strand[0], subsample);

	for (int i=0; i<all_reads.size(); i++) {
		all_reads[i]->strip_leftright_tag() ;
	}
	
	/* filter reads
	 * **************/	
	int left = 0;
	int right = 0;
	
	vector<CRead*> reads;
	for (int i=0; i<all_reads.size(); i++) {
		if (all_reads[i]->left)
			left++;
		if (all_reads[i]->right)
			right++;
		//if (all_reads[i]->max_intron_len()<intron_len_filter && all_reads[i]->min_exon_len()>exon_len_filter && all_reads[i]->get_mismatches()<=filter_mismatch && all_reads[i]->multiple_alignment_index==0)
		if (all_reads[i]->max_intron_len()<intron_len_filter && all_reads[i]->min_exon_len()>exon_len_filter && all_reads[i]->get_mismatches()<=filter_mismatch && (only_clipped==0 || all_reads[i]->is_clipped) && (pair_flag_filter < 0 || (((pair_flag_filter & 2) > 0 && all_reads[i]->left) || ((pair_flag_filter & 1) > 0 && all_reads[i]->right) || (pair_flag_filter == 0 && !(all_reads[i]->right || all_reads[i]->left)) )))
			reads.push_back(all_reads[i]);
	}

 
	/* prepare output
	 * **************/	
	int num_rows = reads.size();
	int num_pos = to_pos-from_pos+1;
 
	if (pair_cov==1 && nlhs>=3) {
		// sort reads by read_id
		//printf("\n\nleft:%i right:%i \n\n", left, right);
		//printf("\nreads[0]->read_id: %s\n", reads[0]->read_id);
		sort(reads.begin(), reads.end(), CRead::compare_by_read_id);
	}
	
	// read coverages collapsed 
	if (collapse) {
		plhs[0] = mxCreateNumericMatrix(1, num_pos, mxUINT32_CLASS, mxREAL);
		uint32_t *mask_ret = (uint32_t*) mxGetData(plhs[0]);
		if (num_pos>0 && mask_ret==NULL)
			mexErrMsgTxt("Error allocating memory\n");
		if (mapped && spliced) {
			for (int i=0; i<reads.size(); i++) {
				reads[i]->get_coverage(from_pos, to_pos, mask_ret);
			}
		} else {
			for (int i=0; i<reads.size(); i++) {
				ssize_t num_exons = reads[i]->block_starts.size();
				if ((num_exons==1 && mapped) || (num_exons>1 && spliced))
					reads[i]->get_coverage(from_pos, to_pos, mask_ret);
			}
		}
	}
	// reads not collapsed
	else {
		uint32_t nzmax = 0; // maximal number of nonzero elements 
		int len = to_pos-from_pos+1;
		for (uint i=0; i<reads.size(); i++) 
		{
			ssize_t num_exons = reads[i]->block_starts.size();
			if (!((mapped && spliced) || (num_exons==1 && mapped) || (num_exons>1 && spliced)))
			{
				continue;
			}
			for (uint n = 0; n < reads[i]->block_starts.size(); n++) 
			{
				uint32_t from, to;
				if (reads[i]->block_starts[n]+reads[i]->start_pos-from_pos >= 0)
					from = reads[i]->block_starts[n]+reads[i]->start_pos-from_pos;
				else
					from = 0;
				if (reads[i]->block_starts[n]+reads[i]->start_pos-from_pos+reads[i]->block_lengths[n] >= 0)
					to = reads[i]->block_starts[n]+reads[i]->start_pos-from_pos+reads[i]->block_lengths[n];
				else
					to = 0;
				for (int bp=from; bp<to&bp<len; bp++)
				{
					nzmax++;
				}
			}
		}
		// 1st row: row indices
		// 2nd row: column indices
		plhs[0] = mxCreateDoubleMatrix(2, nzmax, mxREAL);
		double *mask_ret = (double*) mxGetData(plhs[0]);
		if (nzmax>0 && mask_ret==NULL)
			mexErrMsgTxt("Error allocating memory\n");
		uint32_t mask_ret_c = 0; // counter
		for (uint i=0; i<reads.size(); i++) 
		{
			ssize_t num_exons = reads[i]->block_starts.size();
			if (!((mapped && spliced) || (num_exons==1 && mapped) || (num_exons>1 && spliced)))
			{
				continue;
			}
			reads[i]->get_reads_sparse(from_pos, to_pos, mask_ret, mask_ret_c, i);
		}
		if (mask_ret_c!=2*nzmax)
			mexErrMsgTxt("Error filling index arrays for sparse matrix\n");
	}
	// introns
	if (maxminlen==0 && nlhs>=2) {
			vector<int> intron_list;
			for (int i=0; i<reads.size(); i++) {
				reads[i]->get_introns(&intron_list);
			}
			
			plhs[1] = mxCreateNumericMatrix(2, intron_list.size()/2, mxUINT32_CLASS, mxREAL);
			uint32_t *p_intron_list = (uint32_t*) mxGetData(plhs[1]);
			for (int p = 0; p<intron_list.size(); p++) {
				p_intron_list[p] = intron_list[p];
			}
			intron_list.clear();	
		} else if (nlhs>=2) {
			vector<uint32_t> intron_starts;
			vector<uint32_t> intron_ends;
			vector<uint32_t> block_len1;
			vector<uint32_t> block_len2;
			for (int i=0; i<reads.size(); i++) {
				reads[i]->get_introns(&intron_starts, &intron_ends, &block_len1, &block_len2);
			}
			
			plhs[1] = mxCreateNumericMatrix(4, intron_starts.size(), mxINT32_CLASS, mxREAL);
			uint32_t *p_intron_list = (uint32_t*) mxGetData(plhs[1]);
			for (int p = 0; p<intron_starts.size(); p++) {	
				p_intron_list[4*p] = intron_starts[p];
				p_intron_list[(4*p)+1] = intron_ends[p];
				p_intron_list[(4*p)+2] = block_len1[p];
				p_intron_list[(4*p)+3] = block_len2[p];
			}	 
			intron_starts.clear() ; 
			intron_ends.clear() ;
			block_len1.clear() ;
			block_len2.clear() ;
	}
	if (pair_cov==1 && nlhs>=3) {
		plhs[2] = mxCreateNumericMatrix(1, num_pos, mxUINT32_CLASS, mxREAL);
		uint32_t *p_pair_map = (uint32_t*) mxGetData(plhs[2]);
		if (num_pos>0 && p_pair_map==NULL)
			mexErrMsgTxt("Error allocating memory\n");
		
		vector<int> pair_ids;
		
		int take_cnt = 0;
		int discard_cnt = 0; 
		int discard_strand_cnt = 0; 
		//printf("reads.size(): %i\n", reads.size());
		// find consecutive reads with the same id
		for (int i=0; i<((int) reads.size())-1; i++) 
		{
			//printf("reads[%i]->read_id: %s\n", i, reads[i]->read_id);
			int j = i+1;
			while(j<reads.size() && strcmp(reads[i]->read_id, reads[j]->read_id) == 0) 
			{
				if ((reads[i]->left && reads[j]->right) || (reads[j]->left && reads[i]->right) && (reads[i]->reverse != reads[j]->reverse)) 
				{
					if (reads[i]->get_last_position()==-1 || reads[j]->get_last_position()==-1)
						break;
					if (false)//(reads[i]->strand[0]=='0' && reads[j]->strand[0]=='0' )
					{ 
						// discard pairs without strand information
						discard_strand_cnt++;
					}
					else if (reads[i]->get_last_position()<reads[j]->start_pos && reads[j]->start_pos-reads[i]->get_last_position()<60000) 
					{
						int from = std::max(0, reads[i]->get_last_position()-from_pos);
						int to = std::min(num_pos-1, reads[j]->start_pos-from_pos);
						pair_ids.push_back(i);
						pair_ids.push_back(j);
						for (int k=from; k<to; k++)
							p_pair_map[k]++;
						take_cnt++;
					}
					else if (reads[i]->start_pos>reads[j]->get_last_position() && reads[j]->get_last_position()-reads[i]->start_pos<60000)
					{
						int from = std::max(0, reads[j]->get_last_position()-from_pos);
						int to = std::min(num_pos-1, reads[i]->start_pos-from_pos);
						pair_ids.push_back(i);
						pair_ids.push_back(j);
						for (int k=from; k<to; k++)
							p_pair_map[k]++;
						take_cnt++;
					}
					else
					{
						if (no_pair_filter>0 && reads[i]->start_pos<reads[j]->start_pos)
						{
							pair_ids.push_back(i);
							pair_ids.push_back(j);
							take_cnt++;
						}
						else if (no_pair_filter>0)
						{
							pair_ids.push_back(j);
							pair_ids.push_back(i);
							take_cnt++;
						}
						else
							discard_cnt++;
						//printf("istart:%i, ilast:%i  jstart:%i, jlast: %i\n", reads[i]->start_pos, reads[i]->get_last_position(), reads[j]->start_pos, reads[j]->get_last_position());
					}
				}
				else
					discard_cnt++;
				j++;
			}
		}
		//printf("take:%i, discard:%i discard_strand_cnt:%i\n", take_cnt, discard_cnt+discard_strand_cnt, discard_strand_cnt);
		
		if (nlhs>=4) {
			plhs[3] = mxCreateNumericMatrix(2, pair_ids.size()/2, mxUINT32_CLASS, mxREAL);
			uint32_t *pair_ids_ret = (uint32_t*) mxGetData(plhs[3]);
			if (pair_ids.size()>0 && pair_ids_ret==NULL)
				mexErrMsgTxt("Error allocating memory\n");
			for (int i=0; i<pair_ids.size(); i++) {
				pair_ids_ret[i] = pair_ids[i];
			}
		}
	}
	for (int i=0; i<all_reads.size(); i++)
		delete all_reads[i];
}

