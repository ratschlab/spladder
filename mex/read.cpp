/* written by Jonas Behr, Regina Bohnert and Gunnar Raetsch, FML Tuebingen, Germany, 2010 */

#include "read.h"

CRead::CRead() {
	read_id = NULL;
	sam_line = NULL;
	start_pos = 0;
	matches = 0;
	mismatches = 0;
	multiple_alignment_index = 0;
	strand = NULL;
	left = false;
	right = false;
	reverse = false;
	is_clipped = false;
}

CRead::~CRead() {
	delete[] read_id;
	delete[] sam_line;
	delete[] strand;
}

/*
 * Augments 'coverage' array at the positions covered by the read in the queried interval.
 */
void CRead::get_coverage(int p_start_pos, int p_end_pos, uint32_t* coverage)
{
  //                    block1                block2
  //   |=====|======|============|===========|======|====|
  //         ^      ^                                    ^
  //  p_start_pos   |                            p_end_pos
  //            start_pos                               
  //         |0000001111111111111000000000000111111100000|
  //                           *coverage
  int len = p_end_pos-p_start_pos+1;
  for (uint32_t n = 0; n < block_starts.size(); n++) {
    int32_t from, to;
    from = block_starts[n]+start_pos-p_start_pos;
    to = block_starts[n]+start_pos-p_start_pos+block_lengths[n];
    if (from < 0)
      from = 0;
    if (to < 0)
		continue;
    else if (to > len)
		to = len;
    for (int bp=from; bp<to; bp++) {
      coverage[bp]++;
    }
  }
}
int CRead::get_last_position()
{
	if (block_starts.size()>0) // this is for some reason zero in case of softclips
		return start_pos+block_starts.back()+block_lengths.back();	
	return -1;
}

/*
 * Adds the column indices (= positions) covered by the read to 'reads' array in current row (= read).
 * These indices can be used to build up a sparse matrix of reads x positions.
 */
void CRead::get_reads_sparse(int p_start_pos, int p_end_pos, double* reads, uint32_t & reads_c, uint32_t row_idx) {
  uint32_t len = p_end_pos-p_start_pos+1;
  for (uint32_t n = 0; n < block_starts.size(); n++) {
    uint32_t from, to;
    if (block_starts[n]+start_pos-p_start_pos >= 0)
      from = block_starts[n]+start_pos-p_start_pos;
    else
      from = 0;
    if (block_starts[n]+start_pos-p_start_pos+block_lengths[n] >= 0)
      to = block_starts[n]+start_pos-p_start_pos+block_lengths[n];
    else
      to = 0;
    for (uint32_t bp=from; (bp<to)&(bp<len); bp++) {
      reads[reads_c] = row_idx+1; // row indices for sparse matrix
      reads[reads_c+1] = bp+1; // column indices for sparse matrix
      reads_c += 2;
    }
  }
}

void CRead::get_acc_splice_sites(vector<int>* acc_pos)
{
	if (strand[0]=='+')
	{
		for (uint32_t k=1;k<block_starts.size(); k++)
			acc_pos->push_back(start_pos+block_starts[k]-1);
	}
	else if (strand[0]=='-')
	{
		for (uint32_t k=1;k<block_starts.size(); k++)
			acc_pos->push_back(start_pos+block_starts[k-1]+block_lengths[k-1]-2);
	}
}

void CRead::get_don_splice_sites(vector<int>* don_pos)
{
	
	if (strand[0]=='+')
	{
		for (uint32_t k=1;k<block_starts.size(); k++)
			don_pos->push_back(start_pos+block_starts[k-1]+block_lengths[k-1]-2);
	}
	else if (strand[0]=='-')
	{
		for (uint32_t k=1;k<block_starts.size(); k++)
			don_pos->push_back(start_pos+block_starts[k]-1);
	}
}

int CRead::min_exon_len()
{
	int min = 1e8;
	for (uint32_t k=0;k<block_starts.size(); k++)
		if (block_lengths[k]<min)
			min = block_lengths[k];
	return min;
}

int CRead::max_intron_len()
{
	int max = 0;
	for (uint32_t k=1;k<block_starts.size(); k++)
		if (block_starts[k]-(block_starts[k-1]+block_lengths[k-1])>max)
			max = block_starts[k]-(block_starts[k-1]+block_lengths[k-1]);
	return max;
}

/*
 * Adds start and end of introns in the read consecutively to the 'introns' vector.
 */
void CRead::get_introns(vector<int>* introns) 
{
  for (uint32_t i=1; i<block_starts.size(); i++) 
  {
    int istart = block_starts[i-1]+block_lengths[i-1]+start_pos;
    int iend = block_starts[i]+start_pos-1;
    introns->push_back(istart);
    introns->push_back(iend);
    //fprintf(stdout, "%i intron: %d->%d\n", i, istart, iend);
  }
}
void CRead::get_introns(vector<uint32_t>* intron_starts, vector<uint32_t>* intron_ends, vector<uint32_t>* block_len1, vector<uint32_t>* block_len2) 
{
  for (uint32_t i=1; i<block_starts.size(); i++) 
  {
    uint32_t istart = block_starts[i-1]+block_lengths[i-1]+start_pos;
    uint32_t iend = block_starts[i]+start_pos-1;
    intron_starts->push_back(istart);
    intron_ends->push_back(iend);
    block_len1->push_back(block_lengths[i-1]) ;
    block_len2->push_back(block_lengths[i]) ;
  }
}

bool CRead::operator==(const CRead& read) const
{
	if (block_starts.size()!=read.block_starts.size())
		return false;
	if (block_lengths.size()!=read.block_lengths.size())
		return false;
	if (start_pos!=read.start_pos)
		return false;
	if (strand[0] != read.strand[0])
		return false;
	for (uint32_t i=0; i<block_starts.size(); i++)
		if (block_starts[i]!=read.block_starts[i])
			return false;
	for (uint32_t i=0; i<block_lengths.size(); i++)
		if (block_lengths[i]!=read.block_lengths[i])
			return false;
	return true;
}

void CRead::print()
{
	fprintf(stdout, "start_pos: %d\n", start_pos);
	fprintf(stdout, "starts:");
	for (uint32_t i=0; i<block_starts.size(); i++)
	{
		fprintf(stdout, " %d", block_starts[i]);
	}
	fprintf(stdout, "\n");

	fprintf(stdout, "lengths:");
	for (uint32_t i=0; i<block_starts.size(); i++)
	{
		fprintf(stdout, " %d", block_lengths[i]);
	}
	fprintf(stdout, "\n");
}

void CRead::set_strand(char s)
{
	delete[] strand;
	strand = new char [2];
	strand[0] = s;
	strand[1] = '0';
}

int CRead::get_mismatches()
{
	return mismatches ;
}
