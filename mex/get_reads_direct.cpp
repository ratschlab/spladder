/* written by Jonas Behr, Regina Bohnert and Gunnar Raetsch, FML Tuebingen, Germany, 2010 */

#include <stdio.h>
#include <assert.h>
#include "sam.h"
#include "get_reads_direct.h"

#include <vector>
  using std::vector;
#include <string>
  using std::string;

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	int beg, end;
	samfile_t *in;
} tmpstruct_t;

typedef struct {
    uint64_t u, v;
} pair64_t;

static inline int is_overlap(uint32_t beg, uint32_t end, const bam1_t *b)
{
    uint32_t rbeg = b->core.pos;
    uint32_t rend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + 1;
    return (rend > beg && rbeg < end);
}

pair64_t * get_chunk_coordinates(const bam_index_t *idx, int tid, int beg, int end, int* cnt_off);

  int bam_fetch_reads(bamFile fp, const bam_index_t *idx, int tid, int beg, int end, void *data, bam_header_t* header, vector<CRead*>* reads, char strand, bool var_aware = false);

// callback for bam_plbuf_init()
static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
	//tmpstruct_t *tmp = (tmpstruct_t*)data;
	//if ((int)pos >= tmp->beg && (int)pos < tmp->end)
	//	printf("%s\t%d\t%d\n", tmp->in->header->target_name[tid], pos + 1, n);
	return 0;
}
#ifdef __cplusplus
}
#endif
int parse_sam_line(char* line, CRead* read);
//int set_strand(char c);
//void parse_cigar(bam1_t* b, CRead* read);


int get_reads_from_bam(char* filename, char* region, vector<CRead*>* reads, char strand, int lsubsample, bool var_aware)
{
	subsample = lsubsample;
	//set_strand(strand);

	srand (time(NULL));
	//srand (1234);
	tmpstruct_t tmp;
	tmp.in = samopen(filename, "rb", 0);
	if (tmp.in == 0) {
		fprintf(stderr, "Fail to open BAM file %s\n", filename);
		return 1;
	}
	int ref;
	bam_index_t *idx;
	bam_plbuf_t *buf;
	idx = bam_index_load(filename); // load BAM index
	if (idx == 0) {
		fprintf(stderr, "BAM indexing file is not available.\n");
		samclose(tmp.in);
		return 1;
	}
	bam_parse_region(tmp.in->header, region, &ref,
	                 &tmp.beg, &tmp.end); // parse the region
	if (ref < 0) {
		fprintf(stderr, "Invalid region %s\n", region);
		bam_index_destroy(idx);
		samclose(tmp.in);
		return 1;
	}

	buf = bam_plbuf_init(pileup_func, &tmp); // initialize pileup

	bam_fetch_reads(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, buf, tmp.in->header, reads, strand, var_aware);
	//fprintf(stdout, "intron_list: %d \n", intron_list->size());

	bam_plbuf_push(0, buf); // finalize pileup
	bam_index_destroy(idx);
	bam_plbuf_destroy(buf);
	samclose(tmp.in);
	return 0;
}


int bam_fetch_reads(bamFile fp, const bam_index_t *idx, int tid, int beg, int end, void *data, bam_header_t* header, vector<CRead*>* reads, char strand, bool var_aware)
{
	int n_off;
	pair64_t *off = get_chunk_coordinates(idx, tid, beg, end, &n_off);
	if (off == 0) return 0;
	{
		// retrive alignments
		uint64_t curr_off;
		int i, ret, n_seeks;
		n_seeks = 0; i = -1; curr_off = 0;
		bam1_t *b = (bam1_t*)calloc(1, sizeof(bam1_t));
		for (;;) {
			if (curr_off == 0 || curr_off >= off[i].v) { // then jump to the next chunk
				if (i == n_off - 1) break; // no more chunks
				if (i >= 0) assert(curr_off == off[i].v); // otherwise bug
				if (i < 0 || off[i].v != off[i+1].u) { // not adjacent chunks; then seek
					bam_seek(fp, off[i+1].u, SEEK_SET);
					curr_off = bam_tell(fp);
					++n_seeks;
				}
				++i;
			}
			if ((ret = bam_read1(fp, b)) > 0) {
				curr_off = bam_tell(fp);
				if (b->core.tid != tid || b->core.pos >= end) break; // no need to proceed
				else if (is_overlap(beg, end, b)) 
				{
					int rr = rand();
					if ((rr%1000 < subsample))
					{
						CRead* read = new CRead();
						parse_cigar(b, read, header, var_aware);

						if (strand == '0' || strand==read->strand[0] || read->strand[0]=='0')
						{
								reads->push_back(read);
						}
						else 
						{
							delete read;
						}
						//else if (read->strand[0]=='0'&&((b->core.flag & g_flag_off) >0))
						//{
						//	//fprintf(stdout, "(-)-strand; read->strand[0]==0, num_exons: %i \n", read->block_starts.size());
						//	// this flag means that the read has been reversed for alignment
						//	// flag bit set and (-)-strand requested
						//	reads->push_back(read);
						//}  
						//else if (read->strand[0]=='0'&&(g_flag_on>0&&(b->core.flag & g_flag_on)==0))
						//{
						//	//fprintf(stdout, "(+)-strand; read->strand[0]==0, num_exons: %i \n", read->block_starts.size());
						//	// (+)-strand requested and flag bit not set
						//	reads->push_back(read);
						//}  
					}
				}
			} else break; // end of file
		}
//		fprintf(stderr, "[bam_fetch] # seek calls: %d\n", n_seeks);
		bam_destroy1(b);
	}
	free(off);
	return 0;
}

void parse_cigar(bam1_t* b, CRead* read, bam_header_t* header, bool var_aware)
{
	read->left = (b->core.flag & left_flag_mask) >0;
	read->right = (b->core.flag & right_flag_mask) >0;
	read->reverse = (b->core.flag & reverse_flag_mask) >0;
	read->secondary = (b->core.flag & secondary_flag_mask) >0;

	read->start_pos = b->core.pos+1;
	read->set_strand('0');
	read->read_id = new char[100];
	//sprintf(read->read_id, "%s\0", bam1_qname(b));
	sprintf(read->read_id, "%s", bam1_qname(b));

	for (int k = 0; k < b->core.n_cigar; ++k) 
	{
		int op = bam1_cigar(b)[k] & BAM_CIGAR_MASK; // operation
		int l = bam1_cigar(b)[k] >> BAM_CIGAR_SHIFT; // length
		//int op = bam_cigar_op(bam1_cigar(b)[k]); // operation
		//int l = bam_cigar_oplen(bam1_cigar(b)[k]); // length
		//fprintf(stdout, "op:%d l:%d\n", op, l);
		if (op == BAM_CMATCH) 
		{
			if (k==0)
			{
				read->block_lengths.push_back(l);
				read->block_starts.push_back(0);
			}
			else
			{
				int op_prev = bam1_cigar(b)[k-1] & BAM_CIGAR_MASK; 
				int l_prev = bam1_cigar(b)[k-1] >> BAM_CIGAR_SHIFT;
				if (op_prev==BAM_CREF_SKIP)// intron before
				{
					if (read->block_lengths.size()>=1)
					{
						int last_block_start = (*(read->block_starts.end()-1));
						int intron_start = last_block_start+(*(read->block_lengths.end()-1));
						read->block_lengths.push_back(l);
						read->block_starts.push_back(intron_start+l_prev);
					}
					else
					{
						// start of first block was not a match
						read->block_lengths.push_back(l);
						read->block_starts.push_back(0);
					}
				}
				else
				{
					if (read->block_lengths.size()>=1)
						(*(read->block_lengths.end()-1))+=l;
					else
					{
						read->block_lengths.push_back(l);
						read->block_starts.push_back(0);
					}
				}
			}
		}
		else if (op == BAM_CDEL) 
		{
			if (k>0 && read->block_lengths.size()>=1)
				(*(read->block_lengths.end()-1))+=l;
		} 
		else if (op == BAM_CREF_SKIP)//intron
		{}
		else if (op == BAM_CINS)
		{}
		else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)
		{
			read->is_clipped = true;
		}
	}
	// parse auxiliary data
    uint8_t* s = bam1_aux(b);
	uint8_t* end = b->data + b->data_len; 
	while (s < end) 
	{
		 uint8_t type, key[2];
		 key[0] = s[0]; key[1] = s[1];
		 s += 2; type = *s; ++s;
		 //fprintf(stdout, "\n%c%c:%c\n", key[0], key[1], type);
		 if (type == 'A')
		 {
			if ( key[0] =='X' && key[1] == 'S')
			{
				read->set_strand((char) *s);
			}
		 	++s;
		 }
		else if (type == 'C')
		{ 
			if ( key[0] =='H' && key[1] == '0')
			{
				uint8_t matches = *s;
				read->matches = (int) matches;
			}
            if ( var_aware ) {
                if ( key[0] =='X' && key[1] == 'M')
                {
                    uint8_t mismatches = *s;
                    read->mismatches += (int) mismatches;
                }
                if ( key[0] =='X' && key[1] == 'G')
                {
                    uint8_t mismatches = *s;
                    read->mismatches += (int) mismatches;
                }
            } else {
                if ( key[0] =='N' && key[1] == 'M')
                {
                    uint8_t mismatches = *s;
                    read->mismatches = (int) mismatches;
                }
            }
			if ( key[0] =='H' && key[1] == 'I')
			{
				uint8_t mai = *s;
				read->multiple_alignment_index = (int) mai;
			}

			++s;
		}
		else if (type == 'c') { ++s; }
		else if (type == 'S') { s += 2;	}
		else if (type == 's') { s += 2;	}
		else if (type == 'I') { s += 4; }
		else if (type == 'i') { s += 4; }
		else if (type == 'f') { s += 4;	}
		else if (type == 'd') { s += 8;	}
		else if (type == 'Z') { ++s; }
		else if (type == 'H') { ++s; }
	}
}

//int set_strand(char c)
//{
//	if (c=='+')
//	{
//		char* fl = (char*) "0x0010";
//		g_flag_on = strtol(fl, 0, 0);
//		g_flag_off = 0;
//	}
//	else if (c=='-')
//	{
//		char* fl = (char*) "0x0010";
//		g_flag_off = strtol(fl, 0, 0);
//		g_flag_on = 0;
//	}
//	return 0;
//}

