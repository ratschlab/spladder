/* written by Jonas Behr, Regina Bohnert and Gunnar Raetsch, FML Tuebingen, Germany, 2010 */

#ifndef __GET_READS_DIRECT_H__
#define __GET_READS_DIRECT_H__

#include <vector>
  using std::vector;
#include "read.h"
#include <stdlib.h>
#include "sam.h"

//static int g_flag_on = 0, g_flag_off = 0;
static int left_flag_mask = strtol((char*) "0x40", 0, 0);
static int right_flag_mask = strtol((char*) "0x80", 0, 0);
static int reverse_flag_mask = strtol((char*) "0x10", 0, 0);

static int subsample = 1000; 
//static int collapse = 0;

int get_reads_from_bam(char* filename, char* region, vector<CRead*>* reads, char strand, int lsubsample);
void parse_cigar(bam1_t* b, CRead* read, bam_header_t* header);

#endif
