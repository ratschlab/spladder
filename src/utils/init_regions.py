import os
import sys
import pysam
import scipy as sp

def init_regions(fn_bams):
    # regions=init_regions(fn_bams)

    regions = []
    processed = []
    cnt = 0

    for i in range(len(fn_bams)):
        if not os.path.exists(fn_bams[i]):
            continue
        else:
            ### load bamfile
            IN = pysam.Samfile(fn_bams[i], 'rb')
            header_info = IN.header['SQ']
            strands = ['+', '-']
            for c in range(len(header_info)):
                if not (header_info[c]['SN'], header_info[c]['LN']) in processed:
                    for s in strands:
                        region = Region()
                        region.chr = header_info[c]['SN'])
                        region.chr_num = c
                        region.strand = s
                        region.start = 1
                        region.stop = header_info[c]['LN']
                        region.id = cnt
                        regions.append(region)
                        cnt += 1
                    processed.append((header_info[c]['SN'], header_info[c]['LN']))
            IN.close()
        break

    return sp.array(regions, dtype = 'object')
