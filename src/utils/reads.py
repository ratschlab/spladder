import pysam
import re
import scipy as sp
import scipy.sparse

def get_reads(fname, chr_name, start, stop, strand = None, filter = None, mapped = True, spliced = True, var_aware = None, collapse = False):
    
    infile = pysam.Samfile(fname, 'rb')

    ### vectors to build sparse matrix
    i = []
    j = []

    read_cnt = 0
    introns = []

    ### pysam query is zero based in position (results are as well), all intervals are pythonic half open
    for read in infile.fetch(chr_name, start, stop):
        if read.is_unmapped:
            continue
        if 3 in [x[0] for x in read.cigar]:
            if not spliced:
                continue
        else:
            if not mapped:
                continue

        if filter is not None:
            ### handle mismatches
            if var_aware:
                mm = sum([x[1] for x in read.tags if x[0] in ['XM', 'XG']])
            else:
                mm = sum([x[1] for x in read.tags if x[0] == 'NM'])
            if filter['mismatch'] < mm:
                continue
            ### handle min segment length
            ### remove all elements from CIGAR sting that do not 
            ### contribute to the segments (hard- and softclips and insertions)
            cig = re.sub(r'[0-9]*[HSI]', '', read.cigarstring)
            ### split the string at the introns and sum the remaining segment elements, compare to filter
            if min([sum([int(y) for y in re.split('[MD]', x)[:-1]]) for x in re.split('[0-9]*N', cig)]) < filter['min_exon_len']:
                continue
        
        ### get introns
        if 3 in [x[0] for x in read.cigar]:
            p = read.pos 
            for o in read.cigar:
                if o[0] == 3:
                    introns.append([p, p + o[1]])
                if o[0] in [0, 2, 3]:
                    p += o[1]

        ### get coverage
        for p in read.positions:
            if p - start >= 0:
                if p < start or p >= stop:
                    break
                else:
                    i.append(read_cnt)
                    j.append(p - start)

        read_cnt += 1

    i = sp.array(i)
    j = sp.array(j)

    ### construct sparse matrix
    read_matrix = scipy.sparse.coo_matrix((sp.ones(i.shape), (i, j)), shape = (read_cnt, stop - start), dtype = 'int')

    ### construct introns
    introns = sp.array(introns, dtype = 'int')

    if collapse:
        return (read_matrix.sum(axis = 0), introns)
    else:
        return (read_matrix, introns)


def add_reads_from_bam(blocks, filenames, types, filter = None, var_aware = False)
    # blocks coordinates are assumed to be in closed intervals

    if filter is None:
        filter = dict()
        filter['intron'] = 20000
        filter['exon_len'] = 8
        filter['mismatch']= 1

    if not types: 
        print 'add_reads_from_bam: nothing to do'
        return

    verbose = False
    pair = False

    pair = ('pair_coverage' in types)
    clipped = False

    for b in range(blocks.shape[0]):

        introns = None

        if verbose and  b % 10 == 0:
            print '\radd_exon_track_from_bam: %i(%i)' % (b, blocks.shape[0])
        block_len = blocks[b].stop - blocks[b].start + 1

        ## get data from bam
        if 'exon_track' in types:
            (introns, coverage) = get_all_data(blocks[b], filenames, filter=filter, var_aware=var_aware) 
        if 'mapped_exon_track' in types:
            (introns, mapped_coverage) = get_all_data(blocks[b], filenames, spliced=False, filter=filter, var_aware=var_aware) 
        if 'spliced_exon_track' in types:
            (introns, spliced_coverage) = get_all_data(blocks[b], filenames, mapped=False, filter=filter, var_aware=var_aware) 
        if 'polya_signal_track' in types:
            (introns, polya_signals) = get_all_data_uncollapsed(blocks[b], filenames, filter=filter, clipped=True, var_aware=var_aware)
        if 'end_signal_track' in types:
            (introns, read_end_signals) = get_all_data_uncollapsed(blocks[b], filenames, filter=filter, var_aware=var_aware)

        if introns is None:
            # no exon coverage needed at all
            (introns, spliced_coverage) = get_all_data(blocks[b], mapped=False, filenames, filter=filter, var_aware=var_aware)

        # add requested data to block
        for j = 1:length(types),
        tracks = sp.zeros((0, block_len))
        for type in types:
            ## add exon track to block
            ##############################################################################
            if type == 'exon_track':
                tracks = tracks.r_[tracks, coverage] 
            ## add mapped exon track to block
            ##############################################################################
            elif type == 'mapped_exon_track':
                tracks = tracks.r_[tracks, mapped_coverage] 
            ## add spliced exon track to block
            ##############################################################################
            elif type == 'spliced_exon_track':
                tracks = tracks.r_[tracks, spliced_coverage] 
            ## add intron coverage track to block
            ##############################################################################
            elif type == 'intron_track':
                intron_coverage = sp.zeros((1, block_len))
                if introns.shape[0] > 0:
                    for k in range(introns.shape[0]):
                        from_pos = max(0, introns[k, 0])
                        to_pos = min(block_len, intron_list[k, 1])
                        intron_coverage(from_pos:to_pos) += 1
                tracks = tracks.r_[tracks, intron_coverage] 
            ## compute intron list
            ##############################################################################
            elif type == 'intron_list':
                ### compute number of occurences
                introns = sort_rows(introns)
                num_introns = introns.shape[0]
                (introns, fidx) = unique_rows(introns, index=True)
                lidx = sp.r_[fidx, num_introns - 1]
                introns = sp.c_[introns, lidx - fidx + 1]
            
                if introns.shape[0] > 0:
                    s_idx = sp.argsort(introns[:, 0])
                    introns = introns[s_idx, :]
                
                if 'mincount' in filter:
                    take_idx = sp.where(introns[:, 2] >= filter['mincount']
                    introns = introns[take_idx, :]
                intron_list.append(introns)
            ## add polya signal track
            ##############################################################################
            elif type == 'polya_signal_track':
                ### get only end positions of reads
                shp = polya_signals
                end_idx = shp[0] - 1 - polya_signals[:, ::-1].argmax(axis = 1)
                polya_signals = scipy.sparse.coo_matrix((sp.ones((shp[1],)), (sp.array(range(shp[1])), end_idx)), shape = shp)
                tracks = sp.r_[tracks, polya_signals.sum(axis = 0)]
            ## add end signal track
            ##############################################################################
            elif type == 'end_signal_track':
                ### get only end positions of reads
                shp = end_signals
                end_idx = shp[0] - 1 - end_signals[:, ::-1].argmax(axis = 1)
                end_signals = scipy.sparse.coo_matrix((sp.ones((shp[1],)), (sp.array(range(shp[1])), end_idx)), shape = shp)
                tracks = sp.r_[tracks, end_signals.sum(axis = 0)]
            else: 
                print >> sys.error, 'ERROR: unknown type of data requested: %s' % type
    
    if len(types) == 1 and types[0] == 'intron_list':
        return intron_list
    elif 'intron_list' in types:
        return (tracks, intron_list)
    else:
        return tracks

#function [introns, coverage, pair_cov] = get_all_data(block, mapped, spliced, filenames, filter, clipped, var_aware) 
def get_all_data(block, filenames, mapped=True, spliced=True, filter=None, clipped=False, var_aware=False):

	block_len = block.stop - block.start
	# get all data from bam file
	coverage = sp.zeros((1, block_len))
	introns = None
	for j in range(len(filenames)):
		fname = filenames[j]
        if not os.path.exists(fname):
			print >> sys.stderr, 'add_reads_from_bam: did not find file %s' % fname
			continue

        ### TODO: implement subsampling, if needed
		contig_name = block.chr
		strand = block.strand
        ### check for filter maps -> requires uncollapsed reads
        if 'maps'in filter:
            collapse = False
        else:
            collapse = True
        ### get reads from bam file
        (coverage_tmp, introns_tmp) = get_reads(fname, contig_name, block.start, block.stop, strand, filter, mapped, spliced, var_aware, collapse)

        ### compute total coverages
        if 'maps' in filter:
            ### TODO re-implement these filters !!!
            ### apply filters
            if 'repeat_map' in filter['maps']:
                curr_idx = filter['maps']['repeat_map'][block.chr_num][block.start:block.stop]
                keep_idx = sp.where(sp.sum(coverage_tmp[:, ~curr_idx], axis = 1) > 0)[0]
                coverage_tmp = coverage_tmp[keep_idx, :]
            if 'indel_map' in filter['maps']:
                curr_idx = filter['maps']['indel_map'][block.chr_num][block.start:block.stop]
                #rm_idx = curr_idx(max(coverage_tmp, [], 2)' == 0)' == 1;
                keep_idx = ~curr_idx[coverage_tmp.max(axis = 1) == 0]
                coverage_tmp = coverage_tmp[keep_idx, :]
            if 'gene_overlap_map' in filter['maps']:
                curr_idx = filter['maps']['gene_overlap_map'][block.chr_num][block.start:block.stop]
                #rm_idx = sum(coverage_tmp(:, curr_idx)) > 0;
                keep_idx = coverage_tmp[:, curr_idx].sum(axis = 0) == 0
                coverage_tmp = coverage_tmp[keep_idx, :]
            coverage += coverage_tmp.sum(axis = 0)
        else:
            coverage += coverage_tmp

        if introns is None:
            introns = introns_tmp
        else:
            introns = sp.c_[introns, introns_tmp]

    return (introns, coverage)

#function [introns, coverage, pair_cov] = get_all_data_uncollapsed(block, mapped, spliced, filenames, filter, var_aware) 
def get_all_data_uncollapsed(block,filenames, mapped=True, spliced=True, filter=None, var_aware=False) 

	block_len = block.stop - block.start
	# get all data from bam file
	coverage = sp.zeros((1, block_len))
	introns = None

	for j in range(lenfilenames)):
		fname = filenames[j]
        if not ot.path.exists(fname):
			print >> sys.stderr, 'add_reads_from_bam: did not find file %s' % fname
			continue

        ### TODO: implement subsampling, if needed
		contig_name = block.chr
		strand = block.strand
        collapse = False
        (coverage_tmp, introns_tmp) = get_reads(fname, contig_name, block.start, block.stop, strand, filter, mapped, spliced, var_aware, collapse)
		coverage = sp.r_[coverage, coverage_tmp]

    return (introns, coverage)

