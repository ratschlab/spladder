import pysam
import re
import scipy as sp
import scipy.sparse
import copy
import pdb
import time

from utils import *
from init import *

def get_reads(fname, chr_name, start, stop, strand = None, filter = None, mapped = True, spliced = True, var_aware = None, collapse = False):
    
    infile = pysam.Samfile(fname, 'rb')

    ### vectors to build sparse matrix
    i = []
    j = []

    read_cnt = 0
    introns = []

    ### pysam query is zero based in position (results are as well), all intervals are pythonic half open
    for read in infile.fetch(chr_name, start, stop, until_eof=True):

        if read.is_unmapped:
            continue
        is_spliced = ('N' in read.cigarstring)
        if is_spliced:
            if not spliced:
                continue
        elif not mapped:
            continue

        tags = dict(read.tags)

        if filter is not None:
            ### handle mismatches
            if var_aware:
                if filter['mismatch'] < (tags['XM'] + tags['XG']):
                    continue
            else:
                if filter['mismatch'] < tags['NM']:
                    continue

            if is_spliced:
                ### handle min segment length
                ### remove all elements from CIGAR sting that do not 
                ### contribute to the segments (hard- and softclips and insertions)
                cig = re.sub(r'[0-9]*[HSI]', '', read.cigarstring)
                ### split the string at the introns and sum the remaining segment elements, compare to filter
                if min([sum([int(y) for y in re.split('[MD]', x)[:-1]]) for x in re.split('[0-9]*N', cig)]) <= filter['exon_len']:
                    continue
        
        ### check strand information
        if strand is not None:
            try:
                if tags['XS'] != strand:
                    continue
            except KeyError:
                pass

        ### get introns and covergae
        p = read.pos 
        for o in read.cigar:
            if o[0] == 3:
                introns.append([p, p + o[1]])
            if o[0] in [0, 2]:
                r = range(max(p-start, 0), min(p + o[1] - start, stop - start))
                i.extend([read_cnt] * len(r))
                j.extend(r)
                #for pp in range(p, p + o[1]):
                #    if pp - start >= 0 and pp < stop:
                #        i.append(read_cnt)
                #        j.append(pp - start)
            if o[0] in [0, 2, 3]:
                p += o[1]

        ### the follwoing is new behavior and gonne come in the next version --> deletions are not counted towards coverage
        #### get coverage
        #for p in read.positions:
        #    if p - start >= 0:
        #        if p >= stop:
        #            break
        #        else:
        #            i.append(read_cnt)
        #            j.append(p - start)

        read_cnt += 1

    i = sp.array(i)
    j = sp.array(j)

    ### construct sparse matrix
    read_matrix = scipy.sparse.coo_matrix((sp.ones(i.shape), (i, j)), shape = (read_cnt, stop - start), dtype = 'int')

    ### construct introns
    if len(introns) > 0:
        introns = sp.array(introns, dtype = 'int')
    else:
        introns = sp.zeros((0, 2), dtype='int')

    if collapse:
        return (read_matrix.sum(axis = 0), introns)
    else:
        return (read_matrix, introns)


def add_reads_from_bam(blocks, filenames, types, filter = None, var_aware = False):
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

    if type(blocks).__module__ != 'numpy':
        blocks = sp.array([blocks])

    for b in range(blocks.shape[0]):

        introns = None

        if verbose and  b % 10 == 0:
            print '\radd_exon_track_from_bam: %i(%i)' % (b, blocks.shape[0])
        block_len = blocks[b].stop - blocks[b].start

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
            (introns, spliced_coverage) = get_all_data(blocks[b], filenames, mapped=False, filter=filter, var_aware=var_aware)

        introns = sort_rows(introns)

        # add requested data to block
        tracks = sp.zeros((0, block_len))
        intron_list = []
        for ttype in types:
            ## add exon track to block
            ##############################################################################
            if ttype == 'exon_track':
                tracks = sp.r_[tracks, coverage] 
            ## add mapped exon track to block
            ##############################################################################
            elif ttype == 'mapped_exon_track':
                tracks = sp.r_[tracks, mapped_coverage] 
            ## add spliced exon track to block
            ##############################################################################
            elif ttype == 'spliced_exon_track':
                tracks = sp.r_[tracks, spliced_coverage] 
            ## add intron coverage track to block
            ##############################################################################
            elif ttype == 'intron_track':
                intron_coverage = sp.zeros((1, block_len))
                if introns.shape[0] > 0:
                    for k in range(introns.shape[0]):
                        from_pos = max(0, introns[k, 0])
                        to_pos = min(block_len, intron_list[k, 1])
                        intron_coverage[from_pos:to_pos] += 1
                tracks = sp.r_[tracks, intron_coverage] 
            ## compute intron list
            ##############################################################################
            elif ttype == 'intron_list':
                if introns.shape[0] > 0:
                    ### compute number of occurences
                    introns = sort_rows(introns)
                    num_introns = introns.shape[0]
                    (introns, fidx) = unique_rows(introns, index=True)
                    lidx = sp.r_[fidx[1:] - 1, num_introns - 1]
                    introns = sp.c_[introns, lidx - fidx + 1]

                    ### filter introns for location relative to block
                    ### this is legacy behavior for matlab versions!
                    k_idx = sp.where((introns[:, 0] > blocks[0].start) & (introns[:, 1] < blocks[0].stop))[0]
                    introns = introns[k_idx, :]

                    if introns.shape[0] > 0:
                        s_idx = sp.argsort(introns[:, 0])
                        introns = introns[s_idx, :]
                    
                    if 'mincount' in filter:
                        take_idx = sp.where(introns[:, 2] >= filter['mincount'])[0]
                        if take_idx.shape[0] > 0:
                            intron_list.append(introns[take_idx, :])
                        else:
                            intron_list.append(sp.zeros((0, 3), dtype='int'))
                else:
                    intron_list.append(sp.zeros((0, 3), dtype='int'))
            ## add polya signal track
            ##############################################################################
            elif ttype == 'polya_signal_track':
                ### get only end positions of reads
                shp = polya_signals
                end_idx = shp[0] - 1 - polya_signals[:, ::-1].argmax(axis = 1)
                polya_signals = scipy.sparse.coo_matrix((sp.ones((shp[1],)), (sp.array(range(shp[1])), end_idx)), shape = shp)
                tracks = sp.r_[tracks, polya_signals.sum(axis = 0)]
            ## add end signal track
            ##############################################################################
            elif ttype == 'end_signal_track':
                ### get only end positions of reads
                shp = end_signals
                end_idx = shp[0] - 1 - end_signals[:, ::-1].argmax(axis = 1)
                end_signals = scipy.sparse.coo_matrix((sp.ones((shp[1],)), (sp.array(range(shp[1])), end_idx)), shape = shp)
                tracks = sp.r_[tracks, end_signals.sum(axis = 0)]
            else: 
                print >> sys.stderr, 'ERROR: unknown type of data requested: %s' % ttype
    
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

    if isinstance(filenames, str):
        filenames = [filenames]

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
            introns = sp.r_[introns, introns_tmp]

    return (introns, coverage)

#function [introns, coverage, pair_cov] = get_all_data_uncollapsed(block, mapped, spliced, filenames, filter, var_aware) 
def get_all_data_uncollapsed(block,filenames, mapped=True, spliced=True, filter=None, var_aware=False):

    block_len = block.stop - block.start
    # get all data from bam file
    coverage = sp.zeros((1, block_len))
    introns = None

    for j in range(lenfilenames):
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

def get_intron_list(genes, CFG):

    # (author) Georg Zeller & Gunnar Raetsch, Mpi Tuebingen, Germany, 2009
    # (author) Andre Kahles, MSKCC NYC, USA, 2013 

    #function introns = get_intron_list(genes, CFG)

    ### form chunks for quick sorting
    strands = ['+', '-']
    chunks = sp.array([[x.chr_num, strands.index(x.strand), x.start, x.stop] for x in genes], dtype = 'int')
    #chunks = sp.c_[sp.array([x.chr_num for x in genes], dtype = 'int'), sp.array([x.strand for x in genes], dtype = 'int'), sp.array([x.start for x in genes], dtype = 'int'), sp.array([x.stop for x in genes], dtype = 'int')]
    (chunks, chunk_idx) = sort_rows(chunks, index=True)

    introns = sp.zeros((chunks.shape[0], 2), dtype = 'object')
    introns[:] = None

    ### collect all possible combinations of contigs and strands
    regions = init_regions(CFG['bam_fnames'])
    keepidx = sp.where(sp.in1d(sp.array([x.chr_num for x in regions]), sp.unique(sp.array([x.chr_num for x in genes]))))[0]
    regions = regions[keepidx]

    c = 0
    num_introns_filtered = 0
    t0 = time.time()

    for j in range(regions.shape[0]):
        chr = regions[j].chr
        chr_num = regions[j].chr_num
        s = strands.index(regions[j].strand)
        
        # fill the chunks on the corresponding chromosome
        while c < chunks.shape[0]:
            if chunks[c, 0] > chr_num or chunks[c, 1] > s:
                break
            if chunks[c, 0] != chr_num:
                print >> sys.stderr, 'ERROR: c logic seems wrong' 
                sys.exit(1)

            if CFG['verbose'] and (c+1) % 100 == 0:
                t1 = time.time()
                print >> sys.stdout, '%i (%i) genes done (%i introns taken) ... took %i secs' % (c+1, chunks.shape[0], num_introns_filtered, t1 - t0)
                t0 = t1

            gg = sp.array([copy.copy(genes[chunk_idx[c]])], dtype='object')
            gg[0].strand = strands[s]
            gg[0].start = max(gg[0].start - 5000, 1)
            gg[0].stop = gg[0].stop + 5000
            assert(gg[0].chr == chr)

            intron_list_tmp = add_reads_from_bam(gg, CFG['bam_fnames'], ['intron_list'], CFG['read_filter'], CFG['var_aware'])
            num_introns_filtered += intron_list_tmp[0].shape[0]
            introns[chunk_idx[c], s] = intron_list_tmp[0]

            c += 1

    for j in range(introns.shape[0]):
        if introns[j, 0] is None:
            introns[j, 0] = sp.zeros((0, 3), dtype='int')
        if introns[j, 1] is None:
            introns[j, 1] = sp.zeros((0, 3), dtype='int')

    return introns
