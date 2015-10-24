if __package__ is None:
    __package__ ='modules'

import cPickle
import math
import h5py
import scipy as sp

from .classes.segmentgraph import Segmentgraph
from .classes.counts import Counts
from reads import *
import rproc as rp

def count_graph_coverage(genes, fn_bam=None, CFG=None, fn_out=None):
# [counts] = count_graph_coverage(genes, fn_bam, CFG, fn_out)

    if fn_bam is None and isinstance(genes, dict):
        PAR = genes
        genes = PAR['genes']
        fn_bam = PAR['fn_bam']
        if 'fn_out' in PAR:
            fn_out = PAR['fn_out'] 
        CFG = PAR['CFG']

    if not isinstance(fn_bam, list):
        fn_bam = [fn_bam]
    counts = sp.zeros((len(fn_bam), genes.shape[0]), dtype='object')

    intron_tol = 0 

    sys.stdout.write('genes: %i\n' % genes.shape[0])
    for f in range(counts.shape[0]):
        sys.stdout.write('sample %i/%i\n' % (f + 1, counts.shape[0])) 

        bam_cache = None

        ### iterate over all genes and generate counts for
        ### the segments in the segment graph
        ### and the splice junctions in the splice graph
        for i in range(genes.shape[0]):
            sys.stdout.write('.')
            if i > 0 and i % 50 == 0:
                sys.stdout.write('%i\n' % i)
            gg = genes[i]
            if gg.segmentgraph.is_empty():
                gg.segmentgraph = Segmentgraph(gg)
            gg.start = gg.segmentgraph.segments.ravel().min()
            gg.stop = gg.segmentgraph.segments.ravel().max()

            counts[f, i] = Counts(gg.segmentgraph.segments.shape[1])

            if CFG['bam_to_sparse'] and (fn_bam[f].endswith('npz') or os.path.exists(re.sub(r'bam$', '', fn_bam[f]) + 'npz')):
                ### load counts from summary file
                if bam_cache is None:
                    bam_cache = dict()
                    if fn_bam[f].endswith('npz'):
                        tmp = sp.load(fn_bam[f])
                    else:
                        tmp = sp.load(re.sub(r'bam$', '', fn_bam[f]) + 'npz')
                    ### re-built sparse matrix
                    for c in sp.unique([re.sub(r'_reads_dat$', '', x) for x in tmp if x.endswith('_reads_dat')]):
                        bam_cache[c + '_reads'] = scipy.sparse.coo_matrix((tmp[c + '_reads_dat'], (tmp[c + '_reads_row'], tmp[c + '_reads_col'])), shape=tmp[c + '_reads_shp'], dtype='uint32').tocsc()
                        bam_cache[c + '_introns_m'] = tmp[c + '_introns_m']
                        bam_cache[c + '_introns_p'] = tmp[c + '_introns_p']
                    del tmp

                if bam_cache[gg.chr + '_reads'].shape[0] == 0:
                    tracks = sp.zeros((1, gg.stop - gg.start), dtype='int')
                elif bam_cache[gg.chr + '_reads'].shape[0] > 1:
                    tracks = bam_cache[gg.chr + '_reads'][[0, 1 + int(gg.strand == '-')], gg.start:gg.stop].todense() 
                else:
                    tracks = bam_cache[gg.chr + '_reads'][:, gg.start:gg.stop].todense() 

                if bam_cache[c + '_introns_m'].shape[0] > 0:
                    if gg.strand == '-':
                        intron_list = get_intron_range(bam_cache[gg.chr + '_introns_m'], gg.start, gg.stop)
                    else:
                        intron_list = get_intron_range(bam_cache[gg.chr + '_introns_p'], gg.start, gg.stop)
                else:
                    intron_list = get_intron_range(bam_cache[gg.chr + '_introns_p'], gg.start, gg.stop)
            else:
                ### add RNA-seq evidence to the gene structure
                #(tracks, intron_list) = add_reads_from_bam(gg, fn_bam[f], ['exon_track','intron_list'], CFG['read_filter'], CFG['var_aware'], CFG['primary_only']);
                (tracks, intron_list) = add_reads_from_bam(gg, fn_bam[f], ['exon_track','intron_list'], None, CFG['var_aware'], CFG['primary_only']);
                intron_list = intron_list[0] ### TODO

            ### extract mean exon coverage for all segments
            for j in range(gg.segmentgraph.segments.shape[1]):
                idx = sp.arange(gg.segmentgraph.segments[0, j], gg.segmentgraph.segments[1, j]) - gg.start
                counts[f, i].segments[j] = sp.mean(sp.sum(tracks[:, idx], axis=0))
                counts[f, i].seg_pos[j] = sp.sum(sp.sum(tracks[:, idx], axis=0) > 0)

            k, l = sp.where(gg.segmentgraph.seg_edges == 1)

            ### there are no introns to count
            if intron_list.shape[0] == 0:
                for m in range(k.shape[0]):
                    counts[f, i].edges = sp.atleast_2d(sp.array([sp.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), 0]))
                continue

            ### extract intron counts 
            for m in range(k.shape[0]):
                idx = sp.where((sp.absolute(intron_list[:, 0] - gg.segmentgraph.segments[1, k[m]]) <= intron_tol) & (sp.absolute(intron_list[:, 1] - gg.segmentgraph.segments[0, l[m]]) <= intron_tol))[0]
                if counts[f, i].edges.shape[0] == 0:
                    if idx.shape[0] > 0:
                        counts[f, i].edges = sp.atleast_2d(sp.array([sp.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), sp.sum(intron_list[idx, 2])]))
                    else:
                        counts[f, i].edges = sp.atleast_2d(sp.array([sp.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), 0]))
                else:
                    if idx.shape[0] > 0:
                        counts[f, i].edges = sp.r_[counts[f, i].edges, sp.atleast_2d(sp.array([sp.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), sp.sum(intron_list[idx, 2])]))]
                    else:
                        counts[f, i].edges = sp.r_[counts[f, i].edges, sp.atleast_2d(sp.array([sp.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), 0]))]

    if fn_out is not None:
        cPickle.dump(counts, open(fn_out, 'w'), -1)
    else:
        return counts



def count_graph_coverage_wrapper(fname_in, fname_out, CFG):

    (genes, inserted) = cPickle.load(open(fname_in, 'r'))
    
    if genes[0].segmentgraph.is_empty():
        for g in genes:
            g.segmentgraph = Segmentgraph(g)
        cPickle.dump((genes, inserted), open(fname_in, 'w'), -1)

    counts = dict()
    counts['segments'] = []
    counts['seg_pos'] = []
    counts['gene_ids_segs'] = []
    counts['edges'] = []
    counts['gene_ids_edges'] = []
    if not CFG['rproc']:
        for s_idx in range(CFG['strains'].shape[0]):
            print '\n%i/%i' % (s_idx + 1, CFG['strains'].shape[0])
            if s_idx == 0:
                counts_tmp = count_graph_coverage(genes, CFG['bam_fnames'][s_idx], CFG)
            else:
                counts_tmp = sp.r_[sp.atleast_2d(counts_tmp), count_graph_coverage(genes, CFG['bam_fnames'][s_idx], CFG)]

        for c in range(counts_tmp.shape[1]):
            counts['segments'].append(sp.hstack([sp.atleast_2d(x.segments).T for x in counts_tmp[:, c]]))
            counts['seg_pos'].append(sp.hstack([sp.atleast_2d(x.seg_pos).T for x in counts_tmp[:, c]]))
            counts['gene_ids_segs'].append(sp.ones((sp.atleast_2d(counts_tmp[0, c].seg_pos).shape[1], 1), dtype='int') * c)
            tmp = [sp.atleast_2d(x.edges) for x in counts_tmp[:, c] if x.edges.shape[0] > 0]
            if len(tmp) == 0:
                continue
            tmp = sp.hstack(tmp)
            if tmp.shape[0] > 0:
                counts['edges'].append(sp.c_[tmp[:, 0], tmp[:, range(1, tmp.shape[1], 2)]])
                counts['gene_ids_edges'].append(sp.ones((tmp.shape[0], 1), dtype='int') * c)
    else:
        ### have an adaptive chunk size, that takes into account the number of strains (take as many genes as it takes to have ~10K strains)
        chunksize = int(max(1, math.floor(10000 / len(CFG['strains']))))

        jobinfo = []

        PAR = dict()
        PAR['CFG'] = CFG
        for c_idx in range(0, genes.shape[0], chunksize):
            cc_idx = min(genes.shape[0], c_idx + chunksize)
            fn = fname_out.replace('.pickle', '.chunk_%i_%i.pickle' % (c_idx, cc_idx))
            if os.path.exists(fn):
                continue
            else:
                print 'submitting chunk %i to %i' % (c_idx, cc_idx)
                PAR['genes'] = genes[c_idx:cc_idx]
                PAR['fn_bam'] = CFG['bam_fnames']
                PAR['fn_out'] = fn
                PAR['CFG'] = CFG
                jobinfo.append(rp.rproc('count_graph_coverage', PAR, 6000, CFG['options_rproc'], 60*48))

        rp.rproc_wait(jobinfo, 30, 1.0, -1)

        ### merge results from count chunks
        if 'verbose' in CFG and CFG['verbose']:
            print '\nCollecting count data from chunks ...\n'

        for c_idx in range(0, genes.shape[0], chunksize):
            cc_idx = min(genes.shape[0], c_idx + chunksize)
            if 'verbose' in CFG and CFG['verbose']:
                print 'collecting chunk %i-%i (%i)' % (c_idx, cc_idx, genes.shape[0])
            fn = fname_out.replace('.pickle', '.chunk_%i_%i.pickle' % (c_idx, cc_idx))
            if not os.path.exists(fn):
                print >> sys.stderr, 'ERROR: Not all chunks in counting graph coverage completed!'
                sys.exit(1)
            else:
                counts_tmp = cPickle.load(open(fn, 'r'))
                for c in range(counts_tmp.shape[1]):
                    counts['segments'].append(sp.hstack([sp.atleast_2d(x.segments).T for x in counts_tmp[:, c]]))
                    counts['seg_pos'].append(sp.hstack([sp.atleast_2d(x.seg_pos).T for x in counts_tmp[:, c]]))
                    counts['gene_ids_segs'].append(sp.ones((sp.atleast_2d(counts_tmp[0, c].seg_pos).shape[1], 1), dtype='int') * (c_idx + c))
                    tmp = [sp.atleast_2d(x.edges) for x in counts_tmp[:, c] if x.edges.shape[0] > 0]
                    if len(tmp) == 0:
                        continue
                    tmp = sp.hstack(tmp)
                    if tmp.shape[0] > 0:
                        counts['edges'].append(sp.c_[tmp[:, 0], tmp[:, range(1, tmp.shape[1], 2)]])
                        counts['gene_ids_edges'].append(sp.ones((tmp.shape[0], 1), dtype='int') * (c_idx + c))

    for key in counts:
        counts[key] = sp.vstack(counts[key]) if len(counts[key]) > 0 else counts[key]
    if len(counts['edges']) > 0:
        counts['edge_idx'] = counts['edges'][:, 0]
        counts['edges'] = counts['edges'][:, 1:]
    else:
        counts['edge_idx'] = sp.array([])
        counts['edges'] = sp.array([])
    counts['seg_len'] = sp.hstack([x.segmentgraph.segments[1, :] - x.segmentgraph.segments[0, :] for x in genes]).T

    ### write result data to hdf5
    h5fid = h5py.File(fname_out, 'w')
    h5fid.create_dataset(name='gene_names', data=sp.array([x.name for x in genes], dtype='str'))
    h5fid.create_dataset(name='strains', data=CFG['strains'])
    for key in counts:
        h5fid.create_dataset(name=key, data=counts[key])
    h5fid.close()

def get_intron_range(introns, start, stop):
    """Given a sorted list of introns, return the subset of introns that
       overlaps that start stop interval"""

    if introns.shape[0] == 0:
        return introns

    idx = sp.where(((introns[:, 0] < stop) & (introns[:, 1] > start)) |
                   ((introns[:, 0] < start) & (introns[:, 1] > stop)))[0]

    return introns[idx, :]

