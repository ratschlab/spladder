import cPickle
import math
import h5py
import scipy as sp

from .classes.segmentgraph import Segmentgraph
from .classes.counts import Counts
from reads import *

def count_graph_coverage(genes, fn_bam=None, CFG=None, fn_out=None):
# [counts] = count_graph_coverage(genes, fn_bam, CFG, fn_out)

    if fn_bam is None and instanceof(genes, dict):
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

    for f in range(counts.shape[0]):
        ### iterate over all genes and generate counts for
        ### the segments in the segment graph
        ### and the splice junctions in the splice graph
        for i in range(genes.shape[0]):
            print '.',
            if i > 0 and i % 50 == 0:
                print '%i' % i
            gg = genes[i]
            if gg.segmentgraph is None:
                gg.segmentgraph = Segmentgraph(gg)
            gg.start = gg.segmentgraph.segments.ravel().min()
            gg.stop = gg.segmentgraph.segments.ravel().max()

            ### add RNA-seq evidence to the gene structure
            (tracks, intron_list) = add_reads_from_bam(gg, fn_bam[f], ['exon_track','intron_list'], CFG['read_filter'], CFG['var_aware']);
            intron_list = intron_list[0] ### TODO

            ### extract mean exon coverage for all segments
            counts[f, i] = Counts(gg.segmentgraph.segments.shape[1])

            for j in range(gg.segmentgraph.segments.shape[1]):
                idx = sp.arange(gg.segmentgraph.segments[0, j], gg.segmentgraph.segments[1, j]) - gg.start
                counts[f, i].segments[j] = sp.mean(sp.sum(tracks[:, idx], axis=0))
                counts[f, i].seg_pos[j] = sp.sum(sp.sum(tracks[:, idx], axis=0) > 0)

            ### extract intron counts 
            k, l = sp.where(gg.segmentgraph.seg_edges == 1)
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
    
    if genes[0].segmentgraph is None:
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
            tmp = sp.hstack([sp.atleast_2d(x.edges) for x in counts_tmp[:, c]])
            if tmp.shape[0] > 0:
                counts['edges'].append(sp.c_[tmp[:, 0], tmp[:, range(1, tmp.shape[1], 2)]])
                counts['gene_ids_edges'].append(sp.ones((tmp.shape[0], 1), dtype='int') * c)
    else:
        ### have an adaptive chunk size, that takes into account the number of strains (take as many genes as it takes to have ~10K strains)
        chunksize = max(1, math.floor(10000 / len(CFG['strains'])));

        jobinfo = rproc_empty(0)

        PAR = dict()
        PAR['CFG'] = CFG
        for c_idx in range(0, genes.shape[0], chunksize):
            cc_idx = min(genes.shape[0], c_idx + chunksize - 1)
            fn = fname_out.replace('.mat', '.chunk_%i_%i.mat' % (c_idx, cc_idx))
            if os.path.exists(fn):
                continue
            else:
                print 'submitting chunk %i to %i' % (c_idx, cc_idx)
                PAR['genes'] = genes[c_idx:cc_idx]
                PAR['fn_bam'] = CFG['bam_fnames']
                PAR['fn_out'] = fn
                PAR['CFG'] = CFG
                jobinfo.append(rproc('count_graph_coverage', PAR, 30000, CFG['options_rproc'], 60))

        jobinfo = rproc_wait(jobinfo, 30, 1, -1)

        ### merge results
        counts = sp.zeros((CFG['strains'].shape[0], genes.shape[0]), dtype='object')
        for c_idx in range(0, genes.shape[0], chunksize):
            cc_idx = min(genes.shape[0], c_idx + chunksize - 1)
            fn = fname_out.replace('.mat', '.chunk_%i_%i.mat' % (c_idx, cc_idx))
            if not os.path.exists(fn):
                print >> sys.stderr, 'ERROR: Not all chunks in counting graph coverage completed!'
                sys.exit(1)
            else:
                counts_tmp = cPickle.load(open(fn, 'r'))
                for c in range(counts_tmp.shape[1]):
                    counts['segments'].append(sp.hstack([sp.atleast_2d(x.segments).T for x in counts_tmp[:, c]]))
                    counts['seg_pos'].append(sp.hstack([sp.atleast_2d(x.seg_pos).T for x in counts_tmp[:, c]]))
                    counts['gene_ids_segs'].append(sp.ones((sp.atleast_2d(counts_tmp[0, c].seg_pos).shape[1], 1), dtype='int') * c)
                    tmp = sp.hstack([sp.atleast_2d(x.edges) for x in counts_tmp[:, c]])
                    if tmp.shape[0] > 0:
                        counts['edges'].append(sp.c_[tmp[:, 0], tmp[:, range(1, tmp.shape[1], 2)]])
                    counts['gene_ids_edges'].append(sp.ones((tmp.shape[0], 1), dtype='int') * c)

    for key in counts:
        if len(counts[key]) > 0:
            counts[key] = sp.vstack(counts[key])
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

