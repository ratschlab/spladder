#if __package__ is None:
#    __package__ ='modules'

import pickle
import math
import h5py
import numpy as np
import re

from .classes.segmentgraph import Segmentgraph
from .classes.counts import Counts
from .helpers import log_progress, decodeUTF8, codeUTF8, get_filename
from .reads import *
from .hdf5 import appendToHDF5

### intermediate fix to load pickle files stored under previous version
from .classes import gene as cgene
from .classes import splicegraph as csplicegraph
from .classes import segmentgraph as csegmentgraph
sys.modules['modules.classes.gene'] = cgene
sys.modules['modules.classes.splicegraph'] = csplicegraph
sys.modules['modules.classes.segmentgraph'] = csegmentgraph
### end fix

def count_graph_coverage(genes, fn_bam=None, options=None, fn_out=None):

    if fn_bam is None and isinstance(genes, dict):
        PAR = genes
        genes = PAR['genes']
        fn_bam = PAR['fn_bam']
        if 'fn_out' in PAR:
            fn_out = PAR['fn_out'] 
        options = PAR['options']

    if hasattr(genes[0], 'splicegraph_edges_data'):
        for gg in genes:
            gg.from_sparse()

    if not isinstance(fn_bam, list):
        fn_bam = [fn_bam]
    counts = np.zeros((len(fn_bam), genes.shape[0]), dtype='object')

    intron_tol = 0 

    sys.stdout.write('genes: %i\n' % genes.shape[0])
    for f in range(counts.shape[0]):
        sys.stdout.write('\nsample %i/%i\n' % (f + 1, counts.shape[0])) 

        ### iterate over all genes and generate counts for
        ### the segments in the segment graph
        ### and the splice junctions in the splice graph
        ### iterate per contig, so the bam caching works better
        contigs = np.array([x.chr for x in genes])
        for contig in np.unique(contigs):
            contig_idx = np.where(contigs == contig)[0]
            bam_cache = dict()
            print('\ncounting %i genes on contig %s' % (contig_idx.shape[0], contig))
            for ii,i in enumerate(contig_idx):
                log_progress(ii, contig_idx.shape[0], 50)
                gg = genes[i]
                if gg.segmentgraph.is_empty():
                    gg.segmentgraph = Segmentgraph(gg)
                gg.start = gg.segmentgraph.segments.ravel().min()
                gg.stop = gg.segmentgraph.segments.ravel().max()

                counts[f, i] = Counts(gg.segmentgraph.segments.shape[1])

                if options.sparse_bam and \
                  (fn_bam[f].endswith('npz') or \
                   os.path.exists(re.sub(r'\.[bB][aA][mM]|\.[cC][rR][aA][mM]$', '', fn_bam[f]) + '.npz') or \
                   fn_bam[f].endswith('hdf5') or \
                   os.path.exists(re.sub(r'\.[bB][aA][mM]|\.[cC][rR][aA][mM]$', '', fn_bam[f]) + '.hdf5')):
                    ### make sure that we query the right contig from cache
                    assert(gg.chr == contig)
                    (tracks, intron_list) = add_reads_from_sparse_bam(gg, fn_bam[f], contig, options.confidence, types=['exon_track','intron_list'], filter=None, cache=bam_cache)
                else:
                    ### add RNA-seq evidence to the gene structure
                    (tracks, intron_list) = add_reads_from_bam(gg, fn_bam[f], ['exon_track','intron_list'], None, options.var_aware, options.primary_only, mm_tag=options.mm_tag, cram_ref=options.ref_genome);
                    intron_list = intron_list[0] ### TODO

                ### extract mean exon coverage for all segments
                for j in range(gg.segmentgraph.segments.shape[1]):
                    idx = np.arange(gg.segmentgraph.segments[0, j], gg.segmentgraph.segments[1, j]) - gg.start
                    counts[f, i].segments[j] = np.mean(np.sum(tracks[:, idx], axis=0))
                    counts[f, i].seg_pos[j] = np.sum(np.sum(tracks[:, idx], axis=0) > 0)

                k, l = np.where(gg.segmentgraph.seg_edges == 1)

                ### there are no introns to count
                if intron_list.shape[0] == 0:
                    for m in range(k.shape[0]):
                        if counts[f, i].edges.shape[0] == 0:
                            counts[f, i].edges = np.atleast_2d(np.array([np.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), 0]))
                        else:
                            counts[f, i].edges = np.r_[counts[f, i].edges, np.atleast_2d(np.array([np.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), 0]))]
                    continue

                ### extract intron counts 
                for m in range(k.shape[0]):
                    idx = np.where((np.absolute(intron_list[:, 0] - gg.segmentgraph.segments[1, k[m]]) <= intron_tol) & (np.absolute(intron_list[:, 1] - gg.segmentgraph.segments[0, l[m]]) <= intron_tol))[0]
                    if counts[f, i].edges.shape[0] == 0:
                        if idx.shape[0] > 0:
                            counts[f, i].edges = np.atleast_2d(np.array([np.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), np.sum(intron_list[idx, 2])]))
                        else:
                            counts[f, i].edges = np.atleast_2d(np.array([np.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), 0]))
                    else:
                        if idx.shape[0] > 0:
                            counts[f, i].edges = np.r_[counts[f, i].edges, np.atleast_2d(np.array([np.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), np.sum(intron_list[idx, 2])]))]
                        else:
                            counts[f, i].edges = np.r_[counts[f, i].edges, np.atleast_2d(np.array([np.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), 0]))]

    if fn_out is not None:
        pickle.dump(counts, open(fn_out, 'wb'), -1)
    else:
        return counts


def count_graph_coverage_wrapper(fname_in, fname_out, bam_fnames, options, sample_idx=None, qmode='all'):

    (genes, inserted) = pickle.load(open(fname_in, 'rb')) 
    for g in genes:
        g.from_sparse()
    
    if genes[0].segmentgraph is None or genes[0].segmentgraph.is_empty():
        for g in genes:
            g.segmentgraph = Segmentgraph(g)
            g.to_sparse()
        pickle.dump((genes, inserted), open(fname_in, 'wb'), -1)
        for g in genes:
            g.from_sparse()

    counts = dict()
    counts['segments'] = []
    counts['seg_pos'] = []
    counts['gene_ids_segs'] = []
    counts['edges'] = []
    counts['gene_ids_edges'] = []
    counts['seg_len'] = np.hstack([x.segmentgraph.segments[1, :] - x.segmentgraph.segments[0, :] for x in genes]).T
    counts['gene_names'] = np.array([x.name for x in genes], dtype='str')

    if options.merge == 'single':
        print('\nprocessing %s' % (options.samples[sample_idx]))
        counts_tmp = count_graph_coverage(genes, bam_fnames[sample_idx], options)
    elif options.merge == 'merge_graphs' and qmode == 'single':
        print('\nquantifying merged graph in single mode (first file only) on %s' % options.samples[0])
        counts_tmp = count_graph_coverage(genes, bam_fnames[0], options)
    else:
        for s_idx in range(options.samples.shape[0]):
            print('\n%i/%i' % (s_idx + 1, options.samples.shape[0]))
            if s_idx == 0:
                counts_tmp = count_graph_coverage(genes, bam_fnames[s_idx], options)
            else:
                counts_tmp = np.r_[np.atleast_2d(counts_tmp), count_graph_coverage(genes, bam_fnames[s_idx], options)]

    for c in range(counts_tmp.shape[1]):
        counts['segments'].append(np.hstack([np.atleast_2d(x.segments).T for x in counts_tmp[:, c]]))
        counts['seg_pos'].append(np.hstack([np.atleast_2d(x.seg_pos).T for x in counts_tmp[:, c]]))
        counts['gene_ids_segs'].append(np.ones((np.atleast_2d(counts_tmp[0, c].seg_pos).shape[1], 1), dtype='int') * c)
        tmp = [np.atleast_2d(x.edges) for x in counts_tmp[:, c] if x.edges.shape[0] > 0]
        if len(tmp) == 0:
            continue
        tmp = np.hstack(tmp)
        if tmp.shape[0] > 0:
            counts['edges'].append(np.c_[tmp[:, 0], tmp[:, np.arange(1, tmp.shape[1], 2)]])
            counts['gene_ids_edges'].append(np.ones((tmp.shape[0], 1), dtype='int') * c)

    ### write result data to hdf5
    for key in counts:
        counts[key] = np.vstack(counts[key]) if len(counts[key]) > 0 else counts[key]
    counts['edge_idx'] = counts['edges'][:, 0] if len(counts['edges']) > 0 else np.array([])
    counts['edges'] = counts['edges'][:, 1:] if len(counts['edges']) > 0 else np.array([])
    h5fid = h5py.File(fname_out, 'w')
    h5fid.create_dataset(name='samples', data=codeUTF8(options.samples))
    h5fid['strains'] = h5py.SoftLink('/samples')
    for key in counts:
        if np.issubdtype(counts[key].dtype, np.str_):
            h5fid.create_dataset(name=key, data=codeUTF8(counts[key]))
        else:
            h5fid.create_dataset(name=key, data=counts[key])
    h5fid.close()

def collect_single_quantification_results(fname_out, sample_idxs, options):

        ### merge results from single count files
        if options.verbose:
            print('\nCollecting count data from files quantified in single mode ...\n')
            print('writing data to %s' % fname_out)

        ### write data to hdf5 continuously
        h5fid = h5py.File(fname_out, 'w')
        for i, idx in enumerate(sample_idxs):

            fname = get_filename('fn_collect_in', options, options.samples[idx])
            if options.verbose:
                print('collecting counts from %s (%i/%i)' % (fname, i+1, len(sample_idxs)))
            if not os.path.exists(fname):
                print('ERROR: Counts for sample %i cannot be loaded - %s is not present' % (i + 1, fname), file=sys.stderr)
                sys.exit(1)

            CIN = h5py.File(fname, 'r')
            if i == 0:
                for k in ['gene_ids_segs', 'edge_idx', 'gene_ids_edges', 'gene_names', 'seg_len']:
                    h5fid.create_dataset(name=k, data=CIN[k][:], chunks=True, compression='gzip')
                h5fid.create_dataset(name='edges', data=CIN['edges'][:], chunks=True, compression='gzip', maxshape=(CIN['edges'].shape[0], None))
                h5fid.create_dataset(name='segments', data=CIN['segments'][:], chunks=True, compression='gzip', maxshape=(CIN['segments'].shape[0], None))
                h5fid.create_dataset(name='seg_pos', data=CIN['seg_pos'][:], chunks=True, compression='gzip', maxshape=(CIN['seg_pos'].shape[0], None))
                h5fid.create_dataset(name='samples', data=CIN['samples'][:], dtype='|S255', chunks=True, compression='gzip', maxshape=(None,))
            else:
                appendToHDF5(h5fid, CIN['edges'][:], 'edges', faxis=1, daxis=1)
                appendToHDF5(h5fid, CIN['segments'][:], 'segments', faxis=1, daxis=1)
                appendToHDF5(h5fid, CIN['seg_pos'][:], 'seg_pos', faxis=1, daxis=1)
                appendToHDF5(h5fid, CIN['samples'][:], 'samples', faxis=0, daxis=0)
        h5fid.close() 

def compute_gene_expression(options, fname_genes, fname_count_in, fn_out=None, sample_idx=None):

    if options.verbose:
        sys.stdout.write('Quantifying gene expression ...\n')

    ### load gene information
    genes = pickle.load(open(fname_genes, 'rb'), encoding='latin1')[0]
    numgenes = genes.shape[0]

    ### open hdf5 file containing graph count information
    IN = h5py.File(fname_count_in, 'r')
    samples = IN['samples'][:].astype('str')
    ### sort by sample ID
    if sample_idx is None:
        sample_idx = np.arange(samples.shape[0])
    gene_counts = np.zeros((numgenes, len(sample_idx)), dtype='float')
    gene_counts_non_alt = np.zeros((numgenes, len(sample_idx)), dtype='float')
    gene_ids = np.array([x.name for x in genes], dtype='str')
    gene_symbols = np.array([x.symbol if (not x is None) and hasattr(x, 'symbol') else 'NA' for x in genes], dtype='str')

    seg_lens = IN['seg_len'][:]
    gene_ids_segs = IN['gene_ids_segs'][:].astype('int')

    ### no longer assume that the gene_ids_segs are sorted by gene ID
    s_idx = np.argsort(gene_ids_segs[:, 0], kind='mergesort')
    _, u_idx = np.unique(gene_ids_segs[s_idx, 0], return_index=True)
    s_idx = s_idx[u_idx]

    ### buffer seg_counts
    cs = IN['segments'].chunks[0] if not IN['segments'].chunks is None else 1  # chunk size
    co = 0 # chunk offset
    cn = 10000 // cs # chunk number
    seg_buffer = IN['segments'][co:(co+(cn*cs)), :]

    ### iterate over genes
    for gidx, iidx in enumerate(s_idx):

        if options.verbose:
            log_progress(gidx, numgenes, 100)

        seg_idx = np.arange(iidx, iidx + genes[gidx].segmentgraph.seg_edges.shape[0])
        ### do we need a new chunk?
        if seg_idx.max() >= co+(cn*cs):
            co = seg_idx.min()
            seg_buffer = IN['segments'][co:(co+(cn*cs)), :]

        gene_idx = gene_ids_segs[seg_idx, 0]
        if len(gene_idx.shape) > 0:
            gene_idx = gene_idx[0]

        assert(decodeUTF8(IN['gene_names'][:][gene_idx]) == genes[gidx].name)
        assert(genes[gidx].name == gene_ids[gidx])

        if seg_idx.shape[0] > 1:
            gene_counts[gidx, :] = np.squeeze(np.dot(seg_buffer[seg_idx - co, :][:, sample_idx].T, seg_lens[seg_idx])) / options.readlen
        else:
            gene_counts[gidx, :] = seg_buffer[seg_idx[0] - co, :][sample_idx] * seg_lens[seg_idx] / options.readlen

        ### get idx of non alternative segments
        non_alt_idx = genes[gidx].get_non_alt_seg_ids()
        seg_idx = seg_idx[non_alt_idx]

        if seg_idx.shape[0] > 1:
            gene_counts_non_alt[gidx, :] = np.squeeze(np.dot(seg_buffer[seg_idx - co, :][:, sample_idx].T, seg_lens[seg_idx])) / options.readlen
        else:
            gene_counts_non_alt[gidx, :] = seg_buffer[seg_idx[0] - co, :][sample_idx] * seg_lens[seg_idx] / options.readlen

    IN.close()

    if options.verbose:
        sys.stdout.write('\n... done.\n')

    ### compute size factors
    sf = get_size_factors(gene_counts, options)
    sf_non_alt = get_size_factors(gene_counts_non_alt, options)

    ### write results to hdf5
    if fn_out is not None:
        OUT = h5py.File(fn_out, 'w')
        OUT.create_dataset(name='samples', data=codeUTF8(samples[sample_idx]), compression='gzip')
        OUT['strains'] = h5py.SoftLink('/samples')
        OUT.create_dataset(name='gene_ids', data=codeUTF8(gene_ids), compression='gzip')
        OUT.create_dataset(name='gene_symbols', data=codeUTF8(gene_symbols), compression='gzip')
        OUT.create_dataset(name='raw_count', data=gene_counts, compression='gzip', chunks=True)
        OUT.create_dataset(name='raw_count_non_alt', data=gene_counts_non_alt, compression='gzip', chunks=True)
        OUT.create_dataset(name='size_factors', data=sf, compression='gzip')
        OUT.create_dataset(name='size_factors_non_alt', data=sf_non_alt, compression='gzip')
        OUT.close()


def get_size_factors(gene_counts, options, kind='geomean'):

    if options.verbose:
        print('Estimating size factors ...')

    size_factors = []

    ### take geometric mean of counts
    if kind == 'geomean':
        gmean = np.exp(np.mean(np.log(gene_counts + 1), axis=1))
        for i in range(gene_counts.shape[1]):
            idx = gene_counts[:, i] > 0
            size_factors.append(np.median(gene_counts[idx, i] / gmean[idx]))
    elif kind == 'tc': # total count
        for i in range(gene_counts.shape[1]):
            size_factors.append(gene_counts[:, i].sum() / 100000000) 
    elif kind == 'uq': # upper quartile
        for i in range(gene_counts.shape[1]):
            size_factors.append(scoreatpercentile(gene_counts[:, i], 75) + 1) 
    else:
        sys.stderr.write('ERROR: unknown normalization method %s\n' % kind)
        sys.exit(1)

    if options.verbose:
        sys.stdout.write('\n... done.\n')

    return np.array(size_factors, dtype='float')

