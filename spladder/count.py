#if __package__ is None:
#    __package__ ='modules'

import pickle
import math
import h5py
import scipy as sp
import re

from .classes.segmentgraph import Segmentgraph
from .classes.counts import Counts
from .helpers import *
from .reads import *
from .hdf5 import appendToHDF5
from . import rproc as rp

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
    counts = sp.zeros((len(fn_bam), genes.shape[0]), dtype='object')

    intron_tol = 0 

    sys.stdout.write('genes: %i\n' % genes.shape[0])
    for f in range(counts.shape[0]):
        sys.stdout.write('\nsample %i/%i\n' % (f + 1, counts.shape[0])) 

        ### iterate over all genes and generate counts for
        ### the segments in the segment graph
        ### and the splice junctions in the splice graph
        ### iterate per contig, so the bam caching works better
        contigs = sp.array([x.chr for x in genes])
        for contig in sp.unique(contigs):
            contig_idx = sp.where(contigs == contig)[0]
            bam_cache = dict()
            print('\ncounting %i genes on contig %s' % (contig_idx.shape[0], contig))
            for ii,i in enumerate(contig_idx):
                sys.stdout.write('.')
                if ii > 0 and ii % 50 == 0:
                    sys.stdout.write('%i/%i\n' % (ii, contig_idx.shape[0]))
                sys.stdout.flush()
                gg = genes[i]
                if gg.segmentgraph.is_empty():
                    gg.segmentgraph = Segmentgraph(gg)
                gg.start = gg.segmentgraph.segments.ravel().min()
                gg.stop = gg.segmentgraph.segments.ravel().max()

                counts[f, i] = Counts(gg.segmentgraph.segments.shape[1])

                if options.sparse_bam and \
                  (fn_bam[f].endswith('npz') or \
                   os.path.exists(re.sub(r'bam$', '', fn_bam[f]) + 'npz') or \
                   fn_bam[f].endswith('hdf5') or \
                   os.path.exists(re.sub(r'bam$', '', fn_bam[f]) + 'hdf5')):
                    ### make sure that we query the right contig from cache
                    assert(gg.chr == contig)
                    (tracks, intron_list) = add_reads_from_sparse_bam(gg, fn_bam[f], contig, options.confidence, types=['exon_track','intron_list'], filter=None, cache=bam_cache)
                else:
                    ### add RNA-seq evidence to the gene structure
                    (tracks, intron_list) = add_reads_from_bam(gg, fn_bam[f], ['exon_track','intron_list'], None, options.var_aware, options.primary_only);
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
                        if counts[f, i].edges.shape[0] == 0:
                            counts[f, i].edges = sp.atleast_2d(sp.array([sp.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), 0]))
                        else:
                            counts[f, i].edges = sp.r_[counts[f, i].edges, sp.atleast_2d(sp.array([sp.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), 0]))]
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
        pickle.dump(counts, open(fn_out, 'wb'), -1)
    else:
        return counts



def count_graph_coverage_wrapper(fname_in, fname_out, options, sample_idx=None, qmode='all'):

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
    counts['seg_len'] = sp.hstack([x.segmentgraph.segments[1, :] - x.segmentgraph.segments[0, :] for x in genes]).T
    counts['gene_names'] = sp.array([x.name for x in genes], dtype='str')

    if not options.pyproc:
        if options.merge == 'single':
            print('\nprocessing %s' % (options.samples[sample_idx]))
            counts_tmp = count_graph_coverage(genes, options.bam_fnames[sample_idx], options)
        elif options.merge == 'merge_graphs' and qmode == 'single':
            print('\nquantifying merged graph in single mode (first file only) on %s' % options.samples[0])
            counts_tmp = count_graph_coverage(genes, options.bam_fnames[0], options)
        else:
            for s_idx in range(options.strains.shape[0]):
                print('\n%i/%i' % (s_idx + 1, options.strains.shape[0]))
                if s_idx == 0:
                    counts_tmp = count_graph_coverage(genes, options.bam_fnames[s_idx], options)
                else:
                    counts_tmp = sp.r_[sp.atleast_2d(counts_tmp), count_graph_coverage(genes, options.bam_fnames[s_idx], options)]

        for c in range(counts_tmp.shape[1]):
            counts['segments'].append(sp.hstack([sp.atleast_2d(x.segments).T for x in counts_tmp[:, c]]))
            counts['seg_pos'].append(sp.hstack([sp.atleast_2d(x.seg_pos).T for x in counts_tmp[:, c]]))
            counts['gene_ids_segs'].append(sp.ones((sp.atleast_2d(counts_tmp[0, c].seg_pos).shape[1], 1), dtype='int') * c)
            tmp = [sp.atleast_2d(x.edges) for x in counts_tmp[:, c] if x.edges.shape[0] > 0]
            if len(tmp) == 0:
                continue
            tmp = sp.hstack(tmp)
            if tmp.shape[0] > 0:
                counts['edges'].append(sp.c_[tmp[:, 0], tmp[:, sp.arange(1, tmp.shape[1], 2)]])
                counts['gene_ids_edges'].append(sp.ones((tmp.shape[0], 1), dtype='int') * c)

        ### write result data to hdf5
        for key in counts:
            counts[key] = sp.vstack(counts[key]) if len(counts[key]) > 0 else counts[key]
        counts['edge_idx'] = counts['edges'][:, 0] if len(counts['edges']) > 0 else sp.array([])
        counts['edges'] = counts['edges'][:, 1:] if len(counts['edges']) > 0 else sp.array([])
        h5fid = h5py.File(fname_out, 'w')
        h5fid.create_dataset(name='strains', data=codeUTF8(options.strains))
        for key in counts:
            if sp.issubdtype(counts[key].dtype, sp.str_):
                h5fid.create_dataset(name=key, data=codeUTF8(counts[key]))
            else:
                h5fid.create_dataset(name=key, data=counts[key])
        h5fid.close()
    else:
        ### have an adaptive chunk size, that takes into account the number of strains (take as many genes as it takes to have ~10K strains)
        if options.sparse_bam:
            chunksize = int(max(1, math.floor(1000000 / len(options.strains))))
        else:
            chunksize = int(max(1, math.floor(100000 / len(options.strains))))

        jobinfo = []

        PAR = dict()
        PAR['options'] = options
        if options.merge == 'single':
            PAR['options'].bam_fnames = PAR['options'].bam_fnames[sample_idx]
            PAR['options'].samples = PAR['options'].samples[sample_idx]
            PAR['options'].strains = PAR['options'].strains[sample_idx]

        #s_idx = sp.argsort([x.chr for x in genes]) # TODO
        s_idx = sp.arange(genes.shape[0])
        for c_idx in range(0, s_idx.shape[0], chunksize):
            cc_idx = min(s_idx.shape[0], c_idx + chunksize)
            fn = re.sub(r'.hdf5$', '', fname_out) + '.chunk_%i_%i.pickle' % (c_idx, cc_idx)
            if os.path.exists(fn):
                continue
            else:
                print('submitting chunk %i to %i (%i)' % (c_idx, cc_idx, s_idx.shape[0]))
                PAR['genes'] = genes[s_idx][c_idx:cc_idx]
                for gg in PAR['genes']:
                    gg.to_sparse()
                PAR['fn_bam'] = options.bam_fnames
                PAR['fn_out'] = fn
                PAR['options'] = options
                jobinfo.append(rp.rproc('count_graph_coverage', PAR, 15000, options.options_rproc, 60*48))

        rp.rproc_wait(jobinfo, 30, 1.0, -1)
        del genes

        ### merge results from count chunks
        if options.verbose:
            print('\nCollecting count data from chunks ...\n')
            print('writing data to %s' % fname_out)

        ### write data to hdf5 continuously
        h5fid = h5py.File(fname_out, 'w')
        h5fid.create_dataset(name='gene_names', data=codeUTF8(counts['gene_names']))
        h5fid.create_dataset(name='seg_len', data=counts['seg_len'])
        h5fid.create_dataset(name='strains', data=codeUTF8(options.strains))
        for c_idx in range(0, s_idx.shape[0], chunksize):
            cc_idx = min(s_idx.shape[0], c_idx + chunksize)
            if options.verbose:
                print('collecting chunk %i-%i (%i)' % (c_idx, cc_idx, s_idx.shape[0]))
            fn = re.sub(r'.hdf5$', '', fname_out) + '.chunk_%i_%i.pickle' % (c_idx, cc_idx)
            if not os.path.exists(fn):
                print('ERROR: Not all chunks in counting graph coverage completed!', file=sys.stderr)
                sys.exit(1)
            else:
                counts_tmp = pickle.load(open(fn, 'rb'))
                for c in range(counts_tmp.shape[1]):
                    if 'segments' in h5fid:
                        appendToHDF5(h5fid, sp.hstack([sp.atleast_2d(x.segments).T for x in counts_tmp[:, c]]), 'segments')
                        appendToHDF5(h5fid, sp.hstack([sp.atleast_2d(x.seg_pos).T for x in counts_tmp[:, c]]), 'seg_pos') 
                        appendToHDF5(h5fid, sp.ones((sp.atleast_2d(counts_tmp[0, c].seg_pos).shape[1], 1), dtype='int') * (s_idx[c_idx + c]), 'gene_ids_segs')
                    else:
                        h5fid.create_dataset(name='segments', data=sp.hstack([sp.atleast_2d(x.segments).T for x in counts_tmp[:, c]]), chunks=True, compression='gzip', maxshape=(None, len(options.strains)))
                        h5fid.create_dataset(name='seg_pos', data=sp.hstack([sp.atleast_2d(x.seg_pos).T for x in counts_tmp[:, c]]), chunks=True, compression='gzip', maxshape=(None, len(options.strains)))
                        h5fid.create_dataset(name='gene_ids_segs', data=sp.ones((sp.atleast_2d(counts_tmp[0, c].seg_pos).shape[1], 1), dtype='int') * (s_idx[c_idx + c]), chunks=True, compression='gzip', maxshape=(None, 1))

                    tmp = [sp.atleast_2d(x.edges) for x in counts_tmp[:, c] if x.edges.shape[0] > 0]
                    if len(tmp) == 0:
                        continue
                    tmp = sp.hstack(tmp)
                    if tmp.shape[0] > 0:
                        if 'edges' in h5fid:
                            appendToHDF5(h5fid, tmp[:, sp.arange(1, tmp.shape[1], 2)], 'edges')
                            appendToHDF5(h5fid, tmp[:, 0], 'edge_idx')
                            appendToHDF5(h5fid, sp.ones((tmp.shape[0], 1), dtype='int') * (s_idx[c_idx + c]), 'gene_ids_edges')
                        else:
                            h5fid.create_dataset(name='edges', data=tmp[:, sp.arange(1, tmp.shape[1], 2)], chunks=True, compression='gzip', maxshape=(None, tmp.shape[1] / 2))
                            h5fid.create_dataset(name='edge_idx', data=tmp[:, 0], chunks=True, compression='gzip', maxshape=(None,))
                            h5fid.create_dataset(name='gene_ids_edges', data=sp.ones((tmp.shape[0], 1), dtype='int') * (s_idx[c_idx + c]), chunks=True, compression='gzip', maxshape=(None, 1))
                del tmp, counts_tmp
        h5fid.close()


def collect_single_quantification_results(fname_out, sample_idxs, options):

        ### merge results from single count files
        if options.verbose:
            print('\nCollecting count data from files quantified in single mode ...\n')
            print('writing data to %s' % fname_out)

        ### write data to hdf5 continuously
        h5fid = h5py.File(fname_out, 'w')
        for i, idx in enumerate(sample_idxs):

            fname = get_filename('fn_count_in', options, sample_idx=idx)
            if options.verbose:
                print('collecting counts from %s (%i/%i)' % (fname, i+1, len(sample_idxs)))
            if not os.path.exists(fname):
                print('ERROR: Counts for sample %i cannot be loaded - %s is not present' % (i + 1, fname), file=sys.stderr)
                sys.exit(1)

            CIN = h5py.File(fname, 'r')
            if i == 0:
                for k in ['segments', 'seg_pos', 'gene_ids_segs', 'edge_idx', 'gene_ids_edges', 'gene_names', 'seg_len']:
                    h5fid.create_dataset(name=k, data=CIN[k][:], chunks=True, compression='gzip')
                h5fid.create_dataset(name='edges', data=CIN['edges'][:], chunks=True, compression='gzip', maxshape=(CIN['edges'].shape[0], None))
                h5fid.create_dataset(name='strains', data=CIN['strains'][:], chunks=True, compression='gzip', maxshape=(None,))
            else:
                appendToHDF5(h5fid, CIN['edges'][:], 'edges', faxis=1, daxis=1)
                appendToHDF5(h5fid, CIN['strains'][:], 'strains', faxis=0, daxis=0)
        h5fid.close() 

