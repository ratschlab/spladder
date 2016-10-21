""" This module is concerned with functions that are specific to our identity (matlab vs python)"""

import copy
import cPickle
import os
import h5py
import scipy as sp
import re

### return identity
def identity():
    return "python"

### return a gene object
def get_gene(gene):
    """This function returns a deep copy of the gene object"""

    return copy.deepcopy(gene)

### load the gene list in the right format
def load_genes(options):
    """This is a helper function to load the gene data from file"""

    if options.validate_sg:
        gene_fname = os.path.join(options.outdir, 'spladder', 'genes_graph_conf%s.merge_graphs.validated.pickle' % options.confidence)
    else:
        gene_fname = os.path.join(options.outdir, 'spladder', 'genes_graph_conf%s.merge_graphs.pickle' % options.confidence)

    if options.verbose:
        print 'loading annotation information from %s' % gene_fname
    (genes, events) = cPickle.load(open(gene_fname, 'r'))

    return genes


### load a single event from file
def load_events(options, event_info):

    event_list = []
    for event_type in sp.unique(event_info[:, 0]):
        events = cPickle.load(open(os.path.join(options.outdir, 'merge_graphs_%s_C%s.pickle' % (event_type, options.confidence)),'r'))
        s_idx = sp.where(event_info[:, 0] == event_type)[0]
        for e in s_idx:
            event_list.append(events[int(event_info[e, 1])])

    return event_list

def get_gene_ids(options):

    gids = []

    if options.event_id is None:
        for event_type in options.event_types:
            IN = h5py.File(os.path.join(options.outdir, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, options.confidence)), 'r')
            if 'conf_idx' in IN and IN['conf_idx'].shape[0] > 0:
                gids.extend(IN['gene_idx'][:][IN['conf_idx'][:]])
            IN.close()
    else:
        IN = h5py.File(os.path.join(options.outdir, 'merge_graphs_%s_C%i.counts.hdf5' % (re.sub(r'_[0-9]*$', '', options.event_id), options.confidence)), 'r')
        gids.append(IN['gene_idx'][int(options.event_id.split('_')[-1])].astype('int'))
        IN.close()

    return sp.unique(gids)

def get_conf_events(options, gid):

    event_info = []

    for event_type in options.event_types:
	if event_type not in ["exon_skip","intron_retention","alt_3prime","alt_5prime","mult_exon_skip","mutex_exons"]:
		raise Exception('Unknown event type: %s' % event_type)
        IN = h5py.File(os.path.join(options.outdir, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, options.confidence)), 'r')
        if 'conf_idx' in IN and IN['conf_idx'].shape[0] > 0:
            conf_idx = IN['conf_idx'][:]
            k_idx = sp.where(sp.in1d(IN['gene_idx'][:][conf_idx], gid))[0]
            if k_idx.shape[0] > 0:
                event_info.extend([[event_type, x] for x in conf_idx[k_idx]])
        IN.close()

    return sp.array(event_info, dtype='str')

def get_seg_counts(options, gid):

    IN = h5py.File(os.path.join(options.outdir, 'spladder', 'genes_graph_conf%i.merge_graphs.count.hdf5' % (options.confidence)), 'r')
    idx = sp.where(IN['gene_ids_edges'][:] == gid)[0]
    edges = IN['edges'][:][idx, :]
    edge_idx = IN['edge_idx'][:][idx]
    idx = sp.where(IN['gene_ids_segs'][:] == gid)[0]
    segments = IN['segments'][:][idx]
    IN.close()

    return (segments, edges, edge_idx)

