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
def load_genes(CFG):
    """This is a helper function to load the gene data from file"""

    if CFG['validate_splicegraphs']:
        gene_fname = os.path.join(CFG['out_dirname'], 'spladder', 'genes_graph_conf%s.merge_graphs.validated.pickle' % CFG['confidence_level'])
    else:
        gene_fname = os.path.join(CFG['out_dirname'], 'spladder', 'genes_graph_conf%s.merge_graphs.pickle' % CFG['confidence_level'])

    if CFG['verbose']:
        print 'loading annotation information from %s' % gene_fname
    (genes, events) = cPickle.load(open(gene_fname, 'r'))

    return genes


### load a single event from file
def load_events(CFG, event_info):

    event_list = []
    for event_type in sp.unique(event_info[:, 0]):
        events = cPickle.load(open(os.path.join(CFG['out_dirname'], 'merge_graphs_%s_C%s.pickle' % (event_type, CFG['confidence_level'])),'r'))
        s_idx = sp.where(event_info[:, 0] == event_type)[0]
        for e in s_idx:
            event_list.append(events[int(event_info[e, 1])])

    return event_list

def get_gene_ids(CFG, gene_names=None):

    gids = []
    
    if gene_names is not None:
        for g in gene_names:
            event_type = re.sub(r'_[0-9]*$', '', g[1])
            IN = h5py.File(os.path.join(CFG['out_dirname'], 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, CFG['confidence_level'])), 'r')
            gids.append([sp.where(IN['gene_names'][:] == g[0])[0][0], g[1]])
            IN.close()
    elif CFG['event_id'] is None:
        for event_type in CFG['event_types']:
            IN = h5py.File(os.path.join(CFG['out_dirname'], 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, CFG['confidence_level'])), 'r')
            if 'conf_idx' in IN and IN['conf_idx'].shape[0] > 0:
                gids.append([IN['gene_idx'][:][IN['conf_idx'][:]], None])
            IN.close()
    else:
        IN = h5py.File(os.path.join(CFG['out_dirname'], 'merge_graphs_%s_C%i.counts.hdf5' % (re.sub(r'_[0-9]*$', '', CFG['event_id']), CFG['confidence_level'])), 'r')
        gids.append([IN['gene_idx'][int(CFG['event_id'].split('_')[-1])].astype('int'), CFG['event_id']])
        IN.close()

    return gids

def get_conf_events(CFG, gid):

    event_info = []

    for event_type in CFG['event_types']:
        if event_type not in ["exon_skip","intron_retention","alt_3prime","alt_5prime","mult_exon_skip","mutex_exons"]:
            raise Exception('Unknown event type: %s' % event_type)
        IN = h5py.File(os.path.join(CFG['out_dirname'], 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, CFG['confidence_level'])), 'r')
        if 'conf_idx' in IN and IN['conf_idx'].shape[0] > 0:
            conf_idx = IN['conf_idx'][:]
            k_idx = sp.where(sp.in1d(IN['gene_idx'][:][conf_idx], gid))[0]
            if k_idx.shape[0] > 0:
                event_info.extend([[event_type, x] for x in conf_idx[k_idx]])
        IN.close()

    return sp.array(event_info, dtype='str')

def get_seg_counts(CFG, gid):

    IN = h5py.File(os.path.join(CFG['out_dirname'], 'spladder', 'genes_graph_conf%i.merge_graphs.count.hdf5' % (CFG['confidence_level'])), 'r')
    idx = sp.where(IN['gene_ids_edges'][:] == gid)[0]
    edges = IN['edges'][:][idx, :]
    edge_idx = IN['edge_idx'][:][idx]
    idx = sp.where(IN['gene_ids_segs'][:] == gid)[0]
    segments = IN['segments'][:][idx]
    strains = IN['strains'][:]
    IN.close()

    return (segments, edges, edge_idx, strains)

