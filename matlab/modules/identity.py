""" This module is concerned with functions that are specific to our identity (matlab vs python)"""

import os
import scipy.io as scio
import scipy as sp
import h5py

from .classes.gene import Gene
from .classes.event import Event

### return identity
def identity():
    return "matlab"

### return a gene object
def get_gene(matlab_gene):
    """This function takes a gene in matlab format and returns 
    the corresponding python object"""

    gene = Gene()
    gene.splicegraph.vertices = matlab_gene.splicegraph[0, 0]
    gene.splicegraph.edges = matlab_gene.splicegraph[0, 1]
    gene.segmentgraph.segments = matlab_gene.segmentgraph[0, 0]
    gene.segmentgraph.seg_edges = matlab_gene.segmentgraph[0, 2]
    gene.exons = [matlab_gene.exons[0, x] for x in range(matlab_gene.exons.shape[1])]
    gene.chr = str(matlab_gene.chr[0])
    gene.name = matlab_gene.name[0]

    return gene


### load the gene list in the right format
def load_genes(options):
    """This is a helper function to load the gene data from file"""

    if options.validate_sg:
        genes = scio.loadmat(os.path.join(options.outdir, 'spladder', 'genes_graph_conf%s.merge_graphs.validated.mat' % options.confidence), struct_as_record=False)['genes'][0, :]
    else:
        genes = scio.loadmat(os.path.join(options.outdir, 'spladder', 'genes_graph_conf%s.merge_graphs.mat' % options.confidence), struct_as_record=False)['genes'][0, :]

    return genes

### load a single event from file
def load_events(options, event_info):

    event_list = []
    for event_type in sp.unique(event_info[:, 0]):
        events_matlab = scio.loadmat(os.path.join(options.outdir, 'merge_graphs_%s_C%s.mat' % (event_type, options.confidence)), struct_as_record=False)['events_all']
        s_idx = sp.where(event_info[:, 0] == event_type)[0]
        for e in s_idx:
            event_matlab = events_matlab[0, int(event_info[e, 1])]
            event = Event(event_matlab.event_type[0])
            if event.event_type == 'exon_skip':
                event.exons1 = sp.r_[event_matlab.exon_pre, event_matlab.exon_aft]
                event.exons2 = sp.r_[event_matlab.exon_pre, event_matlab.exon, event_matlab.exon_aft]
            elif event.event_type == 'intron_retention':
                event.exons1 = sp.r_[event_matlab.exon1, event_matlab.exon2]
                event.exons2 = sp.array([event_matlab.exon1[0, 0], event_matlab.exon2[0, 1]])
            elif event.event_type in ['alt_3prime', 'alt_5prime']:
                event.exons1 = sp.r_[event_matlab.exon_const, event_matlab.exon_alt1]
                s_idx = sp.argsort(event.exons1[:, 0])
                event.exons1 = event.exons1[s_idx, :]
                event.exons2 = sp.r_[event_matlab.exon_const, event_matlab.exon_alt2]
                s_idx = sp.argsort(event.exons2[:, 0])
                event.exons2 = event.exons2[s_idx, :]
            elif event.event_type[0] == 'mutex_exons':
                event.exons1 = sp.r_[event_matlab.exon_pre, event_matlab.exon1, event_matlab.exon_aft]
                event.exons2 = sp.r_[event_matlab.exon_pre, event_matlab.exon2, event_matlab.exon_aft]
            elif event.event_type[0] == 'mult_exon_skip':
                event.exons1 = sp.r_[event_matlab.exon_pre, event_matlab.exon_aft]
                event.exons2 = sp.r_[event_matlab.exon_pre, event_matlab.exons.reshape(event_matlab.exons.shape[1] / 2, 2), event_matlab.exon_aft]
            event_list.append(event)

    return event_list

def get_gene_ids(options):

    gids = []

    for event_type in options.event_types:
        IN = h5py.File(os.path.join(options.outdir, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, options.confidence)), 'r')
        if 'conf_idx' in IN and IN['conf_idx'].shape[0] > 0 and IN['conf_idx'][0] != -1:
            gids.extend(IN['gene_idx'][:][IN['conf_idx'][:].astype('int') - 1].astype('int') - 1)
        IN.close()

    return sp.unique(gids)

def get_conf_events(options, gid):

    event_info = []

    for event_type in options.event_types:
        IN = h5py.File(os.path.join(options.outdir, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, options.confidence)), 'r')
        if 'conf_idx' in IN and IN['conf_idx'].shape[0] > 0 and IN['conf_idx'][0] != -1:
            conf_idx = IN['conf_idx'][:].astype('int') - 1
            k_idx = sp.where(sp.in1d(IN['gene_idx'][:][conf_idx].astype('int') - 1, gid))[0]
            if k_idx.shape[0] > 0:
                event_info.extend([[event_type, x] for x in conf_idx[k_idx]])
        IN.close()

    return sp.array(event_info, dtype='str')

def get_seg_counts(options, gid):

    IN = h5py.File(os.path.join(options.outdir, 'spladder', 'genes_graph_conf%i.merge_graphs.count.mat' % (options.confidence)), 'r')
    idx = sp.where(IN['gene_ids_edges'][:] == gid + 1)[1]
    edges = IN['edges'][:][:, idx].T
    edge_idx = IN['edge_idx'][:][0, idx].astype('int') - 1
    idx = sp.where(IN['gene_ids_segs'][:] == gid + 1)[1]
    segments = IN['segments'][:][:, idx].T
    IN.close()

    return (segments, edges, edge_idx)

