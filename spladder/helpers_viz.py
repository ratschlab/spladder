""" This module is concerned with functions that are specific to our identity (matlab vs python)"""

import copy
import pickle
import os
import h5py
import scipy as sp
import re
import sys

from . import helpers as hp

### return a gene object
def get_gene(gene):
    """This function returns a deep copy of the gene object"""

    return copy.deepcopy(gene)


def _get_gene_fname(outdir, confidence, validate_sg):

    val_tag = ''
    if validate_sg:
        val_tag = '.validated'
    return os.path.join(outdir, 'spladder', 'genes_graph_conf%s.merge_graphs%s.pickle' % (confidence, val_tag))


def get_gene_names(outdir, confidence, validate_sg, verbose):
    """Return the list of gene names used in this run"""
    
    gene_file = _get_gene_fname(outdir, confidence, validate_sg)
    gene_name_file = re.sub(r'.pickle$', '', gene_file) + '.names.pickle'
    if not os.path.exists(gene_name_file):
        if verbose:
            print('Generating list of gene names for easy access')
        tmp_genes = load_genes(outdir, confidence, validate_sg, verbose)
        gene_names = sp.array([x.name.split('.')[0] for x in tmp_genes])
        pickle.dump(gene_names, open(gene_name_file, 'wb'), -1)
    else:
        gene_names = pickle.load(open(gene_name_file, 'rb'), encoding='latin1')

    return gene_names


### load the gene list in the right format
def load_genes(outdir, confidence, validate_sg, verbose, idx=None, genes=None):
    """This is a helper function to load the gene data from file"""

    if not genes is None:
        if not idx is None:
            return copy.deepcopy(genes[idx])
    else:
        gene_file = _get_gene_fname(outdir, confidence, validate_sg)

        if verbose:
            print('loading annotation information from %s' % gene_file)
        if idx is None:
            (genes, events) = pickle.load(open(gene_file, 'rb'), encoding='latin1')
        else:
            gene_db_file = re.sub(r'.pickle$', '', gene_file) + '.db.pickle'
            gene_idx_file = re.sub(r'.pickle$', '', gene_file) + '.idx.pickle'
            if os.path.exists(gene_idx_file):
                genes = []
                offsets = pickle.load(open(gene_idx_file, 'rb'))
                gene_handle = open(gene_db_file, 'rb')
                if not hasattr(idx, '__iter__'):
                    idx = [idx]
                for e in idx:
                    gene_handle.seek(offsets[e], 0)
                    genes.append(pickle.load(gene_handle), encoding='latin1')
                genes = sp.array(genes)
            else:
                (genes, events) = pickle.load(open(gene_file, 'rb'), encoding='latin1')
                genes = genes[idx]

    return genes


### load a single event from file
def load_events(event_info, outdir, confidence, verbose):

    event_list = [] 
    for event_type in sp.unique(event_info[:, 0]):
        event_file = os.path.join(outdir, 'merge_graphs_%s_C%s.pickle' % (event_type, confidence))
        event_db_file = re.sub(r'.pickle$', '', event_file) + '.db.pickle'
        event_idx_file = re.sub(r'.pickle$', '', event_file) + '.idx.pickle'
        s_idx = sp.where(event_info[:, 0] == event_type)[0]
        if not os.path.exists(event_db_file):
            if verbose:
                print('Indexing event files for faster future access')
            events = pickle.load(open(os.path.join(outdir, 'merge_graphs_%s_C%s.pickle' % (event_type, confidence)),'rb'), encoding='latin1')
            out_db = open(event_db_file, 'wb')
            out_idx = open(event_idx_file, 'wb')
            offsets = []
            for o, obj in enumerate(events):
                if o > 0 and o % 1000 == 0:
                    sys.stdout.write('.')
                    if o % 10000 == 0:
                        sys.stdout.write('%i\n' % o)
                    sys.stdout.flush()
                offsets.append(out_db.tell())
                pickle.dump(obj, out_db, -1)
            pickle.dump(offsets, out_idx, -1)
            out_db.close()
            out_idx.close()
            for e in s_idx:
                event_list.append(events[int(event_info[e, 1])])
        else:
            offsets = pickle.load(open(event_idx_file, 'rb'))
            events_handle = open(event_db_file, 'rb')
            for e in s_idx:
                events_handle.seek(offsets[int(event_info[e, 1])], 0)
                event_list.append(pickle.load(events_handle))

    return event_list

def get_event_ids_from_gene(gene_id, event_type, outdir, confidence):
    
    eids = []
    IN = h5py.File(os.path.join(outdir, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, confidence)), 'r')
    if 'conf_idx' in IN and IN['conf_idx'].shape[0] > 0:
        cidx = IN['conf_idx'][:]
        gidx = sp.where(IN['gene_idx'][:] == gene_id)[0]
        eids.extend(sp.intersect1d(cidx, gidx))
    IN.close()

    return eids


def get_gene_ids(options, gene_names=None):

    gids = []
    
    if gene_names is not None:
        for g in gene_names:
            event_type = re.sub(r'_[0-9]*$', '', g[1])
            IN = h5py.File(os.path.join(options.outdir, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, options.confidence)), 'r')
            gids.append([sp.where(IN['gene_names'][:] == g[0])[0][0], g[1]])
            IN.close()
    elif options.event_id is None:
        for event_type in options.event_types:
            IN = h5py.File(os.path.join(options.outdir, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, options.confidence)), 'r')
            if 'conf_idx' in IN and IN['conf_idx'].shape[0] > 0:
                gids.append([IN['gene_idx'][:][IN['conf_idx'][:]], None])
            IN.close()
    else:
        IN = h5py.File(os.path.join(options.outdir, 'merge_graphs_%s_C%i.counts.hdf5' % (re.sub(r'_[0-9]*$', '', options.event_id), options.confidence)), 'r')
        gids.append([IN['gene_idx'][int(options.event_id.split('_')[-1])].astype('int'), options.event_id])
        IN.close()

    return gids

def get_conf_events(options, gid):

    event_info = []

    for event_type in options.event_types:
        if event_type not in ["exon_skip","intron_retention","alt_3prime","alt_5prime","mult_exon_skip","mutex_exons"]:
            raise Exception('Unknown event type: %s' % event_type)
        IN = h5py.File(os.path.join(options.outdir, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, options.confidence)), 'r')
        if 'conf_idx' in IN and IN['conf_idx'].shape[0] > 0:
            conf_idx = IN['conf_idx'][:]
            k_idx = sp.where(sp.isin(IN['gene_idx'][:][conf_idx], gid))[0]
            if k_idx.shape[0] > 0:
                event_info.extend([[event_type, x] for x in conf_idx[k_idx]])
        IN.close()

    return sp.array(event_info, dtype='str')


def get_seg_counts(gid, outdir, confidence, validate_sg):

    if validate_sg:
        IN = h5py.File(os.path.join(outdir, 'spladder', 'genes_graph_conf%i.merge_graphs.validated.count.hdf5' % confidence), 'r')
    else:
        IN = h5py.File(os.path.join(outdir, 'spladder', 'genes_graph_conf%i.merge_graphs.count.hdf5' % confidence), 'r')
    idx = sp.where(IN['gene_ids_edges'][:] == gid)[0]
    edges = IN['edges'][idx, :]
    edge_idx = IN['edge_idx'][:][idx].astype('int')
    idx = sp.where(IN['gene_ids_segs'][:] == gid)[0]
    segments = IN['segments'][idx, :]
    samples = hp.decodeUTF8(IN['samples'][:])
    IN.close()

    return (segments, edges, edge_idx, samples)

def stack_exons(exons1, exons2):
    
    return sp.r_[exons1, exons2]
