import random
import os
import sys
import pickle
import math
import glob
import numpy as np

if __package__ is None:
    __package__ = 'modules'

from .utils import *
from .classes.segmentgraph import Segmentgraph


def merge_duplicate_exons(genes, options):

    num_removed = 0
    for i in range(genes.shape[0]):
        if options.verbose and i % 1000 == 0: 
            print('%i\r' % i)

        ### check if there are non unique exons
        if unique_rows(genes[i].splicegraph.vertices.T).shape[0] == genes[i].splicegraph.vertices.shape[1]:
            continue

        genes[i].splicegraph.sort()

        exons = genes[i].splicegraph.vertices.copy()
        admat = genes[i].splicegraph.edges.copy()
        initial = genes[i].splicegraph.terminals[0, :]
        terminal = genes[i].splicegraph.terminals[1, :]

        remove = np.zeros((exons.shape[1],), dtype='bool')
        for j in range(exons.shape[1]):
            if remove[j]:
                continue
            idx = np.where((exons[0, j] == exons[0, :]) & (exons[1, j] == exons[1, :]))[0]
            idx = idx[np.where(~np.in1d(idx, j))[0]]
            if idx.shape[0] > 0:
                remove[idx] = True
                for k in idx:
                    initial[j] = (initial[j] | initial[k])
                    terminal[j] = (terminal[j] | terminal[k])
                    admat[j, :] = (admat[j, :] | admat[k, :])
                    admat[:, j] = (admat[:, j] | admat[:, k])

        if np.sum(remove) == 0:
            continue

        genes[i].splicegraph.edges = admat.copy()
        genes[i].splicegraph.terminals[0, :] = initial
        genes[i].splicegraph.terminals[1, :] = terminal

        k_idx = np.where(~remove)[0]
        genes[i].splicegraph.vertices = genes[i].splicegraph.vertices[:, k_idx]
        genes[i].splicegraph.edges = genes[i].splicegraph.edges[k_idx, :][:, k_idx]
        genes[i].splicegraph.terminals = genes[i].splicegraph.terminals[:, k_idx]
        num_removed += sum(remove)

    if options.verbose:
        print('\n... removed %i duplicate exons ...' % num_removed)

    return genes


def merge_genes_by_splicegraph(options, merge_list=None, fn_out=None):
    #
    #   This script takes several gene structures and merges them into one. 
    #   the merge is based on the splicegraphs within the genes struct
    
    ### if we are running with rproc we only get one parameter struct
    if merge_list is None and 'options' in options:
        if 'merge_list' in options:
            merge_list = options['merge_list']
        if 'fn_out' in options:
            fn_out = options['fn_out']
        options = options['options']

    ### add also graph of all bam files combined
    if options.do_merge_all and os.path.exists('%s/spladder/genes_graph_conf%i.merge_bams.pickle' % (options.outdir, options.confidence)):
        merge_list.append('%s/spladder/genes_graph_conf%i.merge_bams.pickle' % (options.outdir, options.confidence))

    ### iterate over merge list 
    appended = False
    for i in range(len(merge_list)):
        ### load gene structure from sample i
        print('Loading %s ...' % merge_list[i])
        (genes, inserted) = pickle.load(open(merge_list[i], 'rb'))
        for g in genes:
            g.from_sparse()
        print('... done (%i / %i)' % (i + 1, len(merge_list)))

        ### sort genes by name
        name_list = np.array([x.name for x in genes])
        s_idx = np.argsort(name_list)
        genes = genes[s_idx]

        ### make sure, that splicegraph is unique
        for g in genes:
            g.splicegraph.uniquify()

        ### skip first sample - nothig to add yet 
        if i == 0:
            genes2 = genes.copy()
            if not hasattr(genes2[0], 'edge_count'):
                for j in range(genes.shape[0]):
                    genes2[j].edge_count = genes2[j].splicegraph.edges.copy()
            del genes
            continue

        ### did we append genes in the last round? --> re-sort
        if appended:
            name_list = np.array([x.name for x in genes2])
            s_idx = np.argsort(name_list)
            genes2 = genes2[s_idx]
            appended = False

        ### iterate over current genes
        g_idx = 0
        print('Processing ...')
        for j in range(genes.shape[0]):
            if j % 100 == 0:
                print('.', end=' ')
                if j % 1000 == 0:
                    print('%i/%i' % (j, genes.shape[0]))
            g_idx_ = g_idx
            while (g_idx <= genes2.shape[0] and genes2[g_idx].name < genes[j].name):
                g_idx += 1

            # same gene
            if g_idx <= genes2.shape[0] and genes2[g_idx].name == genes[j].name:
                
                tmp, s_idx = sort_rows(genes[j].splicegraph.vertices.T, index=True) 
                genes[j].splicegraph.vertices = tmp.T

                splice1 = genes[j].splicegraph.edges[s_idx, :][:, s_idx].copy()
                splice2 = genes2[g_idx].splicegraph.edges.copy()
                try:
                    edgecnt1 = genes[j].edge_count[s_idx, :][:, s_idx].copy()
                except AttributeError:
                    edgecnt1 = splice1
                edgecnt2 = genes2[g_idx].edge_count

                s1_len = genes[j].splicegraph.vertices.shape[1]
                s2_len = genes2[g_idx].splicegraph.vertices.shape[1]

                if s2_len > 10000:
                    print('Do not further merge into gene %i, has more than 10000 vertices!' % g_idx)
                    ### still count edges that can be confirmed
                    tmp, c_idx, a_idx = intersect_rows(genes2[g_idx].splicegraph.vertices.T, genes[j].splicegraph.vertices.T, index=True)
                    if c_idx.shape[0] > 0:
                        genes2[g_idx].edge_count = replace_sub_matrix(genes2[g_idx].edge_count, c_idx, genes2[g_idx].edge_count[c_idx, :][:, c_idx] + edgecnt1[a_idx, :][:, a_idx])
                else:
                    m_graph = np.r_[np.c_[genes[j].splicegraph.vertices.T, np.ones((s1_len, 1), dtype='int')], np.c_[genes2[g_idx].splicegraph.vertices.T, 2 * np.ones((s2_len, 1), dtype='int')]]
                    tmp, s_idx = sort_rows(m_graph[:, 0:3], index=True)
                    m_graph = m_graph[s_idx, :]

                    um_graph, u_f = unique_rows(m_graph[:, 0:2], index=True)
                    u_l = np.r_[u_f[1:] - 1, m_graph.shape[0] - 1]
                    u_graph = u_l - u_f
                    
                    f_idx = np.where(u_graph == 0)[0]

                    if f_idx.shape[0] > 0:
                        splice1_ = np.zeros((u_graph.shape[0], u_graph.shape[0]), dtype='int')
                        splice2_ = np.zeros((u_graph.shape[0], u_graph.shape[0]), dtype='int')
                        edgecnt1_ = np.zeros((u_graph.shape[0], u_graph.shape[0]), dtype='int')
                        edgecnt2_ = np.zeros((u_graph.shape[0], u_graph.shape[0]), dtype='int')
                        idx1_ = np.where(m_graph[u_f, 2] == 1)[0]
                        idx2_ = np.where(m_graph[u_l, 2] == 2)[0]

                        splice1_ = replace_sub_matrix(splice1_, idx1_, splice1)
                        splice2_ = replace_sub_matrix(splice2_, idx2_, splice2)
                        edgecnt1_ = replace_sub_matrix(edgecnt1_, idx1_, edgecnt1)
                        edgecnt2_ = replace_sub_matrix(edgecnt2_, idx2_, edgecnt2)
                    else:
                        splice1_ = splice1
                        splice2_ = splice2
                        edgecnt1_ = edgecnt1
                        edgecnt2_ = edgecnt2

                    if not np.all(splice1_.shape == splice2_.shape):
                        print('ERROR: splice1_ and splice2_ differ in size!', file=sys.stderr)
                        sys.exit(1)

                    genes2[g_idx].splicegraph.edges = (splice1_ | splice2_)
                    genes2[g_idx].splicegraph.vertices = um_graph.T
                    genes2[g_idx].splicegraph.terminals = np.vstack([(np.tril(genes2[g_idx].splicegraph.edges).sum(axis=1) == 0),
                                                                     (np.triu(genes2[g_idx].splicegraph.edges).sum(axis=1) == 0)]).astype('int')
                    genes2[g_idx].edge_count = edgecnt1_ + edgecnt2_


            ### we did not find the gene name --> append new gene to genes2
            elif g_idx > genes2.shape[0] or genes2[g_idx].name > genes[j].name:
                g_idx = g_idx_
                genes2 = np.c_[genes2, genes[j]]
                genes2[-1].edge_count = genes[j].splicegraph.edges
                appended = True
        print('... done\n')
        del genes

    genes = genes2.copy()
    del genes2

    for g in genes:
        g.label_alt()
        if not g.segmentgraph.is_empty():
            g.segmentgraph = Segmentgraph(g)
        g.to_sparse()

    ### re-sort genes by position - makes quantification more efficient
    s1_idx = np.argsort([x.start for x in genes])
    s2_idx = np.argsort([x.chr for x in genes[s1_idx]], kind='mergesort')
    genes = genes[s1_idx[s2_idx]]

    print('Store genes at: %s' % fn_out)
    pickle.dump((genes, inserted), open(fn_out, 'wb'), -1)


def run_merge(samples, options):

    merge_all = (options.merge == 'merge_all')
    merge_all_tag = ''
    if merge_all:
        merge_all_tag = '_merged_bams'

    fn_out = '%s/spladder/genes_graph_conf%i.%s.pickle' % (options.outdir , options.confidence, options.merge)
    if options.validate_sg:
        fn_out_count = '%s/spladder/genes_graph_conf%i.%s.validated.count.hdf5' % (options.outdir, options.confidence, options.merge)
    else:
        fn_out_count = '%s/spladder/genes_graph_conf%i.%s.count.hdf5' % (options.outdir, options.confidence, options.merge)

    if not os.path.exists(fn_out):
        if len(options.chunked_merge) > 0:
            curr_level, max_level, chunk_start, chunk_end = [int(_) for _ in options.chunked_merge[0]]
            if curr_level == 1:
                merge_list = np.array(sorted(['%s/spladder/genes_graph_conf%i.%s.pickle' % (options.outdir, options.confidence, x) for x in samples]))
            else:
                merge_list = np.array(sorted(glob.glob('%s/spladder/genes_graph_conf%i.%s_level%i_chunk*.pickle' % (options.outdir, options.confidence, options.merge, curr_level - 1))))
            chunk_end = min(len(merge_list), chunk_end)

            if curr_level == max_level:
                assert len(merge_list) <= options.chunksize, 'chunksize is %i but merge_list has length %i with: %s' % (options.chunksize, len(merge_list), str(merge_list))
                fn = fn_out
            else:
                fn = '%s/spladder/genes_graph_conf%i.%s_level%i_chunk%i_%i.pickle' % (options.outdir, options.confidence, options.merge, curr_level, chunk_start, chunk_end)
            ### nothing to be done
            if os.path.exists(fn):
                print('%s already exists' % fn)
                return
            print('merging level %i chunk %i to %i' % (curr_level, chunk_start, chunk_end))
            merge_genes_by_splicegraph(options, merge_list=merge_list[chunk_start:chunk_end], fn_out=fn)
        else:
            merge_list = np.array(['%s/spladder/genes_graph_conf%i.%s.pickle' % (options.outdir, options.confidence, x) for x in samples])
            merge_genes_by_splicegraph(options, merge_list=merge_list, fn_out=fn_out)
    else:
        print('File %s already exists!' % fn_out)

