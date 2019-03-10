import random
import os
import sys
import pickle
import math

if __package__ is None:
    __package__ = 'modules'

from .utils import *
from .count import count_graph_coverage_wrapper
from .editgraph import filter_by_edgecount
from . import rproc as rp


def merge_duplicate_exons(genes, options):

    num_removed = 0
    for i in range(genes.shape[0]):
        if options.verbose and i % 100 == 0: 
            print('%i\r' % i)

        ### check if there are non unique exons
        if unique_rows(genes[i].splicegraph.vertices.T).shape[0] == genes[i].splicegraph.vertices.shape[1]:
            continue

        genes[i].splicegraph.sort()

        exons = genes[i].splicegraph.vertices.copy()
        admat = genes[i].splicegraph.edges.copy()
        initial = genes[i].splicegraph.terminals[0, :]
        terminal = genes[i].splicegraph.terminals[1, :]

        remove = sp.zeros((exons.shape[1],), dtype='bool')
        for j in range(exons.shape[1]):
            if remove[j]:
                continue
            idx = sp.where((exons[0, j] == exons[0, :]) & (exons[1, j] == exons[1, :]))[0]
            idx = idx[sp.where(~sp.in1d(idx, j))[0]]
            if idx.shape[0] > 0:
                remove[idx] = True
                for k in idx:
                    initial[j] = (initial[j] | initial[k])
                    terminal[j] = (terminal[j] | terminal[k])
                    admat[j, :] = (admat[j, :] | admat[k, :])
                    admat[:, j] = (admat[:, j] | admat[:, k])

        if sp.sum(remove) == 0:
            continue

        genes[i].splicegraph.edges = admat.copy()
        genes[i].splicegraph.terminals[0, :] = initial
        genes[i].splicegraph.terminals[1, :] = terminal

        k_idx = sp.where(~remove)[0]
        genes[i].splicegraph.vertices = genes[i].splicegraph.vertices[:, k_idx]
        genes[i].splicegraph.edges = genes[i].splicegraph.edges[k_idx, :][:, k_idx]
        genes[i].splicegraph.terminals = genes[i].splicegraph.terminals[:, k_idx]
        num_removed += sum(remove)

    if options.verbose:
        print('\n... removed %i duplicate exons ...' % num_removed)

    return genes


def merge_genes_by_isoform(options):
    ### This script takes several gene structures and merges them into one. If the number of isoforms of an annotated transcript should 
    ### exceed max_num_isoforms, max_num_isoforms many are subsampled
    
    appended = False
    max_num_isoforms = 10

    ### generate merge list
    merge_list = []
    if options.do_prune:
        prune_tag = '_pruned'
    else:
        prune_tag = ''

    ### add all single bam file isoforms
    for i in range(len(options.samples)):
        merge_list.append('%s/spladder/genes_graph_conf%i.%s%s.pickle' % (options.outdir, options.confidence, options.samples[i], prune_tag))

    ### add also isoforms of all bam files combined
    fn_bam_merged = '%s/spladder/genes_graph_conf%i.merge_bams%s.pickle' % (options.outdir, options.confidence, prune_tag)
    if options.merge == 'merge_all' and os.path.exists(fn_bam_merged):
        merge_list.append(fn_bam_merged)

    for i in range(len(merge_list)):
        ### load gene structure from sample i
        print('Loading %s ...' % merge_list[i])
        (genes, inserted) = pickle.load(open(merge_list[i]), 'rb')
        for g in genes:
            g.from_sparse()
        print('... done')

        ### sort
        name_list = sp.array([x.name for x in genes])
        s_idx = sp.argsort(name_list)
        name_list = name_list[s_idx]
        genes = genes[s_idx]

        ### jump over first sample - nothig to add yet 
        if i == 0:
            genes2 = genes.copy()
            del genes
            continue

        ### did we append genes in the last round? --> re-sort
        if appended:
            name_list = sp.array([x.name for x in genes2])
            s_idx = sp.argsort(name_list)
            name_list = name_list[s_idx]
            genes2 = genes2[s_idx]
            appended = False

        ### iterate over current genes
        g_idx = 0
        print('Processing ...')
        for j in range(genes.shape[0]):
            if j % 100 == 0:
                print('.', end=' ')
                if j % 1000 == 0:
                    print('%i/%i' % (j + 1, genes.shape[0]))
            g_idx_ = g_idx

            while genes2[g_idx].name < genes[j].name:
                g_idx += 1
            # same gene
            if genes2[g_idx].name == genes[j].name:
                ### pairwise comparison of all exons
                for k in range(len(genes[j].exons)):
                    broken = False
                    for l in range(len(genes2[g_idx].exons)):
                        ### transcripts are not identical in size --> continue
                        if not sp.all(genes[j].exons[k].shape == genes2[g_idx].exons[l].shape):
                            continue 
                        ### found transcript in genes2 that is identical to current transcript --> break
                        if sp.all(genes[j].exons[k] == genes2[g_idx].exons[l]):
                            broken = True
                            break 
                    ### we did not find any identical transcript in genes2 --> append transcript
                    if not broken:
                        genes2[g_idx].exons.append(genes[j].exons[k])
                        genes2[g_idx].transcripts.append('%s_%i' % (genes[j].transcripts[k], i))
                        genes2[g_idx].start = min(genes2[g_idx].start, genes[j].exons[k][:, 0].min())
                        genes2[g_idx].stop = max(genes2[g_idx].stop, genes[j].exons[k][:, 1].max())
                        genes2[g_idx].splicegraph = Splicegraph(genes2[g_idx])
            ### we did non find the gene name --> append new gene to genes2
            elif genes2[g_idx].name > genes[j].name:
                g_idx = g_idx_
                genes2 = sp.c_[genes2, genes[j]]
                appended = True
        print('... done\n')
        del genes

    genes = genes2
    del genes2

    fn = '%s/spladder/genes_graph_conf%i.%s%s_merge_isoforms.pickle' % (options.outdir, options.confidence, options.merge, prune_tag)
    print('Store genes at: %s' % fn)
    for g in genes:
        g.to_sparse()
    pickle.dump((genes, inserted), open(fn, 'wb'), -1)

    ### subsample transcripts if neccessary 
    print('Subsample genes ...')
    for i in range(genes.shape[0]):
        if len(genes[i].transcripts) > max_num_isoforms:
            print('Subsample for gene %i from %i transcripts.' % (i, len(genes[i].transcripts)))
            for r_idx in sorted(random.sample(sp.arange(len(genes[i].transcripts)), len(genes[i].transcripts) - max_num_isoforms), reverse=True):
                del genes[i].transcripts[r_idx]
                del genes[i].exons[r_idx]
            genes[i].start = min(genes[i].start, genes[i].exons[0][:, 0].min())
            genes[i].stop = max(genes[i].stop, genes[i].exons[-1][:, 1].max())
            genes.splicegraph = Splicegraph(genes[i])
    print('... done\n')

    fn = '%s/spladder/genes_graph_conf%i.%s%s_merge_isoforms_subsampled.pickle' % (options.outdir, options.confidence, options.merge, prune_tag)
    print('Store subsampled genes at: %s' % fn)
    pickle.dump((genes, inserted), open(fn, 'wb'), -1)       


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

    ### generate merge list
    if options.do_prune:
        prune_tag = '_pruned'
    else:
        prune_tag = ''

    ### add also graph of all bam files combined
    if options.do_merge_all and os.path.exists('%s/spladder/genes_graph_conf%i.merge_bams%s.pickle' % (options.outdir, options.confidence, prune_tag)):
        merge_list.append('%s/spladder/genes_graph_conf%i.merge_bams%s.pickle' % (options.outdir, options.confidence, prune_tag))

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
        name_list = sp.array([x.name for x in genes])
        s_idx = sp.argsort(name_list)
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
            name_list = sp.array([x.name for x in genes2])
            s_idx = sp.argsort(name_list)
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
                    m_graph = sp.r_[sp.c_[genes[j].splicegraph.vertices.T, sp.ones((s1_len, 1), dtype='int')], sp.c_[genes2[g_idx].splicegraph.vertices.T, 2 * sp.ones((s2_len, 1), dtype='int')]]
                    tmp, s_idx = sort_rows(m_graph[:, 0:3], index=True)
                    m_graph = m_graph[s_idx, :]

                    um_graph, u_f = unique_rows(m_graph[:, 0:2], index=True)
                    u_l = sp.r_[u_f[1:] - 1, m_graph.shape[0] - 1]
                    u_graph = u_l - u_f
                    
                    f_idx = sp.where(u_graph == 0)[0]

                    if f_idx.shape[0] > 0:
                        splice1_ = sp.zeros((u_graph.shape[0], u_graph.shape[0]), dtype='int')
                        splice2_ = sp.zeros((u_graph.shape[0], u_graph.shape[0]), dtype='int')
                        edgecnt1_ = sp.zeros((u_graph.shape[0], u_graph.shape[0]), dtype='int')
                        edgecnt2_ = sp.zeros((u_graph.shape[0], u_graph.shape[0]), dtype='int')
                        idx1_ = sp.where(m_graph[u_f, 2] == 1)[0]
                        idx2_ = sp.where(m_graph[u_l, 2] == 2)[0]

                        splice1_ = replace_sub_matrix(splice1_, idx1_, splice1)
                        splice2_ = replace_sub_matrix(splice2_, idx2_, splice2)
                        edgecnt1_ = replace_sub_matrix(edgecnt1_, idx1_, edgecnt1)
                        edgecnt2_ = replace_sub_matrix(edgecnt2_, idx2_, edgecnt2)
                    else:
                        splice1_ = splice1
                        splice2_ = splice2
                        edgecnt1_ = edgecnt1
                        edgecnt2_ = edgecnt2

                    if not sp.all(splice1_.shape == splice2_.shape):
                        print('ERROR: splice1_ and splice2_ differ in size!', file=sys.stderr)
                        sys.exit(1)

                    genes2[g_idx].splicegraph.edges = (splice1_ | splice2_)
                    genes2[g_idx].splicegraph.vertices = um_graph.T
                    genes2[g_idx].splicegraph.terminals = sp.vstack([(sp.tril(genes2[g_idx].splicegraph.edges).sum(axis=1) == 0),
                                                                     (sp.triu(genes2[g_idx].splicegraph.edges).sum(axis=1) == 0)]).astype('int')
                    genes2[g_idx].edge_count = edgecnt1_ + edgecnt2_


            ### we did not find the gene name --> append new gene to genes2
            elif g_idx > genes2.shape[0] or genes2[g_idx].name > genes[j].name:
                g_idx = g_idx_
                genes2 = sp.c_[genes2, genes[j]]
                genes2[-1].edge_count = genes[j].splicegraph.edges
                appended = True
        print('... done\n')
        del genes

    genes = genes2.copy()
    del genes2

    for g in genes:
        g.label_alt()
        g.to_sparse()

    ### re-sort genes by position - makes quantification more efficient
    s1_idx = sp.argsort([x.start for x in genes])
    s2_idx = sp.argsort([x.chr for x in genes[s1_idx]], kind='mergesort')
    genes = genes[s1_idx[s2_idx]]

    print('Store genes at: %s' % fn_out)
    pickle.dump((genes, inserted), open(fn_out, 'wb'), -1)


def run_merge(options):

    merge_all = (options.merge == 'merge_all')
    merge_all_tag = ''
    if merge_all:
        merge_all_tag = '_merged_bams'

    prune_tag = ''
    if options.do_prune:
        prune_tag = '_pruned'

    chunksize = 10

    fn_out = '%s/spladder/genes_graph_conf%i.%s%s.pickle' % (options.outdir , options.confidence, options.merge, prune_tag)
    if options.validate_sg:
        fn_out_count = '%s/spladder/genes_graph_conf%i.%s%s.validated.count.hdf5' % (options.outdir, options.confidence, options.merge, prune_tag)
    else:
        fn_out_count = '%s/spladder/genes_graph_conf%i.%s%s.count.hdf5' % (options.outdir, options.confidence, options.merge , prune_tag)

    if not os.path.exists(fn_out):
        if not options.pyproc:
            merge_list = sp.array(['%s/spladder/genes_graph_conf%i.%s%s.pickle' % (options.outdir, options.confidence, x, prune_tag) for x in options.samples])
            merge_genes_by_splicegraph(options, merge_list=merge_list, fn_out=fn_out)
        else:
            jobinfo = []
            PAR = dict()
            PAR['options'] = options
            if chunksize > 0:
                levels = int(math.ceil(math.log(len(options.samples), chunksize)))
                level_files = dict()
                for level in range(1, levels + 1):
                    print('merging files on level %i' % level)
                    if level == 1:
                        merge_list = sp.array(['%s/spladder/genes_graph_conf%i.%s%s.pickle' % (options.outdir, options.confidence, x, prune_tag) for x in options.samples])
                    else:
                        merge_list = sp.array(level_files[level - 1])
                    level_files[level] = []
                    for c_idx in range(0, len(merge_list), chunksize):
                        if level == levels:
                            assert(len(merge_list) <= chunksize)
                            fn = fn_out
                        else:
                            fn = '%s/spladder/genes_graph_conf%i.%s%s_level%i_chunk%i_%i.pickle' % (options.outdir, options.confidence, options.merge, prune_tag, level, c_idx, min(len(merge_list), c_idx + chunksize))
                        level_files[level].append(fn)
                        if os.path.exists(fn):
                            continue
                        else:
                            print('submitting level %i chunk %i to %i' % (level, c_idx, min(len(merge_list), c_idx + chunksize)))
                            chunk_idx = sp.arange(c_idx, min(len(merge_list), c_idx + chunksize))
                            PAR['merge_list'] = merge_list[chunk_idx]
                            PAR['fn_out'] = fn
                            jobinfo.append(rp.rproc('merge_genes_by_splicegraph', PAR, 20000*level, options.options_rproc, 40*60))
                    rp.rproc_wait(jobinfo, 30, 1.0, -1)
            else:
                PAR['merge_list'] = options.samples
                PAR['fn_out'] = fn_out
                jobinfo.append(rp.rproc('merge_genes_by_splicegraph', PAR, 20000, options.options_rproc, 40*60))
                rp.rproc_wait(jobinfo, 30, 1.0, -1)
    else:
        print('File %s already exists!' % fn_out)

    ### generate validated version of splice graph
    #if options.validate_sg and not os.path.exists(fn_out_val):
    #    (genes, inserted) = cPickle.load(open(fn_out, 'r'))
    #    genes = filter_by_edgecount(genes, options)
    #    cPickle.dump((genes, inserted), open(fn_out_val, 'w'), -1)
    #    del genes

    ### count segment graph
    #if options.validate_sg:
    #   count_graph_coverage_wrapper(fn_out_val, fn_out_count, options)
    #else:
    #   count_graph_coverage_wrapper(fn_out, fn_out_count, options)

    if options.do_gen_isoforms:
        fn_out = '%s/spladder/genes_graph_conf%i.%s%s_isoforms.pickle' % (options.outdir, options.confidence, options.merge, prune_tag)
        if not os.path.exists(fn_out):
            if not options.pyproc:
                merge_genes_by_isoform(options.outdir, options.confidence, merge_all, experiment)
            else:
                jobinfo = [rp.rproc('merge_genes_by_isoform', PAR, 10000, options.options_rproc, 40*60)]
                rp.rproc_wait(jobinfo, 30, 1.0, 1)
        else:
            print('File %s already exists!' % fn_out)
