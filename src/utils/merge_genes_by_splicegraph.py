def merge_genes_by_splicegraph(CFG, chunk_idx=None):
    # merge_genes_by_splicegraph(CFG, chunk_idx) 
    #
    #   This script takes several gene structures and merges them into one. 
    #   the merge is based on the splicegraphs within the genes struct
    
    ### if we are running with rproc we only get one parameter struct
    if chunk_idx is None and isinstance(CFG, dict):
        if 'chunk_idx' in CFG:
            chunk_idx = CFG['chunk_idx']
        CFG = CFG['CFG']

    ### generate merge list
    merge_list = []
    if CFG['do_prune']:
        prune_tag = '_pruned'
    else:
        prune_tag = ''

    ### subset samples in case of chunked computation
    if chunk_idx is not None:
        samples = CFG['samples'][chunk_idx]
    else:
        samples = CFG['samples']

    ### add all single bam file graphs
    for i in range(len(samples)):
        merge_list.append('%s/spladder/genes_graph_conf%i.%s%s.mat', CFG['out_dirname'], CFG['confidence_level'], CFG['samples'][i], prune_tag)

    ### add also graph of all bam files combined
    if CFG['do_merge_all'] and os.path.exists('%s/spladder/genes_graph_conf%i.merge_bams%s.mat' % (CFG['out_dirname'], CFG['confidence_level'], prune_tag)):
        merge_list.append('%s/spladder/genes_graph_conf%i.merge_bams%s.mat' % (CFG['out_dirname'], CFG['confidence_level'], prune_tag))

    ### iterate over merge list 
    appended = False
    for i in range(len(merge_list)):
        ### load gene structure from sample i
        print 'Loading %s ...' % merge_list[i]
        genes = cPickle.load(open(merge_list[i], 'r'))
        print '... done (%i / %i)' % (i, len(merge_list))
        assert(genes.splicegraph is not None)

        ### sort genes by name
        name_list = sp.array([x.name for x in genes])
        s_idx = sp.argsort(name_list)
        genes = genes[s_idx]

        ### make sure, that splicegraph is unique
        genes.splicegraph.uniquify()

        ### jump over first sample - nothig to add yet 
        if i == 0:
            genes2 = genes.copy()
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
        print 'Processing ...'
        for j in range(genes.shape[0]):
            if j % 100 == 0:
                print '.',
                if j % 1000 == 0:
                    print '%i/%i' % (j, genes.shape[0])
            g_idx_ = g_idx
            while (g_idx <= genes2.shape[0] and genes2[g_idx].name < genes[j].name):
                g_idx += 1

            # same gene
            if g_idx <= genes2.shape[0] and genes2[g_idx].name == genes[j].name:
                
                tmp, s_idx = sort_rows(genes[j].splicegraph.vertices.T, index=True) 
                genes[j].splicegraph.vertices = tmp.T

                splice1 = genes[j].splicegraph.edges[s_idx, :][:, s_idx]
                splice2 = genes2[g_idx].splicegraph.edges

                s1_len = genes[j].splicegraph.vertices.shape[1]
                s2_len = genes2[g_idx].splicegraph.vertices.shape[1]
          
                if s2_len > 10000:
                    print 'Do not further merge into gene %i, has more than 10000 vertices!' % g_idx
                    ### still count edges that can be confirmed
                    tmp, c_idx, a_idx = intersect_rows(genes2[g_idx].splicegraph.vertices.T, genes[j].splicegraph.vertices.T, index=True)
                    if c_idx.shape[0] > 0:
                        genes2[g_idx].edge_count[c_idx, :][:, c_idx] = genes2[g_idx].edge_count[c_idx, :][:, c_idx] + genes[j].splicegraph.edges[a_idx, :][:, a_idx]
                else:
                    m_graph = sp.r_[sp.c_[genes[j].splicegraph.vertices.T, sp.ones((s1_len, 1))], sp.c_[genes2[g_idx].splicegraph.vertices.T, 2 * sp.ones((s2_len, 1))]]
                    tmp, s_idx = sort_rows(m_graph[:, 0:3], index=True)
                    m_graph = m_graph[s_idx, :]

                    um_graph, u_f = unique(m_graph[:, 0:2], index=True)
                    u_l = sp.r_[u_f[1:] - 1, m_graph.shape[0] - 1]
                    u_graph = u_l - u_f
                    
                    f_idx = sp.where(u_graph == 0)[0]

                    if f_idx.shape[0] > 0:
                        splice1_ = sp.zeros((u_graph.shape[0], u_graph.shape[0]))
                        splice2_ = sp.zeros((u_graph.shape[0], u_graph.shape[0]))
                        edgecnt = sp.zeros((u_graph.shape[0], u_graph.shape[0]))
                        idx1_ = sp.where(m_graph[u_f, 2] == 1)[0]
                        idx2_ = sp.where(m_graph[u_l, 2] == 2)[0]
                        splice1_[idx1_, :][:, idx1_] = splice1
                        splice2_[idx2_, :][:, idx2_] = splice2
                        edgecnt[idx2_, :][:, idx2_] = genes2[g_idx].edge_count
                    else:
                        splice1_ = splice1
                        splice2_ = splice2
                        edgecnt = genes2[g_idx].edge_count

                    if not sp.all(splice1_.shape == splice2_.shape):
                        print >> sys.stderr, 'ERROR: splice1_ and splice2_ differ in size!'
                        sys.exit(1)

                    genes2[g_idx].splicegraph.edges = (splice1_ | splice2_)
                    genes2[g_idx].splicegraph.vertices = um_graph.T
                    genes2[g_idx].splicegraph.terminals = sp.r_[(sp.tril(genes2[g_idx].splicegraph.edges).sum(axis=1) == 0).T.astype('int'), 
                                                                (sp.triu(genes2[g_idx].splicegraph.edges).sum(axis=1) == 0).T.astype('int')]
                    genes2[g_idx].edge_count = edgecnt + splice1_
            ### we did not find the gene name --> append new gene to genes2
            elif g_idx > genes2.shape[0] or genes2[g_idx].name > genes[j].name:
                g_idx = g_idx_
                genes2 = sp.c_[genes2, genes[j]]
                genes2[-1].edge_count = genes[j].splicegraph.edges
                appended = True
        print '... done\n'
        del genes

    genes = genes2.copy()
    del genes2

    fn_out = '%s/spladder/genes_graph_conf%i.%s%s.mat' % (CFG['out_dirname'], CFG['confidence_level'], CFG['merge_strategy'], prune_tag)
    if chunk_idx is not None:
        chunk_tag = '_chunk%i_%i' % (chunk_idx[0], chunk_idx[-1])
        fn_out = fn_out.replace('.mat', '%s.mat' % chunk_tag)
    
    print 'Store genes at: %s' % fn_out
    cPickle.dump(genes, open(fn_out, 'w'), -1)
