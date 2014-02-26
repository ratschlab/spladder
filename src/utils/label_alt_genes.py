def label_alt_genes(genes, CFG):
    # genes = label_alt_genes(genes, CFG) ;
    # 
    # This script labels each gene 'is_alt_spliced' that
    # has any transcript position that can be intronic or
    # exonic.
    #
    # This script labels each gene 'is_alt', that has a
    # simple alternative transcript start or end. Simple meaning
    # that the remaining transcript structure is the same and no 
    # other exons overlap to the alternative start or end.

    tot_exons = 0
    for ix in range(genes.shape[0]):
        if CFG['verbose'] and ix % 100 == 0:
            print >> CFG['fd_log'], '.',
            if ix % 1000 == 0:
                print >> CFG['fd_log'], '%i' % ix
      
        num_exons = genes[idx].splicegraph.vertices.shape[1]
        tot_exons += num_exons
      
        vertices = genes[idx].splicegraph.vertices
        edges = genes[ix].splicegraph.edges
       
        ### no edges in graph --> continue
        if edges.shape[0] == 0:
            genes[idx].is_alt_spliced = 0
            genes[idx].is_alt = 0
            continue 
        
        init_alt_idx = []
        term_alt_idx = []
        for i in range(num_exons):
            ### handle start terminal exons -> no incoming edges
            if not sp.any(edges[:i, i]):
                ### find other exons with same end
                idx = sp.where(~sp.in1d(sp.where(vertices[1, i] == vertices[1, :])[0], i))[0]
                if idx.shape[0] > 0:
                    is_simple = True
                    for j in idx:
                        ### not simple, if different introns follow --> break
                        ii1 = sp.where(edges[j, j+1:])[0] + j
                        ii2 = sp.where(edges[i, i+1:])[0] + i
                        if ii1.shape != ii2.shape or not sp.all(ii1 == ii2):
                            is_simple = False
                            break
                        ### not simple, if exons previous to j overlap the start of i --> break
                        idx_prev = sp.where(edges[:j, j])[0]
                        for k in idx_prev:
                            if vertices[1, k] > vertices[0, i]:
                                is_simple = False
                                break
                    if not is_simple:
                        continue
                    ### if simple, report alternative initial exon
                    init_alt_idx.append(i)
            ### handle end terminal exons -> no outgoing edges
            if not sp.any(edges[i, i+1:]):
                ### find other exons with the same start
                idx = sp.where(~sp.in1d(sp.where(vertices[0, i] == vertices[0, :])[0], i))[0]
                if idx.shape[0] > 0: 
                    is_simple = True
                    for j in idx:
                        ### not simple, if different introns precede --> break
                        ii1 = sp.where(edges[:j, j])[0]
                        ii2 = sp.where(edges[:i, i])[0]
                        if ii1.shape != ii2.shape or not sp.all(ii1 == ii2):
                            is_simple = False
                            break
                        ### not simple, if exons following to j start before i ends --> break
                        idx_next = sp.where(edges[j, j+1:])[0] + j
                        for k in idx_next:
                            if vertices[0, k] < vertices[1, i]:
                                is_simple = False
                                break
                    if not is_simple:
                        continue
                    ### if simple, report alternative terminal exon
                    term_alt_idx.append(i)
      
        ### further only consider exons that are neither init_alt nor term_alt
        take_idx = sp.where(~sp.in1d(range(num_exons), [init_alt_idx, term_alt_idx]))[0]
        vertices = genes[idx].splicegraph.vertices[:, take_idx]
        edges = genes[idx].splicegraph.edges[take_idx, :][:, take_idx]
      
        start = vertices[0, :].min()
        stop = vertices[1, :].max()
        
        exon_loc = sp.zeros((stop - start,))
        
        ### map all edges (introns) to genomic coordinates
        for i in range(edges.shape[0]):
            for j in range(i + 1 : edges.shape[0]):
                if edges[i, j] == 1:
                    cur_edge = sp.array([vertices[1, i], vertices[0, j]]) - start
                    exon_loc[cur_edge[0]:cur_edge[1]] += 1
      
        ### map all vertices (exons) to genomic coordinates 
        for i in range(vertices.shape[1]):
          cur_vertex = vertices[:, i] - start
          exon_loc(cur_vertex[0]:cur_vertex[1]) += 1
      
        ### if at any position more than one exon or intron -> is_alt__spliced
        if max(exon_loc) > 1:
            genes[idx].is_alt_spliced = 1
            genes[idx].is_alt = 1
        else:
            genes[idx].is_alt_spliced = 0
            ### if not alt_spliced but term_alt or init_alt --> is_alt
            if init_alt_idx.shape[0] == 0 and term_alt_idx.shape[0] == 0:
                genes[idx].is_alt = 0 
            else:
                genes[idx].is_alt = 1 
      
    if CFG['verbose']:
        print >> CFG['fd_log'],'\n\nTotal genes:\t\t\t\t\t\t\t%d' % genes.shape[0]
        print >> CFG['fd_log'],'Total exons:\t\t\t\t\t\t\t%d' % tot_exons
        print >> CFG['fd_log'],'Total genes with alternative isoforms:\t\t\t\t%d' % sp.sum([x.is_alt for x in genes])
        print >> CFG['fd_log'],'Total genes alternatively spliced:\t\t\t\t%d' % sp.sum([x.is_alt_spliced for x in genes])
        print >> CFG['fd_log'],'Total constitutively spliced:\t\t\t\t\t%d' % genes.shape[0] - sp.sum([x.is_alt_spliced for x in genes])
