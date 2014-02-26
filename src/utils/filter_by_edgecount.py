def filter_by_edgecount(genes, CFG):

    ### filter splicegraphs by support count over samples
    for i in len(genes):
        k_idx = sp.where(genes[i].splicegraph.edges.sum(axis = 1) == 0)[0]
        genes[i].splicegraph.edges = (genes[i].edge_count >= CFG['sg_min_edge_count'])
        ### remove all exons that have no incoming or outgoing edges (caused by validation, keep single exon transcripts that occured before)
        k_idx2 = sp.where(genes[i].splicegraph.edges.sum(axis = 1) == 0)[0]
        rm_idx = sp.where(~sp.in1d(k_idx2, k_idx))[0]
        keep_idx = sp.where(~sp.in1d(sp.array(range(genes[i].splicegraph.edges.shape[0])), rm_idx))[0]
        if keep_idx.shape[0] > 0:
            genes[i].subset(keep_idx)
        else:
            genes = []

    return genes
