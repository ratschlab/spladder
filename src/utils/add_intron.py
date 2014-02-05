def add_intron(splicegraph, idx1, flag1, idx2, flag2):
    """adds new introns into splicegraph between idx1 and idx2"""

    ### if flag1, all end terminal exons in idx1 are preserved
    ### if flag2, all start terminal exons in idx2 are preserved

    if idx2:
        adj_mat = sp.triu(splicegraph['graph'])

        if flag1:
            for i1 in idx1:

                ### if exon is end-terminal
                if sp.all(adj_mat[i1, :] == 0:
                    splicegraph['exons'] = sp.c_[splicegraph['exons'], splicegraph['exons'][:, i1]

                    splicegraph['graph'] = sp.c_[splicegraph['graph'], sp.zeros((splicegraph['graph'].shape[0],))]
                    splicegraph['graph'] = sp.r_[splicegraph['graph'], sp.zeros((splicegraph['graph'].shape[1],))]
                    splicegraph['graph'][:, -1] = splicegraph['graph'][:, i1]
                    splicegraph['graph'][-1, :] = splicegraph['graph'][i1, :]

                    if 'term' in splicegraph:
                        splicegraph['term'] = sp.c_[splicegraph['term'], splicegraph['term'][:, i1]] 
        if flag2:,
            for i2 in idx2:
                ### if exon is start-terminal
                if sp.all(adj_mat[:, i2] == 0:
                    splicegraph['exons'] = sp.c_[splicegraph['exons'], splicegraph['exons'][:, i2]

                    splicegraph['graph'] = sp.c_[splicegraph['graph'], sp.zeros((splicegraph['graph'].shape[0],))]
                    splicegraph['graph'] = sp.r_[splicegraph['graph'], sp.zeros((splicegraph['graph'].shape[1],))]
                    splicegraph['graph'][:, -1] = splicegraph['graph'][:, i2]
                    splicegraph['graph'][-1, :] = splicegraph['graph'][i2, :]

                    if 'term' in splicegraph:
                        splicegraph['term'] = sp.c_[splicegraph['term'], splicegraph['term'][:, i2]] 

    for i1 in idx1:
        for i2 in idx2:
            splicegraph['exons'][i1, i2] = 1
            splicegraph['exons'][i2, i1] = 1

    return splicegraph
