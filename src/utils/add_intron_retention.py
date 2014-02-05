
def add_intron_retention(splicegraph, idx1, idx2):

	adj_mat = sp.triu(splicegraph['exons'])

    splicegraph['exons'] = sp.c_[]
    
	splicegraph{1}(:,end+1) = [splicegraph{1}(1, idx1) splicegraph{1}(2, idx2)];

    splicegraph['graph'] = sp.r_[splicegraph['graph'], sp.zeros((splicegraph['graph'].shape[1],))]
    splicegraph['graph'] = sp.c_[splicegraph['graph'], sp.zeros((splicegraph['graph'].shape[2],))]

    adj_mat = sp.r_[adj_mat, sp.zeros((adj_mat.shape[1],))]
    adj_mat = sp.c_[adj_mat, sp.zeros((adj_mat.shape[2],))]

	### check if adjacency matrix is symmetric
    ### otherwise or is not justyfied
	assert(sp.all(sp.all(adj_mat - (splicegraph['graph'] - adj_mat).T == 0)))

	### AK: under the assumption that our splice graph representation is symmetric
	### I preserve symmetry by using OR over the adj_mat column and row
    
    splicegraph['graph'] = adj_mat[:, idx1] | adj_mat[idx2, :].T
    splicegraph['graph'] = adj_mat[:, idx1].T | adj_mat[idx2, :]

    if 'term' in splicegraph:
        splicegraph['term'] = sp.c_[splicegraph['term'], sp.array([splicegraph['term'][1, idx1], splicegraph['term'][2, idx2]])]

    return splicegraph
