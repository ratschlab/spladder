from scipy.sparse import lil_matrix
from scipy.matlib import repmat

def detect_multipleskips(genes, idx_alt):
    # [idx_multiple_skips, exon_multiple_skips, id_multiple_skips] = detect_multipleskips(genes, idx_alt) ;

    id = 0
    idx_multiple_skips = []
    id_multiple_skips = []
    exon_multiple_skips = []
    for ix in idx_alt:
        if ix % 50 == 0:
            print '.',
        num_exons = genes[ix].splicegraph.get_len()
        edges = genes[ix].splicegraph.edges
        labels = repmat(sp.arange(num_exons), 1, num_exons)
        
        # adjecency matrix: upper half only
        A = sp.zeros((num_exons, num_exons))
        for i in range(num_exons - 1):
            for j in range(i + 1, num_exons):
                A[i, j] = edges[i, j]
        
        # possible starting and ending exons of a multiple exon skip
        Pairs = lil_matrix((num_exons, num_exons))
        Ai = sp.dot(sp.dot(A, A) * A) #paths of length 3
        while sp.any(Ai.ravel() > 0):
            coords = sp.where((A > 0 ) & (Ai > 0)) # multiple skip
            Pairs[coords[0], coords[1]] = 1
            Ai = sp.dot(Ai, A)  # paths of length ..+1
        
        edge = sp.where(Pairs == 1)
        
        if edge[0].shape[0] > 10000,
            print 'Warning: not processing gene %d, because there are more than 10000 potential hits.' % ix
            continue
        
        for cnt in range(edge[0].shape[0]):
            exon_idx_first = edge[0][cnt]
            exon_idx_last = edge[1][cnt]
      
            if edges[exon_idx_first, exon_idx_last] == 1:
                
                # find all pairs shortest path
                exist_path = triu(edges).astype('double')
                exist_path[exon_idx_first, exon_idx_last] = 0
                exist_path[exist_path == 0] = sp.inf
                # set diagonal to 0
                exist_path[sp.arange(exist_path.shape[0]), sp.arange(exist_path.shape[0])] = 0
                
                long_exist_path = -triu(edges).astype('double')
                long_exist_path[exon_idx_first, exon_idx_last] = 0
                long_exist_path[long_exist_path == 0] = sp.inf
                # set diagonal to 0
                long_exist_path[sp.arange(long_exist_path.shape[0]), sp.arange(long_exist_path.shape[0])] = 0
                
                path = sp.isfinite(exist_path) * labels
                long_path = sp.isfinite(long_exist_path) * labels
                
                for k in range(num_exons):
                    for i in range(num_exons):
                        idx = sp.where((exist_path[i, k] + exist_path[k, :]) < exist_path[i, :])[0]
                        exist_path[i, idx] = exist_path[i, k] + exist_path[k, idx]
                        path[i, idx] = path[k, idx]
                        
                        idx = sp.where((long_exist_path[i,k] + long_exist_path[k, :]) < long_exist_path[i, :])[0]
                        long_exist_path[i, idx] = long_exist_path[i, k] + long_exist_path[k, idx]
                        long_path[i, idx] = long_path[k, idx]
                
                temp_ix = sp.isfinite(long_exist_path)
                long_exist_path[temp_ix] = -long_exist_path[temp_ix]
                
                if (exist_path[exon_idx_first, exon_idx_last] > 2) and sp.isfinite(exist_path[exon_idx_first, exon_idx_last]):
                    backtrace = path[exon_idx_first, exon_idx_last]
                    while backtrace[-1] > exon_idx_first:
                        backtrace = sp.c_[backtrace, path[exon_idx_first, backtrace[-1]]]
                    backtrace = backtrace[:-1]
                    backtrace = backtrace[::-1]
                    idx_multiple_skips.append(repmat(ix, 1, backtrace.shape[1] + 2))
                    exon_multiple_skips.append([exon_idx_first, backtrace, exon_idx_last])
                    id += 1
                    id_multiple_skips.append(id * sp.ones((1, backtrace.shape[1] + 2))
                elif (long_exist_path(exon_idx_first,exon_idx_last) > 2) && ~isinf(long_exist_path(exon_idx_first,exon_idx_last)),
                    backtrace = long_path(exon_idx_first,exon_idx_last);
                    while backtrace[-1] > exon_idx_first:
                        backtrace = sp.c_[backtrace, long_path[exon_idx_first, backtrace[-1]]]
                    backtrace = backtrace[:-1]
                    backtrace = backtrace[::-1]
                    idx_multiple_skips.append(repmat(ix, 1, backtrace.shape[1] + 2))
                    exon_multiple_skips.append(exon_idx_first, backtrace, exon_idx_last])
                    id += 1
                    id_multiple_skips.append(id * sp.ones(1, backtrace.shape[1] + 2))

    print 'Number of multiple exon skips:\t\t\t\t\t%d' % length(idx_multiple_skips)
    return (idx_multiple_skips, exon_multiple_skips, id_multiple_skips)
