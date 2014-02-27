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


def detect_intronreten(genes, idx_alt):
    # [idx_intron_reten,intron_intron_reten] = detect_intronreten(genes,idx_alt) ;

    idx_intron_reten = []
    intron_intron_reten = []
    for ix in idx_alt:
        if ix % 50 == 0:
            print '.',
        num_exons = genes[ix].splicegraph.get_len()
        vertices = genes[ix].splicegraph.vertices
        edges = genes[ix].splicegraph.edges
        introns  = []
        for exon_idx in range(num_exons - 1):  # start of intron
            idx = sp.where(edges[exon_idx, exon_idx + 1 : num_exons] == 1)[0]
            if idx.shape[0]:
                continue
            idx += exon_idx
            for exon_idx2 in idx: # end of intron
                is_intron_reten = False
                for exon_idx1 in range(num_exons): # exon
                    # check that the exon covers the intron
                    if (vertices[1, exon_idx] >= vertices[0, exon_idx1]) and (vertices[0, exon_idx2] <= vertices[1, exon_idx1]):
                        is_intron_reten = True
                        long_exon = exon_idx1 
                        for len in range(len(introns):
                            if (vertices[1, exon_idx] == introns[len][0]) and (vertices[0, exon_idx2] == introns[len][1]):
                                is_intron_reten = False
                if is_intron_reten:
                    idx_intron_reten.append(ix)
                    intron_intron_reten.append([exon_idx, exon_idx2, long_exon])
                    introns.append([vertices[1, exon_idx], vertices[0, exon_idx2]])

    print '\nNumber of intron retentions:\t\t\t\t\t%d', len(idx_intron_reten)
    return (idx_intron_reten, intron_intron_reten)


def detect_exonskips(genes, idx_alt):
    # [idx_exon_skips, exon_exon_skips] = detect_exonskips(genes, idx_alt) ;

    idx_exon_skips = []
    exon_exon_skips = []
    for ix in idx_alt:
        if ix % 50 == 0:
            print '.',
        num_exons = genes[ix].splicegraph.get_len()
        edges = genes[ix].splicegraph.edges
        for exon_idx in range(num_exons - 2): #first exon
            for exon_idx1 in range(exon_idx + 1, num_exons - 1): # middle exon
                for exon_idx2 in range(exon_idx1 + 1, num_exons): # last exon
                    if (edges[exon_idx, exon_idx1] == 1) and edges[exon_idx, exon_idx2] and edges[exon_idx1, exon_idx2]:
                        idx_exon_skips.append(ix)
                        exon_exon_skips.append([exon_idx, exon_idx1, exon_idx2])
    print '\nNumber of single exon skips:\t\t\t\t\t%d', len(idx_exon_skips))
    return (idx_exon_skips, exon_exon_skips)



def detect_altprime(genes, idx_alt):
    # function [idx_alt_5prime,exon_alt_5prime, idx_alt_3prime,exon_alt_3prime] ...
    #    = detect_altprime(genes, idx_alt);
    #
    # detect the alternative 5 and 3 prime ends of the intron. Note that 5 prime refers to the left
    # and 3 prime to the right for a positive strand

    MIN_OVERLAP = 11 # two exons that are alternative should overlap by at least the length

    idx_alt_5prime = []
    idx_alt_3prime = []
    exon_alt_5prime = {}
    exon_alt_3prime = {}

    for ix in idx_alt:
        if idx % 50 == 0:
            print '.',
        num_exons = genes[ix].splicegraph.get_len()
        vertices = genes[ix].splicegraph.vertices
        edges = genes[ix].splicegraph.edges
        strand = genes[ix].strand
        # Find alternative sites on the right of the intron,
        # same site on the left of the intron.
        for exon_idx in range(num_exons - 2):
            rightsites = []
            rightidx = []
            nr_exons = sp.sum(edges[exon_idx, exon_idx + 1 : num_exons])
            if nr_exons >= 2:
                which_exons = sp.where(edges[exon_idx, exon_idx + 1 : num_exons])[0] + exon_idx
                exons = vertices[:, which_exons]
                for i in range(nr_exons - 1):
                    for j in range(i + 1, nr_exons):
                        # check that the left splice site of the exons are different
                        # make sure exons overlap - either:
                        # - left splice site of exon(i) is in exon(j)
                        # - left splice site of exon(j) is in exon(i)
                        # note that the 'overlap' relationship is not transitive
                        if ((exons[0, i] != exons[0, j]) and 
                          (((exons[0, i] >= exons[0, j]) and (exons[0, i] <= exons[1, j])) or ((exons[0, j] >= exons[0, i]) and (exons[0, j] <= exons[1, i]))) and
                          (min(exons[1, i], exons[1, j]) - max(exons[0, i], exons[0, j]) >= MIN_OVERLAP)):

                            assert(not ((exons[0, i] == exons[0, j]) and (exons[1, i] == exons[1, j])));
                            assert(exons[0, i] != exons[0, j])

                            ### add new events to the list
                            if not exons[0, i] in rightsites:
                                rightsites.append(exons[0, i])
                                rightidx.append(which_exons[i])
                            if not exons[0, j] in rightsites:
                                rightsites.append(exons[0, j])
                                rightidx.append(which_exons[j])
           
            # construct the output
            if len(rightsites) >= 2:
                if strand == '+':
                    exon_alt_3prime.append({'fiveprimesite':exon_idx, 'threeprimesites':rightidx})
                    idx_alt_3prime.append(ix)
                if strand == '-':
                    exon_alt_5prime.append({'threeprimesite':exon_idx, 'fiveprimesites':rightidx})
                    idx_alt_5prime.append(ix)
        
        # Find alternative sites on the left of the intron,
        # same site on the right of the intron.
        for exon_idx in range(2, num_exons):
            nr_exons = sp.sum(edges[:exon_idx-1, exon_idx])
            leftsites = []
            leftidx = []
            if nr_exons >= 2:
                which_exons = sp.where(edges[:exon_idx-1, exon_idx])[0]
                exons = vertices[:, which_exons]
                for i in range(nr_exons - 1):
                    for j in range(i + 1, nr_exons):
                        # check that the 5prime sites are different
                        # make sure exons overlap - either:
                        # - right splice site of exon(i) is in exon(j)
                        # - right splice site of exon(j) is in exon(i)
                        # note that the 'overlap' relationship is not transitive
                        if ((exons[1, i] != exons[1, j]) and
                          (((exons[1, i] <= exons[1, j]) and (exons[1, i] >= exons[0,j])) or ((exons[1, j] <= exons[1, i]) and (exons[1, j] >= exons[0, i]))) and
                          (min(exons[1, i], exons[1, j]) - max(exons[1, i], exons[0, j]) >= MIN_OVERLAP)):

                            assert(not((exons[0, i] == exons[0, j]) and (exons[1, i] == exons[1, j])))
                            assert(exons[1, i] != exons[1, j])
                        
                            # add new events to the list
                            if not exons[1, i] in leftsites:
                                leftsites.append(exons[1, i])
                                leftidx.append(which_exons[i])
                            if not exons[1, j] in leftsites:
                                leftsites.append(exons[1, j])
                                leftidx.append(which_exons[j])

            # construct the output
            if len(leftsites) >= 2:
                if strand=='+'
                    exon_alt_5prime.append({'threeprimesite':exon_idx, 'fiveprimesites':leftidx})
                    idx_alt_5prime.append(ix) 
                if strand =='-'
                    exon_alt_3prime.append({'fiveprimesite':exon_idx, 'threeprimesites':leftidx})
                    idx_alt_3prime.append(ix)

        print '\nNumber of alternative 5 prime sites:\t\t\t\t%d' % len(idx_alt_5prime)
        print 'Number of alternative 3 prime sites:\t\t\t\t%d' % len(idx_alt_3prime)

        return (idx_alt_5prime, exon_alt_5prime, idx_alt_3prime, exon_alt_3prime)
