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
