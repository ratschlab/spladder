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

