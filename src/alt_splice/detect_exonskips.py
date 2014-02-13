def detect_exonskips(genes, idx_alt):
# [idx_exon_skips, exon_exon_skips] = detect_exonskips(genes, idx_alt) ;

idx_exon_skips = []
exon_exon_skips = []
for ix in idx_alt:
    if ix % 50 == 0:
        print '.',
    num_exons = genes[ix].splicegraph.get_len()
    edges = genes[ix].splicegraph.edges
    for exon_idx in range(num_exons - 2):
        for exon_idx1 in range(exon_idx + 1, num_exons - 1):
            for exon_idx2 in range(exon_idx1 + 1, num_exons):
                if (edges[exon_idx, exon_idx1] == 1) and edges[exon_idx, exon_idx2] and edges[exon_idx1, exon_idx2]:
                    idx_exon_skips.append(ix)
                    exon_exon_skips.append([exon_idx, exon_idx1, exon_idx2])
print '\nNumber of single exon skips:\t\t\t\t\t%d', len(idx_exon_skips))
return (idx_exon_skips, exon_exon_skips)


