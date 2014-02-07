def remove_short_exons(genes, CFG):
# [genes] = remove_short_exons(genes, terminal_short_extend, terminal_short_len, short_exon_skipped, short_exon_removed)

    short_exon_removed = 0
    short_exon_skipped = 0

    rm_idx = []
    for i in range(len(genes)):
        if CFG['verbose'] and i % 1000 == 0:
            print >> CFG['log_fd'], '%i' % i

        ### remove genes with no exons
        if genes[i].splicegraph.vertices.shape[0] == 0:
            rm_idx.append(i)
            continue
      
        ### extend terminal exons to terminal_short_extend if they are shorter than terminal_short_len
        if genes[i].splicegraph.vertices[1, 0] - genes[i].splicegraph.vertices[0, 0] < CFG['remove_exons']['terminal_short_len']:
            genes[i].splicegraph.vertices[0, 0] = genes[i].splicegraph.vertices[1, 0] - CFG['remove_exons']['terminal_short_extend']
            genes[i].start = min(genes[i].start, genes[i].splicegraph.vertices[0, 0])
      
        if genes[i].splicegraph.vertices[1, -1] - genes[i].splicegraph.vertices[0, -1] < CFG['remove_exons']['terminal_short_len']:
            genes[i].splicegraph.vertices[1, -1] = genes[i].splicegraph.vertices[1, -1] + CFG['remove_exons']['terminal_short_extend']
            genes[i].stop = max(genes[i].stop, genes[i].splicegraph.vertices[1, -1]) 
      
        ### check for very short exons and insert an edge that allows skipping them
        exons_remove_idx = []
        j = 1
        while (j <= genes[i].splicegraph.vertices.shape[1] - 1):
            if (genes[i].splicegraph.vertices[1, j] - genes[i].splicegraph.vertices[0, j]) < CFG['remove_exons']['min_exon_len']:
                foundp = False
                for jp in range(j + 1, genes[i].splicegraph.vertices.shape[1]):
                    if (genes[i].splicegraph.vertices[1, jp] - genes[i].splicegraph.vertices[0, jp] >= CFG['remove_exons']['min_exon_len_remove']) and (genes[i].splicegraph.edges[j, jp] == 1):
                        foundp = True
                        break 
                foundn = False
                for jn in range(j - 2, -1, -1):
                    if (genes[i].splicegraph.vertices[1, jn] - genes[i].splicegraph.vertices[0, jn] >= CFG['remove_exons']['min_exon_len_remove']) and (genes[i].splicegraph.edges[jn, j] == 1):
                        foundn = True
                        break
                if foundp and foundn:
                    genes[i].splicegraph.edges[jn, jp] = 1 
                    genes[i].splicegraph.edges[jp, jn] = 1 
      
                    if genes[i].splicegraph.vertices[1, j] - genes[i].splicegraph.vertices[0, j] < CFG['remove_exons']['min_exon_len_remove']:
                        short_exon_removed += 1
                        exons_remove_idx.append(j)
                    else:
                        short_exon_skipped += 1
            j += 1
      
        keep_idx = sp.where(~sp.in1d(sp.array(range(genes[i].splicegraph.vertices.shape[1])), exons_remove_idx))[0]
        genes[i].splicegraph.vertices = genes[i].splicegraph.vertices[:, keep_idx]
        genes[i].splicegraph{2} = genes[i].splicegraph.edges[keep_idx, :][:, keep_idx]
        genes[i].splicegraph{3} = genes[i].splicegraph.terminals[:, keep_idx]

    
    keep_idx = sp.where(~sp.in1d(sp.array(range(len(genes[i])), rm_idx))[0]
    genes = genes[keep_idx]

    if CFG['verbose']:
        print >> CFG['log_fd'], 'short_exon_removed: %i' % short_exon_removed
        print >> CFG['log_fd'], 'short_exon_skipped: %i' % short_exon_skipped

    return genes
