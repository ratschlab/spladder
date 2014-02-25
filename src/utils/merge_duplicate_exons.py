def merge_duplicate_exons(genes, CFG):
# genes = merge_duplicate_exons(genes, CFG)

    num_removed = 0
    for i in range(genes.shape[0]):
        if CFG['verbose'] and i % 100 == 0: 
            print '%i\r' % i

        ### check if there are non unique exons
        if uniqe_rows(genes[i].splicegraph.vertices.T).shape[0] == genes[i].splicegraph.shape[1]:
            continue

        exon_order = sp.argsort(genes[i].splicegraph.vertices[0, :])
        genes[i].splicegraph.vertices = genes[i].splicegraph.vertices[:, exon_order]
        genes[i].splicegraph.edges = genes[i].splicegraph.edges[exon_order, :][:, exon_order]
        genes[i].splicegraph.terminals = genes[i].splicegraph.terminals[:, exon_order]

        exons = genes[i].splicegraph.vertices.copy()
        admat = genes[i].splicegraph.edges.copy()
        initial = genes[i].splicegraph.terminals[0, :]
        terminal = genes[i].splicegraph.terminals[1, :]

        remove = sp.zeros((1, exons.shape[1]), dtype='bool')
        for j in range(exons.shape[1]):
            if remove[j]:
                continue
            idx = sp.where((exons[0, j] == exons[0, :]) & (exons[1, j] == exons[1, :]))[0]
            idx = sp.where(~sp.in1d(idx, j))[0]
            if idx.shape[0] > 0:
                remove[idx] = True
                for k in idx:
                    initial[j] = (initial[j] | initial[k])
                    terminal[j] = (terminal[j] | terminal[k])
                    admat[j, :] = (admat[j, :] | admat[k, :])
                    admat[:, j] = (admat[:, j] | admat[:, k])

        if sp.sum(remove) == 0:
            continue

        genes[i].splicegraph.edges = admat.copy()
        genes[i].splicegraph.terminals[0, :] = initial
        genes[i].splicegraph.terminals[1, :] = terminal

        k_idx = sp.where(~remove)[0]
        genes[i].splicegraph.vertices = genes[i].splicegraph.vertices[:, k_idx]
        genes[i].splicegraph.edges = genes[i].splicegraph.edges[kidx, :][:, k_idx]
        genes[i].splicegraph.terminals = genes[i].splicegraph.terminals[:, k_idx]
        num_removed += sum(remove)

    if CFG['verbose']:
        print '\n... removed %i duplicate exons ...' % num_removed
