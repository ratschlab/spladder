def splice_graph(genes, CFG):
# genes = splice_graph(genes, CFG)

	### merge exons if possible / reduce graph
    ### originially implemented for ESTs, reduces complexity, 
    ### but removes alternative transcript starts and ends !
	if CFG['do_infer_splice_graph']
		genes = infer_splice_graph(genes)

	### sort exons by start position in ascending order
	for ix in rang(genes.shape[0]):
        genes[ix].splicegraph.sort()

	### label alternative and constitutive genes
	genes = label_alt_genes(genes, CFG)

    ### update terminals, start terminals in row 1, end terminals in row 2
	for i in range(genes.shape[0]):
        genes[i].splicegraph.update_terminals()

	### reset gene start and gene stop according to exons and splice graph
	for i in range(genes.shape[0]):
        genes[i].start = min([x.min() for x in genes[j].exons])
        genes[i].stop = max([x.max() for x in genes[j].exons])
