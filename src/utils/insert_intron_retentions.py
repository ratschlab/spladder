import scipy as sp
import scipy.linalg as spla

# written by Georg Zeller & Gunnar Raetsch, Mpi Tuebingen, Germany, 2009, Andre Kahles, MSKCC, 2013
# [genes, inserted] = insert_intron_retentions(genes, CFG)
def insert_intron_retentions(genes, CFG):

inserted['intron_retention'] = 0

### form chunks for quick sorting
chunks = sp.c_[sp.array(genes.chr_num, dtype = 'int'), sp.array(genes.strand, dtype = 'int'), sp.array(genes.start, dtype = 'int'), sp.array(genes.stop. dtyp = 'int')]
(chunks, chunk_idx) = sort_rows(chunks, index=True)

strands = ['+', '-']

### form all possible combinations of contigs and strands --> regions
regions = init_regions(CFG['bam_fnames'])
### keep only chromosomes found in genes
keepidx = sp.where(sp.in1d(sp.arra([x.chr_num for x in regions]), sp.unique(sp.array([x.chr_num for x in genes]))))[0]
regions = regions[keepidx]

c = 0 
num_introns_added = 0
num_introns = 0
for j in range(regions.shape[0]):
	chr_num = regions[j].chr_num
	s = strands.index(regions[j].strand)
	
	# fill the chunks on the corresponding chromosome
	while c <= chunks.shape[0]:
		if chunks[c, 0] > chr_num or chunks[c, 1] > strands[s]:
			break
		if chunks[c, 0] != chr_num:
            print >> sys.stderr, 'ERROR: c logic seems wrong'
            sys.exit(1)

        if CFG['verbose'] and c % 100 == 0:
			print >> sys.stdout, '\r %i(%i) genes done (found %i new retentions in %i tested introns, %2.1f%%)' % (c, chunks.shape[0], num_introns_added, num_introns, 100 * num_introns_added / float(num_introns))

        gg = genes[chunk_idx[c]]
        gg.strand = strands[s]
        tracks = add_reads_from_bam(gg, CFG['bam_fnames'], ['exon_track'], CFG['read_filter'], CFG['var_aware'])

        exon_coverage = sp.zeros((gg.splicegraph.vertices.shape[1],), dtype = 'int')
		for k in range(gg.splicegraph.vertices.shape[1]):
			idx = sp.arange(gg.splicegraph.vertices[0, k], gg.splicegraph.vertices[1, k]) - gg.start
			exon_coverage[k] = sp.median(sp.sum(tracks[:, idx], axis=0), axis=1) # median coverage for exon k

		### check for all vertex-pairs, if respective intron can be retained
		new_retention = zeros(size(gg.splicegraph{2})) ;
        for k in range(gg.splicegraph.edges.shape[0]):
            for l = range(k + 1, gg.splicegraph.edges.shape[0]):
                if gg.splicegraph.edges[k, l] == 1:
					num_introns += 1
					idx = sp.arange(gg.splicegraph.vertices[1, k], gg.splicegraph.vertices[0, l]) - gg.start
					icov = sp.sum(tracks[:, idx], axis=0) 
					if median(icov) > CFG['intron_retention']['min_retention_cov'] and
						sp.mean(icov > (0.5 * sp.mean(icov))) > CFG['intron_retention']['.in_retention_region'] and  # fraction of covered positions
						max(exon_coverage[k], exon_coverage[l]) / (1e-6 + min(exon_coverage[k], exon_coverage[l])) <= CFG['intron_retention']['min_retention_max_exon_fold_diff'] and
						sp.mean(icov) >= CFG['intron_retention']['min_retention_rel_cov'] * (exon_coverage[k] + exon_coverage[l]) / 2 and
						sp.mean(icov) <= CFG['intron_retention']['max_retention_rel_cov'] * (exon_coverage[k] + exon_coverage[l]) / 2:

						new_retention[k, l] = 1
					#	fprintf(log_fd, '%s\tintron_retention\t%c\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%2.1f\n', gg.chr, gg.strand, gg.splicegraph{1}(1,k), gg.splicegraph{1}(2,k), gg.splicegraph{1}(1,l), gg.splicegraph{1}(2,l), ...
					#			floor(median(icov(1,:)+icov(2,:))), floor(gg.exon_coverage(k)), floor(gg.exon_coverage(l)), 100*mean(icov(1,:)+icov(2,:)>0)) ;
						inserted['intron_retention'] += 1
		any_added = False
		if sp.sum(new_retention.ravel()) > 0:
			new_retention = spla.expm(new_retention)
			while True:
				any_added = False
				for k in range(new_retention.shape[1]):
					for l = range(k + 1, new_retention.shape[1]):
						if new_retention[k, l] > 0:
							gg.splicegraph.add_intron_retention(k, l)
							new_retention = sp.c_[new_retention, sp.zeros((new_retention.shape[0],))]
							new_retention = sp.r_[new_retention, sp.zeros((new_retention.shape[1],))]
							new_retention[k, l] = 0
							any_added = True
							num_introns_added += 1
							#fprintf(log_fd, '%s\tintron_retention\t%i\t%i\t%i\t%i\t%i\t%2.1f\n', gg.chr, gg.splicegraph{1}(2,k), gg.splicegraph{1}(1,l), floor(median(icov(1,:)+icov(2,:))), gg.exon_coverage(k), gg.exon_coverage(l), 100*mean(icov(1,:)+icov(2,:)>0)) ;
							break
					if any_added:
                        break
                exon_order = sp.argsort(gg.splicegraph.vertices[0, :])
				gg.splicegraph.reorder(exon_order)
				new_retention = new_retention[exon_order, :][:, exon_order]
				if not any_added:
                    break
		if any_added:
            exon_order = sp.argsort(gg.splicegraph.vertices[0, :])
            gg.splicegraph.reorder(exon_order)
		genes[chunk_idx[c]] = gg
		c += 1
