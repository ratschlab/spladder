function [genes, inserted] = insert_intron_edges(genes, CFG)

if not 'debug' in CFG:
    CFG['debug'] = False

print_intermediates = False

strands = ['+', '-']
P = []
both_missing = [0 0]
one_missing=[0 0]
multi = 0
next = 0
prev = 0

exon_vicinity_cnt1 = [0 0] 
exon_vicinity_cnt2 = [0 0] 
merge_idx = [] 
intron_tol = 1 

inserted = dict()
inserted['intron_in_exon'] = 0 
inserted['alt_53_prime'] = 0 
inserted['exon_skip'] = 0 
inserted['gene_merge'] = 0 
inserted['new_terminal_exon'] = 0 

num_unused_introns = sp.zeros((1, len(genes))

for i in range(len(genes)):
	if CFG['verbose'] and i % 1000 == 0:
        print >> CFG['fd_log'], '%i of %i genes' % (i, len(genes))

    s = strands.index(genes[i].strand)

    if not genes(i).introns[s]:
        continue

	unused_introns = []
    if CFG['debug']:
        print >> CFG['fd_log'], 'processing gene %i; with %i introns; time since last tic' % (i, len(genes[i].introns[s]))
        ### TODO timing

	for j in range(genes[i].introns[s].shape[1]):
		intron_used = False

        if ( j > 1 and genes[i].splicegraph.vertices.shape[1] > 1):
            genes[i].splicegraph.uniquify()

		### find exons within same gene whose end coincides with intron start
        idx1 = sp.where(sp.absolute(genes[i].splicegraph.vertices[1, :] - genes[i].introns[s][0, j] + 1) <= intron_tol)[0]
		### find exons within same gene whose start coincides with intron end
        idx2 = sp.where(sp.absolute(genes[i].splicegraph.vertices[0, :] - genes[i].introns[s][1, j] - 1) <= intron_tol)[0]

		### intron boundaries do not coincide with any exon boundaries
		if not idx1 and not idx2:
			both_missing[s] += 1

			if CFG['intron_edges']['insert_intron_retention']:
				### find all exons that completely include added introns 
				idx1__ = sp.where((genes[i].introns[s][0, j] > genes[i].splicegraph.vertices[0, :]) & (genes[i].introns[s][1, j] < vertices[1, :]))[0]
				for idx1_ in idx1__:

                    genes[i].splicegraph.vertices = sp.c_[genes[i].splicegraph.vertices, genes[i].splicegraph.vertices[:, idx1_]]
                    genes[i].splicegraph.vertices[1, -1] = genes[i].introns[s][0, j] - 1
							
                    genes[i].splicegraph.vertices = sp.c_[genes[i].splicegraph.vertices, genes[i].splicegraph.vertices[:, idx1_]]
                    genes[i].splicegraph.vertices[0, -1] = genes[i].introns[s][1, j] + 1
							
                    genes[i].splicegraph.new_edge()
                    adj_mat = sp.triu(genes[i].splicegraph.edges)
                    genes[i].splicegraph.edges[:, -1] = adj_mat[:, idx1_]    # incoming edges of idx1_
                    genes[i].splicegraph.edges[-1, :] = adj_mat[:, idx1_].T

                    genes[i].splicegraph.new_edge()
                    adj_mat = sp.triu(genes[i].splicegraph.edges)
                    genes[i].splicegraph.edges[:, -1] = adj_mat[idx1_, :].T    # outgoing edges of idx1_
                    genes[i].splicegraph.edges[-1, :] = adj_mat[idx1_, :]
                    genes[i].splicegraph.edges[-2, -1] = 1
                    genes[i].splicegraph.edges[-1, -2] = 1
							
                    genes[i].splicegraph.terminals = sp.c_[genes[i].splicegraph.terminals, genes[i].splicegraph.terminals[:, idx1_]
                    genes[i].splicegraph.terminals[1, -1] = 0 # cannot be an end
                    genes[i].splicegraph.terminals = sp.c_[genes[i].splicegraph.terminals, genes[i].splicegraph.terminals[:, idx1_]
                    genes[i].splicegraph.terminals[0, -1] = 0 # cannot be a start
							
					inserted['intron_in_exon'] += 1
                    assert(sp.all(genes[i].splicegraph.vertices[0, :] <= genes[i].splicegraph.vertices[1, :]))

                    if CFG['debug']:
                        print >> CFG['fd_log'], '%s\tintron_retention_exon\t%c\t%i\t%i\t%i\t%i\n', % (genes[i].chr, genes[i].strand, genes[i].splicegraph.vertices[0, -2],
                                                                                                      genes[i].splicegraph.vertices[1, -2], genes[i].splicegraph.vertices[0, -1], 
                                                                                                      genes[i].splicegraph.vertices[1, -1])
					intron_used = True

			if not intron_used:
                unused_introns.append(j)
			continue # with next intron

		# did not find exons in same gene sharing boundaries with intron start
		# find first end in previous gene on same strand
        if not idx1 and i > 0 and genes[i - 1].chr_num == genes[i].chr_num and genes[i - 1].strand == genes[i].strand: 
			### find all exon ends in previuos gene that coincide with intron start j
			idx1_ = find(abs(genes(i-1).splicegraph{1}(2,:) - genes(i).introns{s}(1,j) + 1) <= intron_tol) ;
            idx1_ = sp.where(sp.absolute(genes[i-1].splicegraph.vertices[1, :] - genes[i].introns[s][0, j] <= intron_tol))[0]
            if idx1_:
				prev += 1
				# mark the two genes for merging
				if CFG['intron_edges']['gene_merges']:
                    merge_idx = sp.c_[merge_idx, sp.array([i-1, i])]
					intron_used = True
                if not intron_used:
                    unused_introns.append(j)
				continue # with next intron

		# did not find exons in same gene sharing boundaries with intron end
		# find second end in next gene on same strand
        if not idx2 and i < len(genes) and genes[i + 1].chr_num == genes[i]. chr_num and genes[i+1].strand == genes[i].strand:
			### find all exon starts in following gene that coincide with intron end j
            idx2_ = sp.where(sp.asbolute(genes[i+1].splicegraph.vertices[0, :] - genes[i].introns[s][1, j] - 1) <= intron_tol)[0]
            if idx2_:
				next += 1
				# mark the two genes for merging
				if CFG['intron_edges']['gene_merges']: 
                    merge_idx = sp.c_[merge_idx, sp.array([i, i+1])]
					intron_used = True
                if not intron_used:
                    unused_introns.append(j)
				continue # with next intron

		# did not find exons in same gene sharing boundaries with intron start
		# check whether the intron starts in the vicinity of an exon
		if not idx1: 
			### find all exons that overlap intron-start j +/- CFG.intron_edges.vicinity_region
            idx1__ = sp.where((genes[i].splicegraph.vertices[0, :] - CFG['intron_edges']['vicinity_region'] <= genes[i].introns[s][0, j]) & 
                              (genes[i].splicegraph.vertices[1, :] + CFG['intron_edges']['vicinity_region'] >= genes[i].introns[s][0, j]))[0]

            ### check, if we can find an exon after the current intron and there is continuous coverage between intron end and exon
            if not idx1__:
                idx1__ = sp.argmax(genes[i].splicegraph.vertices[0, :] > genes[i].introns[s][1, j])
                if idx1__:
                    gg = genes[i]
                    gg.strand = strands[s]
                    gg.strands = strands[s]
                    gg.start = genes[i].introns[s][1, j] + 1  ### start of presumable exon
                    gg.stop = genes[i].splicegraph.vertices[1, idx1__] ### stop of next exon
                    maxval = inf; 
                    gg = add_reads_from_bam(gg, CFG.bam_fnames, 'exon_track', '', maxval, CFG.read_filter, CFG.var_aware);
                    if gg.strand == '-',
                        gg.tracks = gg.tracks[::-1]
                    ### TODO: make the following a configurable
                    if sp.mean(sp.sum(gg.tracks, axis = 0) > 10) < 0.9:
                        idx1__ = []

			# only take the case closest to an existing splice site
			diff1 = sp.absolute(genes[i].splicegraph.vertices[0, idx1__] - genes[i].introns[s][0, j])
			diff2 = sp.absolute(genes[i].splicegraph.vertices[1, idx1__] - genes[i].introns[s][0, j])
			diff = sp.minimum(diff1, diff2)
			idx1__ = idx1__(sp.argmin(diff))
			for idx1_ in idx1__:
				if genes[i].introns[s][0, j] - 1 - genes[i].splicegraph.vertices[0, idx1_] >= CFG['intron_edges']['min_exon_len']:
					exon_vicinity_cnt1[s] += 1
					genes[i].splicegraph.vertices = sp.c_[genes[i].splicegraph.vertices, genes[i].splicegraph.vertices[:, idx1_]]
					genes[i].splicegraph.vertices[1, -1] = genes[i].introns[s][0, j] - 1  # set exon end to intron start - 1
                    genes[i].splicegraph.new_edge()
					genes[i].splicegraph.edges[:, -1] = genes[i].splicegraph.edges[:, idx1_]
					genes[i].splicegraph.edges[-1, :] = genes[i].splicegraph.edges[idx1_, :]
					genes[i].splicegraph.terminals = sp.c_[genes[i].splicegraph.terminals, genes[i].splicegraph.terminals[:, idx1_]]  # copy from original exon
					genes[i].splicegraph.terminals[1, -1] = 0 # cannot be an end
								
                    assert(sp.all(genes[i].splicegraph.vertices[0, :] <= genes[i].splicegraph.vertices[1, :]))

					# check exons whose start coincides with intron end
					genes[i].splicegraph.add_intron(genes[i].splicegraph.edges.shape[0], 0, idx2, 1)
								
					inserted['alt_53_prime'] += 1
                    if CFG['debug']:
                        for idx2_ in idx2:
                            print >> CFG['fd_log'], '%s\talternative_53_prime1\t%c\t%i\t%i\t%i\n' % (genes[i].chr, genes[i].strand, genes[i].splicegraph,vertices[1, idx1_], 
                                                                                                     genes[i].splicegraph.vertices[1, -1], genes[i].splicegraph.vertices[0, idx2_])
					intron_used = True

			### if no proximal exon was found, insert new terminal exon, if wished
			if  not intron_used and CFG['intron_edges']['append_new_terminal_exons']:
				inserted['new_terminal_exon'] += 1

				iregion = sp.array([[genes[i].introns[s][0, j] - CFG['intron_edges']['append_new_terminal_exons_len']], [genes[i].introns[s][0, j] - 1]])
				idx_iregion = sp.where(genes[i].introns[s][1, :] >= iregion[0] & genes[i].introns[s][1, :] < iregion[1] - 1)[0]
                if idx_region:
                    if not idx_iregion.shape[0] == 1,
						idx_iregion = idx_iregion[sp.argmax(genes[i].introns[s][1, idx_iregion])]
					iregion[0] = genes[i].introns[s][1, idx_iregion] + 1
					assert(iregion[0] < iregion[1])

				genes[i].splicegraph.vertices = sp.c_[genes[i].splicegraph.vertices, iregion]
				genes[i].splicegraph.new_edge()
				genes[i].splicegraph.terminals = sp.c_[genes[i].splicegraph.terminals, sp.array([1, 0])] # can be a start, but cannot be an end

                for tmp_idx in idx2:
                    if genes[i].splicegraph.terminals[0, tmp_idx] == 1 and genes[i].introns[s][1, j] + 1 <= genes[i].splicegraph.vertices[1, tmp_idx]:
                        genes[i].splicegraph.vertices[0, tmp_idx] = genes[i].introns[s][1, j] + 1
                assert(sp.all(genes[i].splicegraph.vertices[1, :] >= genes[i].splicegraph.vertices[0, :]))

				genes[i].splicegraph.add_intron(idx2, 1, genes[i].splicegraph.vertices.shape-0], 0)
				intron_used = True

			if intron_used:
				continue

		# did not find exons in same gene sharing boundaries with intron end
		# check whether the intron ends in the vicinity of an exon
		if not idx2: 
			idx2__ = sp.where(genes[i].splicegraph.vertices[0, :] - CFG['intron_edges']['vicinity_region'] <= genes[i].introns[s][1, j] &
						      genes[i].splicegraph.vertices[1, :] + CFG['intron_edges']['vicinity_region'] >= genes[i].introns[s][1, j])[0]

            ### check, if we can find an exon after the current intron and there is continuous coverage between intron end and exon
            if not idx2__:
                idx2__ = sp.argmax(genes[i].splicegraph.vertices[0, :] > genes[i].introns[s][1, j])
                if idx2__:
                    gg = genes[i]
                    gg.strand = strands[s]
                    gg.strands = strands[s]
                    gg.start = genes[i].introns[s][1, j] + 1  ### start of presumable exon
                    gg.stop = genes[i].splicegraph.vertices[1, idx2__]  ### stop of next exon
                    maxval = inf; 
                    gg = add_reads_from_bam(gg, CFG.bam_fnames, 'exon_track', '', maxval, CFG.read_filter, CFG.var_aware);
                    if gg.strand == '-':
                        gg.tracks = gg.tracks[::-1]
                    ### TODO: make configurable
                    if sp.mean(sp.sum(gg.tracks, axis=1) > 10) < 0.9:
                        idx2__ = []

			# only take the case closest to an existing splice site
			diff1 = sp.absolute(genes[i].splicegraph.vertices[0, idx2__] - genes[i].introns[s][1, j])
			diff2 = sp.absolute(genes[i].splicegraph.vertices[1, idx2__] - genes[i].introns[s][1, j])
			diff = sp.minimum(diff1, diff2) ;
			idx2__ = idx2__(sp.argmin(diff)) 
			for idx2_ in idx2__:
				if genes[i].splicegraph.vertices[1, idx2_] - genes[i].introns[s][1, j] >= CFG['intron_edges']['min_exon_len']:
					exon_vicinity_cnt2[s] += 1
					genes[i].splicegraph.vertices = sp.c_[genes[i].splicegraph.vertices, genes[i].splicegraph.vertices[:, idx2_]]
					genes[i].splicegraph.vertcies[0, -1] = genes[i].introns[s][1, j] + 1
                    genes[i].splicegraph.new_edge()
					genes[i].splicegraph.edges[:, -1] = genes[i].splicegraph.edges[:, idx2_]
					genes[i].splicegraph.edges[-1, :] = genes[i].splicegraph.edges[idx2_, :]
					genes[i].splicegraph.terminals = sp.c_[genes[i].splicegraph.terminals, genes[i].splicegraph.terminals[:, idx2_]   # copy from original exon
					genes[i].splicegraph.terminals[0, -1] = 0   # cannot be a start
						
                    assert(sp.all(genes[i].splicegraph.vertices[0, :] <= genes[i].splicegraph.vertices[1, :]))

					genes[i].splicegraph.add_intron(idx1, 1, genes[i].splicegraph.edges.shape[0], 0)
					
					inserted['alt_53_prime'] += 1

                    if CFG['debug']:
                        for idx1_ in idx1:
                            print >> CFG['fd_log'], '%s\talternative_53_prime2\t%c\t%i\t%i\t%i\n', % (genes[i].chr, genes[i].strand, genes[i].splicegraph.vertices[1, idx1_], 
                                                                                                      genes[i].splicegraph.vertices[0, -1], genes[i].splicegraph.vertices[0, idx2_])
					intron_used = True

			### if no proximal exon was found, insert new terminal exon, if wished
			if not intron_used and CFG['intron_edges']['append_new_terminal_exons']:

                ### define range of new exon
				iregion = sp.array([[genes[i].introns[s][1, j] + 1], [genes[i].introns[s][1, j] + CFG['intron_edges']['append_new_terminal_exons_len']]])
                ### find introns starting within new exon
				idx_iregion = sp.where(genes[i].introns[s][0, :] > iregion[0] + 1 & genes[i].introns[s][0, :] <= iregion[1])[0]

				if idx_iregion:
					if not idx_iregion.shape[0] == 1: 
						idx_iregion = idx_iregion[sp.argmin(genes[i].introns[s][0, idx_iregion])]
                    ### let new exon end at position before next intron starts
					iregion[1] = genes[i].introns[s][0, idx_iregion] - 1 
					assert(iregion[0] < iregion[1])

				inserted['new_terminal_exon'] += 1
				genes[i].splicegraph.vertices = sp.c_[genes[i].splicegraph.vertices, iregion]
				genes[i].splicegraph.new_edge()
				genes[i].splicegraph.terminals = sp.c_[genes[i].splicegraph.terminals, sp.array([0, 1])]  # cannot be a start but can be an end

                ### adapt terminal exon ends if the new intron starts within them
                for tmp_idx in idx1:
                    if genes[i].splicegraph.terminals[1, tmp_idx] and genes[i].introns[s][0, j] - 1 >= genes[i].splicegraph.vertices[0, tmp_idx]:
                        genes(i).splicegraph.vertices[1, tmp_idx] = genes[i].introns[s][0, j] - 1
                assert(sp.all(genes[i].splicegraph.vertices[1, :] >= genes[i].splicegraph.vertices[0, :]))
						
				genes[i].splicegraph.add_intron(idx1, 1, genes[i].splicegraph.edges.shape[0], 0)
				intron_used = True

			if intron_used:
				continue
		
		if not idx1 or not idx2,
			one_missing[s] += 1
			if not intron_used:
                unused_introns.append(j)
			continue

        ### TODO: hard coded limit
		if idx1.shape[0] > 20 or idx2.shape[0] > 20:
			multi += 1
			if not intron_used: 
                unused_introns.append(j)
			continue
		
		### both idx1 and idx2 are not empty and are both shorter than 4
		### insert exon skips
		for idx1_ in idx1:
			for idx2_ in idx2:
				if genes[i].splicegraph.edges[idx1_, idx2_] == 0:
					inserted['exon_skip'] += 1
								
					#adj_mat = triu(genes(i).splicegraph{2}) ;
					#id1 = find(adj_mat(idx1_,:)) ;
					#if length(id1)==1 && adj_mat(id1, idx2_),
					#		fprintf(CFG.fd_log, '%s\texon_skip\t%c\t%i\t%i\t%i\t%i\t%i\t%i\n', genes(i).chr, genes(i).strand, genes(i).splicegraph{1}(1,idx1_), genes(i).splicegraph{1}(2,idx1_), genes(i).splicegraph{1}(1,id1),	...
					#						genes(i).splicegraph{1}(2,id1), genes(i).splicegraph{1}(1,idx2_), genes(i).splicegraph{1}(2,idx2_)) ;

		genes[i].splicegraph.add_intron(idx1, 1, idx2, 1)
		used_intron = True

		#for i1=idx1,
		#	for i2=idx2,
		#		P(end+1,:)=[i1,i2] ;
		#		genes(i).splicegraph{2}(i1,i2)=1 ;
		#		genes(i).splicegraph{2}(i2,i1)=1 ;

	idx_unused = sp.where((genes[i].introns[s][1, unused_introns] >= genes[i].start) & (genes[i].introns[s][0, unused_introns] <= genes[i].stop))[0]
	unused_introns = unused_introns[idx_unused]
	if unused_introns:
		print 'Warning: unused introns: %s' % str(unused_introns)
	num_unused_introns[i] += unused_introns.shape[0]

if print_intermediates:
    print 'one missing: %s' % str(one_missing)
    print 'multi: %s' % str(multi)
    print 'num_unused_introns: $i' % sum(num_unused_introns)

merge_idx = unique_rows(merge_idx)
rm_map = sp.zeros((1, leng(genes)))

for i in range(merge_idx.shape[0]):
	
	while rm_map[merge_idx[i, 0]] == 1: 
		merge_idx[i, 0] -= 1

	# merge transcripts
	for j in range(genes[merge_idx[i, 1]].exons.shape[0]):
		genes[merge_idx[i, 0]].transcripts.append(genes[merge_idx[i, 1]].transcripts[j])
		genes[merge_idx[i, 0]].exons.append(genes[merge_idx[i, 1]].exons[j]
	# merge intron lists
	for k in range(genes[merge_idx[i, 1]].introns.shape[0]):
		for j in range(genes[merge_idx[i, 1]].introns[k].shape[1]):
			genes[merge_idx[i, 0]].introns{k}(:, end + 1) = genes(merge_idx(i, 2)).introns{k}(:, j) ;
	% merge splice graphs
	genes(merge_idx(i, 1)).splicegraph{1} = [genes(merge_idx(i,1)).splicegraph{1} genes(merge_idx(i,2)).splicegraph{1}] ;
	m = size(genes(merge_idx(i, 1)).splicegraph{2}, 1) ;
	n = size(genes(merge_idx(i, 2)).splicegraph{2}, 1) ;
	genes(merge_idx(i, 1)).splicegraph{2}(m + 1 : n + m, m + 1 : n + m) = genes(merge_idx(i, 2)).splicegraph{2} ;
	genes(merge_idx(i, 1)).splicegraph{3} = [genes(merge_idx(i, 1)).splicegraph{3} genes(merge_idx(i, 2)).splicegraph{3}] ;

	% extend start/stop
	genes(merge_idx(i, 1)).start = min(genes(merge_idx(i, 1)).start, genes(merge_idx(i, 2)).start) ;
	genes(merge_idx(i, 1)).stop = max(genes(merge_idx(i, 1)).stop, genes(merge_idx(i, 2)).stop) ;

	%genes(merge_idx(i,1))=build_splice_graph_caller(genes(merge_idx(i,1))) ;
	%genes(merge_idx(i,1))=infer_splice_graph_caller(genes(merge_idx(i,1))) ;

	rm_map(merge_idx(i, 2)) = 1 ;

	inserted.gene_merge = inserted.gene_merge + 1 ;
end ;
genes(rm_map == 1) = [] ;

%size(P,1)/(size(P,1)+both_missing+one_missing)
%both_missing/(size(P,1)+both_missing+one_missing)
%one_missing/(size(P,1)+both_missing+one_missing)

if ( size(genes(i).splicegraph{1},2) > 1) 
    genes(i) = uniquify_splicegraph(genes(i));
end ;

for i = 1:length(genes),
    assert(all(genes(i).splicegraph{1}(1,:) <= genes(i).splicegraph{1}(2, :))) ;
end ;

for ix = 1:length(genes)
	[dummy,exon_order] = sort(genes(ix).splicegraph{1}(1,:),2,'ascend');
	genes(ix).splicegraph{1} = genes(ix).splicegraph{1}(:, exon_order);
	genes(ix).splicegraph{2} = genes(ix).splicegraph{2}(exon_order, exon_order);
	genes(ix).splicegraph{3} = genes(ix).splicegraph{3}(:, exon_order);
end ;

