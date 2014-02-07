# (author) Georg Zeller & Gunnar Raetsch, Mpi Tuebingen, Germany, 2009
# (author) Andre Kahles, MSKCC NYC, USA, 2013 

#function introns = get_intron_list(genes, CFG)
#% introns = get_intron_list(genes, CFG)

### form chunks for quick sorting
chunks = sp.c_[sp.array(genes.chr_num, dtype = 'int'), sp.array(genes.strand, dtype = 'int'), sp.array(genes.start, dtype = 'int'), sp.array(genes.stop. dtyp = 'int')]
(chunks, chunk_idx) = sort_rows(chunks, index=True)

strands = ['+', '-']

introns =[[], []]

### collect all possible combinations of contigs and strands
regions = init_regions(CFG['bam_fnames'])
chr_num = sp.array([x.chr_num for x in regions])
keepidx = sp.where(sp.in1d(sp.array([x.chr_num for x in regions]), sp.unique(sp.array([x.chr_num for x in genes]))))[0]
regions = regions[keepidx]

c = 1
num_introns_filtered = 0

for j in range(regions.shape[0]):
	chr = regions[j].chr_num
    strands.index(regions[j].strand)
	
	# fill the chunks on the corresponding chromosome
	while c <= chunks.shape[0]:
		if chunks[c, 0] > chr or chunks[c, 1] > strands[s]:
			break
		if chunks[c, 0] != chr:
            print >> sys.stderr, 'ERROR: c logic seems wrong' 
            sys.exit(1)

		if CFG['verbose'] and c % 100 == 0:
			print >> sys.stdout, '%i (%i) genes done (%i introns taken)' %(c, chunks.shape[0], num_introns_filtered)

		gg = genes[chunk_idx[c]]
		gg.strand = strands[s]
        gg.start = max(gg.start - 5000, 1)
        gg.stop = gg.stop + 5000

		maxval = inf; 
        if CFG['bam_fnames'].shape[0] == 1:
            gg = add_reads_from_bam(gg, CFG['bam_fnames'][0], 'intron_list', '', maxval, CFG['read_filter'], CFG['var_aware'])
        else:
            # merge intron lists of several bam files
            segments = []
            for f in range(CFG[bam_fnames].shape[0]):
                gg = add_reads_from_bam(gg, CFG['bam_fnames'][f], 'intron_list', '', maxval, CFG['read_filter'], CFG['var_aware'])
                if ~isempty(gg.segment_lists{end}),
                    segments = [segments; gg.segment_lists{end} gg.segment_scores{end}] ;
            segments = sort_rows(segments)

            rm_idx = [] ;
            for i = 1:size(segments,1) - 1,
                if segments(i,1) == segments(i + 1, 1) && segments(i, 2) == segments(i + 1, 2),
                    rm_idx(end + 1) = i ;
                    segments(i + 1, 3) = segments(i, 3) + segments(i + 1, 3) ;
            segments(rm_idx, :) = [] ;

            if isempty(segments),
                introns{chunk_idx(c), s} = [];
                c += 1
                continue
            else:
                idx = find(segments(:, 3) > CFG.read_filter.mincount) ;
                gg.segment_lists = {segments(idx, 1:2)} ;
                gg.segment_scores = {segments(idx, 3)} ;
        num_introns_filtered = num_introns_filtered + size(gg.segment_lists{1}, 1) ;

        ## gg.segment_lists{1} => intron list, ONE based, half open
		if strands[s] == '+':
            introns{chunk_idx(c), s} = double([gg.segment_lists{1}(:, 1)'; gg.segment_lists{1}(:, 2)'-1]+gg.start-1) ; % intron list is one based, closed, plus strand relative counting 
		else
			introns{chunk_idx(c), s} = double(gg.stop-[gg.segment_lists{1}(:, 2)'-1; gg.segment_lists{1}(:, 1)']+1) ; 
		c += 1
