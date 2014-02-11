# (author) Georg Zeller & Gunnar Raetsch, Mpi Tuebingen, Germany, 2009
# (author) Andre Kahles, MSKCC NYC, USA, 2013 

#function introns = get_intron_list(genes, CFG)
def get_intron_list(genes, CFG):

    ### form chunks for quick sorting
    chunks = sp.c_[sp.array(genes.chr_num, dtype = 'int'), sp.array(genes.strand, dtype = 'int'), sp.array(genes.start, dtype = 'int'), sp.array(genes.stop. dtyp = 'int')]
    (chunks, chunk_idx) = sort_rows(chunks, index=True)

    strands = ['+', '-']

    introns = sp.zeros((chunks.shape[0], 2), dtype = 'object')

    ### collect all possible combinations of contigs and strands
    regions = init_regions(CFG['bam_fnames'])
    chr_num = sp.array([x.chr_num for x in regions])
    keepidx = sp.where(sp.in1d(chr_num, sp.unique(sp.array([x.chr_num for x in genes]))))[0]
    regions = regions[keepidx]

    c = 0
    num_introns_filtered = 0

    for j in range(len(regions['chr'])):
        chr = regions[j].chr
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
                print >> sys.stdout, '%i (%i) genes done (%i introns taken)' % (c, chunks.shape[0], num_introns_filtered)

            gg = genes[chunk_idx[c]]
            gg.strand = strands[s]
            gg.start = max(gg.start - 5000, 1)
            gg.stop = gg.stop + 5000
            assert(gg.chr == chr)

            intron_list_tmp = add_reads_from_bam(gg, CFG['bam_fnames'], ['intron_list'], CFG['read_filter'], CFG['var_aware'])
            num_introns_filtered += intron_list_tmp[0].shape[0]
            introns[chunk_idx[c], s] = intron_list_tmp[0]
