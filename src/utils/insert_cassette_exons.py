# written by Andre Kahles, Mpi Tuebingen, Germany, 2012

def insert_cassette_exons(genes, CFG):
    # [genes, inserted] = insert_cassette_exons(genes, CFG)

    inserted['cassette_exon'] = 0

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
    num_exons_added = 0
    num_exons = 0

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
                print '\r %i(%i) genes done (found %i new cassette exons in %i tested intron pairs, %2.1f%%)' % (c, chunks.shape[0], num_exons_added, num_exons, 100*num_exons_added/float(num_exons))

            gg = genes[chunk_idx[c]]
            gg.strand = strands[s]
            tracks = add_reads_from_bam(gg, CFG['bam_fnames'], ['exon_track'], CFG['read_filter'], CFG['var_aware'])

            ### add introns implied by splicegraph to the list
            all_introns = gg.introns[s]
            for k in range(gg.splicegraph.edges.shape[0]):
                for l = range(k+1, gg.splicegraph.edges.shape[0]):
                    if gg.splicegraph.edges[k, l] == 1:
                        all_introns = sp.r_[all_introns, sp.array([gg.splicegraph.vertices[1, k], gg.splicegraph.vertices[0, l]])] # introns are half open
            all_introns = unique_rows(all_introns)
       
            ### use only relevant introns (inside gene boundaries)
            if all_introns.shape[0] > 0:
                keep_idx = sp.where((all_introns[:, 1] > gg.start) & (all_introns[:, 0] < gg.stop))[0]
                all_introns = all_introns[keep_idx, :]

            segment_starts = sp.sort(sp.unique(all_introns[:, 0]))
            segment_ends = sp.sort(sp.unique(all_introns[:, 1]))

            ### check for all intron-pairs, if exon could exist between them
            new_cassette = sp.zeros((all_introns.shape[0],)) 
            for k in range(all_introns.shape[0]):
                for l in range(k + 1, all_introns.shape[0]):
                    if all_introns[k, 1] >= all_introns[l, 0]:
                        continue
                    ### only take intron pair, if outer ends are supported by current exons
                    if (not all_introns[k, 0] in gg.splicegraph.vertices[1, :]) or (not all_introns[l, 1] in gg.splicegraph.vertices[0, :]): 
                        continue
                    curr_exon = [all_introns[k, 1], all_introns[l, 0]]
                    ### do not allow curr_exon to overlap existing exon
                    if sp.sum((gg.splicegraph.vertices[0, :] < curr_exon[1]) & (gg.splicegraph.vertices[1, :] > curr_exon[0])) == 0:
                        continue

                    num_exons += 1

                    if not ismember(curr_exon, gg.splicegraph.vertices.T, rows=True):
                        idx = sp.arange(curr_exon[0], curr_exon[1]) - gg.start
                        exon_cov = sp.sum(tracks[:, idx], axis=0)

                        pre_segment_end = sp.where(segment_ends < curr_exon[0])[0].max()
                        if pre_segment_end.shape[0] > 0:
                            pre_segment_cov = sp.sum(tracks[:, sp.arange(segment_ends[pre_segment_end], curr_exon[0]) - gg.start], axis=0)
                        else:
                            pre_segment_cov = sp.sum(tracks[:, sp.arange(curr_exon[0] - gg.start)], axis=0)
                        min_len_pre = min(pre_segment_cov.shape[0], exon_cov.shape[0])

                        aft_segment_start = sp.where(segment_starts > curr_exon[1])[0].min()
                        if aft_segment_start.shape[0] > 0:
                            aft_segment_cov = sp.sum(tracks[:, sp.arange(curr_exon[1], segment_starts[aft_segment_start]) - gg.start], axis=0)
                        else:
                            aft_segment_cov = sp.sum(tracks[:, (curr_exon[1] - gg.start):], axis=0)
                        min_len_aft = min(aft_segment_cov.shape[0], exon_cov.shape[0])

                        if sp.mean(exon_cov > (0.2 * sp.mean(exon_cov))) > CFG['cassette_exon']['min_cassette_region'] and
                           sp.median(exon_cov) > CFG['cassette_exon']['min_cassette_cov'] and
                           max(sp.median(exon_cov[-min_len_aft:]), sp.median(aft_segment_cov[:min_len_aft])) / min(sp.median(exon_cov[-min_len_aft:]), sp.median(aft_segment_cov[:min_len_aft])) - 1 >= CFG['cassette_exon']['min_cassette_rel_diff'] and
                           max(sp.median(exon_cov[:min_len_pre]), sp.median(pre_segment_cov[-min_len_pre:])) / min(sp.median(exon_cov[:min_len_pre]), sp.median(pre_segment_cov[-min_len_pre:])) - 1 >= CFG['cassette_exon']['min_cassette_rel_diff']:
                            new_cassette[k, l] = 1
                            inserted['cassette_exon'] += 1 
            any_added = False
            if any(new_cassette.ravel()):
                curr_sg = gg.splicegraph.vertices
                for k in range(new_cassette.shape[1]):
                    for l = range(k + 1, new_cassette.shape[1]):
                        if new_cassette[k, l] > 0:
                            exons_pre = sp.where(curr_sg[1, :] == all_introns[k, 0])[0]
                            exons_aft = sp.where(curr_sg[0, :] == all_introns[l, 1])[0]

                            gg.splicegraph.add_cassette_exon(sp.array([all_introns[k, 1], all_introns[l, 0]]), exons_pre, exons_aft)
                            new_cassette[k, l] = 0
                            any_added = True
                            num_exons_added += 1
                exon_order = sp.argsort(gg.splicegraph.vertices[0, :])
                gg.splicegraph.reorder(exon_order)
                if not any_added:
                    break
            if any_added:
                exon_order = sp.argsort(gg.splicegraph.vertices[0, :])
                gg.splicegraph.reorder(exon_order)
            ### clean up gene structure
            genes[chunk_idx[c]] = gg
            c += 1
