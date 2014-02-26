import scipy as sp
import scipy.linalg as spla

def reduce_splice_graph(genes):
    # genes = reduce_splice_graph(genes) ;
    #
    # Iterates over all genes and removes complexity
    # from the splice graph by:
    # - collapsing identical exons into a single node
    # - merging one single exon transcripts containing another 
    #   single exon transcript into one
    # - collapsing alternative transcript starts into a single
    #   start if they share 3' exon boundary
    # - collapsing alternative transcript ends into a single
    #   end if they share 5' exon boundary
    #
    for gene_idx in range(genes.shape[0]):
      
        vertices = genes[gene_idx].splicegraph.vertices
        edges = genes[gene_idx].splicegraph.edges
      
        if gene_idx % 1000 == 0:
            print '%d' % gene_idx
      
        ### no vertices in the splice graph
        if vertices.shape[0] == 0:
          continue
      
        ### find all the intron locations
        exon_order = sp.argsort(vertices[0, :])
        vertices = vertices[:, exon_order]
        edges = edges[exon_order, :][:, exon_order
        intron_loc = []
        for ix1 in range(vertices.shape[1] - 1):
            for ix2 in range(1, vertices.shape[1]):
                if edges[ix1, ix2] == 1:
                    intron_loc.extend(range(vertices[1, ix1] : vertices[0, ix2]))
        intron_loc = sp.unique(intron_loc)
        
        ### if one or two vertices (one or no edge) exist 
        if edges.shape[0] < 2:
            changed = True
            while changed:
                changed = False
                exon_idx = 0
                while exon_idx <= vertices.shape[1]:
                    test_exon_idx = exon_idx + 1
                    while test_exon_idx <= vertices.shape[1]:
                        ###  <------- exon_idx -------> 
                        ###    <-- test_exon_idx -->
                        if (vertices[0, exon_idx] <= vertices[0, test_exon_idx]) and (vertices[1, exon_idx] >= vertices[1, test_exon_idx]):
      
                             ### keep longer exon
                             new_index = sp.r_[range(test_exon_idx), range(test_exon_idx + 1, vertices.shape[1])]
                             vertices = vertices[:, new_index]
                  
                             changed = True
                        test_exon_idx += 1
                    exon_idx += 1
        ### more than two vertices exist
        else:
            changed = True
            while changed:
                changed = False
                exon_idx = 0
                while exon_idx <= vertices.shape[1]:
                    ## re-sort vertices if necessary
                    exon_order = sp.argsort(vertices[0, :])
                    vertices = vertices[:, exon_order]
                    edges = edges[exon_order, :][:, exon_order]
              
                    test_exon_idx = exon_idx + 1
                    while test_exon_idx <= vertices.shape[1]:
                        reduce_now = False
      
                        ### count incoming and outgoing edges for exon and test_exon
                        cur_edge_left = (edges[:exon_idx+1, exon_idx].sum() > 0)
                        test_edge_left = (edges[:test_exon_idx+1, test_exon_idx].sum() > 0)
                        cur_edge_right = (edges[exon_idx:, exon_idx].sum() > 0)
                        test_edge_right = (edges[test_exon_idx:, test_exon_idx].sum() > 0)
                  
                        ### 0000
                        ### no incoming or outgoing edges in exon and test_exon
                        if not cur_edge_left and not cur_edge_right and not test_edge_left and not test_edge_right:
                            ###     <------ exon ------->>>>>>>>>>>>>  OR          <------ exon ------>
                            ###              <--- test_exon ---->           <---- test_exon ---->>>>>>>>>>>>
                            if ((vertices[1, exon_idx] >= vertices[0, test_exon_idx]) and (vertices[0, exon_idx] <= vertices[0, test_exon_idx])) or
                               ((vertices[1, test_exon_idx] >= vertices[0, exon_idx]) and (vertices[0, test_exon_idx] <= vertices[0, exon_idx])) and 
                               (sp.sum(sp.in1d(range(min(vertices[0, exon_idx], vertices[0, test_exon_idx]), max(vertices[1, exon_idx], vertices[1, test_exon_idx])), intron_loc)) == 0):
      
                                ### merge exons if they overlap and they do not span any intronic position
                                vertices[0, exon_idx] = min(vertices[0, exon_idx], vertices[0, test_exon_idx])
                                vertices[1, exon_idx] = max(vertices[1, exon_idx], vertices[1, test_exon_idx])
                                new_index = sp.r_[range(test_exon_idx), range(test_exon_idx + 1, vertices.shape[1])]
                              
                                vertices = vertices[:, new_index]
                                edges = edges[new_index, :][:, new_index] # no need to combine any adges, as both exons have a degree of 0
                              
                                reduce_now = True
                                changed = True
                      
                        ### 0101
                        ### outgoing edges in exon and test_exo, no incoming edges
                        elif not cur_edge_left and cur_edge_right and  not test_edge_left and test_edge_right:
                            ###   ----- exon -----<
                            ###   --- test_exon --<
                            if (vertices[1, exon_idx] == vertices[1, test_exon_idx]) and
                               (sp.sum(sp.in1d(range(min(vertices[0, exon_idx], vertices[0, test_exon_idx]), vertices[1, exon_idx]), intron_loc)) == 0):
                                
                                ### merge exons if they share the same right boundary and do not span intronic positions
                                vertices[0, exon_idx] = min(vertices[0, exon_idx], vertices[0, test_exon_idx])
                                new_index = sp.r_[range(test_exon_idx), range(test_exon_idx + 1, vertices.shape[1])]
                              
                                vertices = vertices[:, new_index]
                                edges[exon_idx, :] = edges[exon_idx, :] | edges[test_exon_idx, :]
                                edges[:, exon_idx] = edges[:, exon_idx] | edges[:, test_exon_idx]
                                edges = edges[new_index, :][:, new_index]
                              
                                reduce_now = True
                                changed = True
      
                        ### 1010
                        ### incoming edges in exon and test_exon, no outgoing edges
                        elif cur_edge_left and not cur_edge_right and test_edge_left and not test_edge_right:
                            ###   >---- exon ------
                            ###   >-- test_exon ---
                            if (vertices[0, exon_idx] == vertices[0, test_exon_idx]) and
                               (sp.sum(sp.in1d(range(vertices[0, exon_idx], max(vertices[1, exon_idx], vertices[1, test_exon_idx])), intron_loc)) == 0):
                                
                                ### merge exons if they share the same left boundary and do not span intronic positions
                                vertices[1, exon_idx] = max(vertices[1, exon_idx], vertices[1, test_exon_idx])
                                new_index = sp.r_[range(test_exon_idx), range(test_exon_idx + 1, vertices.shape[1])]
                              
                                vertices = vertices[:, new_index]
                                edges[exon_idx, :] = edges[exon_idx, :] | edges[test_exon_idx, :]
                                edges[:, exon_idx] = edges[:, exon_idx] | edges[:, test_exon_idx]
                                edges = edges[new_index, :][:, new_index]
                              
                                reduce_now = True
                                changed = True
      
                        ### 1111
                        ### exon and test_exon have both incoming and outgoing edges
                        elif cur_edge_left and cur_edge_right and test_edge_left and test_edge_right:
                            ###  >------ exon -----<
                            ###  >--- test_exon ---<
                            if (vertices[0, exon_idx] == vertices[0, test_exon_idx]) and (vertices[1, exon_idx] == vertices[1, test_exon_idx]):
                              
                                ### collapse identical exons into one node
                                new_index = sp.r_[range(test_exon_idx), range(test_exon_idx + 1, vertices.shape[1])]
                              
                                vertices = vertices[:, new_index]
                                edges[exon_idx, :] = edges[exon_idx, :] | edges[test_exon_idx, :]
                                edges[:, exon_idx] = edges[:, exon_idx] | edges[:, test_exon_idx]
                                edges = edges[new_index, :][:, new_index]
                                
                                reduce_now = True
                                changed = True
      
                        if not reduce_now:
                            test_exon_idx += 1
                    exon_idx += 1
        genes[gene_idx].splicegraph.vertices = vertices
        genes[gene_idx].splicegraph.edges = edges

    return genes


def filter_by_edgecount(genes, CFG):

    ### filter splicegraphs by support count over samples
    for i in len(genes):
        k_idx = sp.where(genes[i].splicegraph.edges.sum(axis = 1) == 0)[0]
        genes[i].splicegraph.edges = (genes[i].edge_count >= CFG['sg_min_edge_count'])
        ### remove all exons that have no incoming or outgoing edges (caused by validation, keep single exon transcripts that occured before)
        k_idx2 = sp.where(genes[i].splicegraph.edges.sum(axis = 1) == 0)[0]
        rm_idx = sp.where(~sp.in1d(k_idx2, k_idx))[0]
        keep_idx = sp.where(~sp.in1d(sp.array(range(genes[i].splicegraph.edges.shape[0])), rm_idx))[0]
        if keep_idx.shape[0] > 0:
            genes[i].subset(keep_idx)
        else:
            genes = []

    return genes


def insert_intron_retentions(genes, CFG):
    # written by Georg Zeller & Gunnar Raetsch, Mpi Tuebingen, Germany, 2009, Andre Kahles, MSKCC, 2013
    # [genes, inserted] = insert_intron_retentions(genes, CFG)

    inserted = 0

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
                            inserted += 1
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
    return (genes, inserted)

def insert_intron_edges(genes, CFG):
    #function [genes, inserted] = insert_intron_edges(genes, CFG)

    if not 'debug' in CFG:
        CFG['debug'] = False

    print_intermediates = False

    strands = ['+', '-']
    P = []
    both_missing = [0 0]
    one_missing = [0 0]
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

    num_unused_introns = sp.zeros((1, genes.shape[0])

    for i in range(genes.shape[0]):
        if CFG['verbose'] and i % 1000 == 0:
            print >> CFG['fd_log'], '%i of %i genes' % (i, genes.shape[0])

        s = strands.index(genes[i].strand)

        if not genes[i].introns[s]:
            continue

        unused_introns = []
        if CFG['debug']:
            print >> CFG['fd_log'], 'processing gene %i; with %i introns; time since last tic' % (i, len(genes[i].introns[s]))
            ### TODO timing

        for j in range(len(genes[i].introns[s])):
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
    rm_map = sp.zeros((1, len(genes)))

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
                genes[merge_idx[i, 0]].introns[k] = sp.c_[genes[merge_idx[i, 0]].introns[k], genes[merge_idx[i, 1]].introns[k][:, j]]
        # merge splice graphs
        genes[merge_idx[i, 0]].splicegraph.vertices = sp.c_[genes[merge_idx[i, 0]].splicegraph.vertices, genes[merge_idx[i, 1]].splicegraph.vertices]
        m = genes[merge_idx[i, 0]].splicegraph.edges.shape[0]
        n = genes[merge_idx[i, 1]].splicegraph.edges.shape[0]
        genes[merge_idx[i, 0]].splicegraph.edges[m + 1 : n + m, :][m + 1 : n + m, :] = genes[merge_idx[i, 1]].splicegraph.edges
        genes[merge_idx[i, 0]].splicegraph.terminals = sp.c_[genes[merge_idx[i, 0]].splicegraph.terminals, genes[merge_idx[i, 1]].splicegraph.terminals]

        # extend start/stop
        genes[merge_idx[i, 0]].start = min(genes[merge_idx[i, 0]].start, genes[merge_idx[i, 1]].start)
        genes[merge_idx[i, 0]].stop = max(genes[merge_idx[i, 0]].stop, genes[merge_idx[i, 1]].stop)

        rm_map[merge_idx[i, 1]] = 1

        inserted['gene_merge'] += 1

    genes[rm_map == 1] = []

    if genes[i].splicegraph.vertices.shape[1] > 1:
        genes[i].splicegraph.uniquify()

    for i in range(len(genes)):
        assert(sp.all(genes[i].splicegraph.vertices[0, :] <= genes[i].splicegraph.vertices[1, :]))

    for ix in range(len(genes)):
        exon_order = sp.argsort(genes[ix].splicegraph.vertices[0, :])
        genes[ix].splicegraph.vertices = genes[ix].splicegraph.vertices[:, exon_order]
        genes[ix].splicegraph.edges = genes[ix].splicegraph.edges[exon_order, :][:, exon_order]
        genes[ix].splicegraph.terminals = genes[ix].splicegraph.terminals[:, exon_order]
    
    return (inserted, genes)
def insert_cassette_exons(genes, CFG):
    # written by Andre Kahles, Mpi Tuebingen, Germany, 2012
    # [genes, inserted] = insert_cassette_exons(genes, CFG)

    inserted = 0

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
                            inserted += 1 
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

    return (genes, inserted)
