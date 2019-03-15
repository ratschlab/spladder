import scipy as sp
import scipy.linalg
import scipy.sparse.linalg
import pdb

if __package__ is None:
    __package__ = 'modules'

from .utils import *
from .init import *
from .reads import *

def remove_short_exons(genes, options):

    short_exon_removed = 0
    short_exon_skipped = 0

    if options.logfile == '-':
        fd_log = sys.stdout
    else:
        fd_log = open(options.logfile, 'w')

    rm_idx = []
    for i in range(len(genes)):
        if options.verbose and (i+1) % 1000 == 0:
            print('%i' % (i+1), file=fd_log)

        ### remove genes with no exons
        if genes[i].splicegraph.vertices.shape[0] == 0:
            rm_idx.append(i)
            continue
      
        ### extend terminal exons to terminal_short_extend if they are shorter than terminal_short_len
        if genes[i].splicegraph.vertices[1, 0] - genes[i].splicegraph.vertices[0, 0] < options.remove_exons['terminal_short_len']:
            genes[i].splicegraph.vertices[0, 0] = genes[i].splicegraph.vertices[1, 0] - options.remove_exons['terminal_short_extend']
            genes[i].start = min(genes[i].start, genes[i].splicegraph.vertices[0, 0])
      
        if genes[i].splicegraph.vertices[1, -1] - genes[i].splicegraph.vertices[0, -1] < options.remove_exons['terminal_short_len']:
            genes[i].splicegraph.vertices[1, -1] = genes[i].splicegraph.vertices[1, -1] + options.remove_exons['terminal_short_extend']
            genes[i].stop = max(genes[i].stop, genes[i].splicegraph.vertices[1, -1]) 
      
        ### check for very short exons and insert an edge that allows skipping them
        exons_remove_idx = []
        j = 1
        while (j <= genes[i].splicegraph.vertices.shape[1] - 1):
            if (genes[i].splicegraph.vertices[1, j] - genes[i].splicegraph.vertices[0, j]) < options.remove_exons['min_exon_len']:
                foundp = False
                for jp in range(j + 1, genes[i].splicegraph.vertices.shape[1]):
                    if (genes[i].splicegraph.vertices[1, jp] - genes[i].splicegraph.vertices[0, jp] >= options.remove_exons['min_exon_len_remove']) and (genes[i].splicegraph.edges[j, jp] == 1):
                        foundp = True
                        break 
                foundn = False
                for jn in range(j - 2, -1, -1):
                    if (genes[i].splicegraph.vertices[1, jn] - genes[i].splicegraph.vertices[0, jn] >= options.remove_exons['min_exon_len_remove']) and (genes[i].splicegraph.edges[jn, j] == 1):
                        foundn = True
                        break
                if foundp and foundn:
                    genes[i].splicegraph.edges[jn, jp] = 1 
                    genes[i].splicegraph.edges[jp, jn] = 1 
      
                    if genes[i].splicegraph.vertices[1, j] - genes[i].splicegraph.vertices[0, j] < options.remove_exons['min_exon_len_remove']:
                        short_exon_removed += 1
                        exons_remove_idx.append(j)
                    else:
                        short_exon_skipped += 1
            j += 1
      
        keep_idx = sp.where(~sp.in1d(sp.array(sp.arange(genes[i].splicegraph.vertices.shape[1])), exons_remove_idx))[0]
        genes[i].splicegraph.subset(keep_idx)
    
    keep_idx = sp.where(~sp.in1d(sp.arange(len(genes[i])), rm_idx))[0]
    genes = genes[keep_idx]

    if options.verbose:
        print('short_exon_removed: %i' % short_exon_removed, file=fd_log)
        print('short_exon_skipped: %i' % short_exon_skipped, file=fd_log)

    if fd_log != sys.stdout:
        fd_log.close()

    return genes


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
            print('%d' % gene_idx)
      
        ### no vertices in the splice graph
        if vertices.shape[0] == 0:
          continue
      
        ### find all the intron locations
        exon_order = sp.argsort(vertices[0, :])
        vertices = vertices[:, exon_order]
        edges = edges[exon_order, :][:, exon_order]
        intron_loc = []
        for ix1 in range(vertices.shape[1] - 1):
            for ix2 in range(1, vertices.shape[1]):
                if edges[ix1, ix2] == 1:
                    intron_loc.extend(sp.arange(vertices[1, ix1], vertices[0, ix2]))
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
                             new_index = sp.r_[sp.arange(test_exon_idx), sp.arange(test_exon_idx + 1, vertices.shape[1])]
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
                            if ((vertices[1, exon_idx] >= vertices[0, test_exon_idx]) and (vertices[0, exon_idx] <= vertices[0, test_exon_idx])) or \
                               ((vertices[1, test_exon_idx] >= vertices[0, exon_idx]) and (vertices[0, test_exon_idx] <= vertices[0, exon_idx])) and \
                               (sp.sum(sp.in1d(sp.arange(min(vertices[0, exon_idx], vertices[0, test_exon_idx]), max(vertices[1, exon_idx], vertices[1, test_exon_idx])), intron_loc)) == 0):
      
                                ### merge exons if they overlap and they do not span any intronic position
                                vertices[0, exon_idx] = min(vertices[0, exon_idx], vertices[0, test_exon_idx])
                                vertices[1, exon_idx] = max(vertices[1, exon_idx], vertices[1, test_exon_idx])
                                new_index = sp.r_[sp.arange(test_exon_idx), sp.arange(test_exon_idx + 1, vertices.shape[1])]
                              
                                vertices = vertices[:, new_index]
                                edges = edges[new_index, :][:, new_index] # no need to combine any adges, as both exons have a degree of 0
                              
                                reduce_now = True
                                changed = True
                      
                        ### 0101
                        ### outgoing edges in exon and test_exo, no incoming edges
                        elif not cur_edge_left and cur_edge_right and  not test_edge_left and test_edge_right:
                            ###   ----- exon -----<
                            ###   --- test_exon --<
                            if (vertices[1, exon_idx] == vertices[1, test_exon_idx]) and \
                               (sp.sum(sp.in1d(sp.arange(min(vertices[0, exon_idx], vertices[0, test_exon_idx]), vertices[1, exon_idx]), intron_loc)) == 0):
                                
                                ### merge exons if they share the same right boundary and do not span intronic positions
                                vertices[0, exon_idx] = min(vertices[0, exon_idx], vertices[0, test_exon_idx])
                                new_index = sp.r_[sp.arange(test_exon_idx), sp.arange(test_exon_idx + 1, vertices.shape[1])]
                              
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
                            if (vertices[0, exon_idx] == vertices[0, test_exon_idx]) and \
                               (sp.sum(sp.in1d(sp.arange(vertices[0, exon_idx], max(vertices[1, exon_idx], vertices[1, test_exon_idx])), intron_loc)) == 0):
                                
                                ### merge exons if they share the same left boundary and do not span intronic positions
                                vertices[1, exon_idx] = max(vertices[1, exon_idx], vertices[1, test_exon_idx])
                                new_index = sp.r_[sp.arange(test_exon_idx), sp.arange(test_exon_idx + 1, vertices.shape[1])]
                              
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
                                new_index = sp.r_[sp.arange(test_exon_idx), sp.arange(test_exon_idx + 1, vertices.shape[1])]
                              
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


def filter_by_edgecount(genes, options):

    ### filter splicegraphs by support count over samples
    keep_genes = []
    for i in range(len(genes)):
        #u_exons = unique_rows(sp.vstack(genes[i].exons))
        #(tmp, tmp, k_idx) = intersect_rows(u_exons, genes[i].splicegraph.vertices.T, index=True)
        genes[i].from_sparse()
        k_idx = sp.where(genes[i].splicegraph.edges.sum(axis = 1) == 0)[0]
        genes[i].splicegraph.edges = (genes[i].edge_count >= options.sg_min_edge_count).astype('int')
        ### remove all exons that have no incoming or outgoing edges (caused by validation, keep single exon transcripts that occured before)
        k_idx2 = sp.where(genes[i].splicegraph.edges.sum(axis = 1) > 0)[0]
        #rm_idx = sp.where(~sp.in1d(k_idx2, k_idx))[0]
        #keep_idx = sp.where(~sp.in1d(sp.arange(genes[i].splicegraph.edges.shape[0]), rm_idx))[0]
        keep_idx = sp.union1d(k_idx, k_idx2).astype('int')
        if keep_idx.shape[0] > 0:
            genes[i].splicegraph.subset(keep_idx)
            keep_genes.append(i)
        genes[i].to_sparse()

    return genes[keep_genes]


def insert_intron_retentions(genes, options):
    # written by Georg Zeller & Gunnar Raetsch, Mpi Tuebingen, Germany, 2009, Andre Kahles, MSKCC, 2013

    inserted = 0
    strands = ['+', '-']

    ### form all possible combinations of contigs and strands --> regions
    (regions, options) = init_regions(options.bam_fnames, options.confidence, options, sparse_bam=options.sparse_bam)

    ### ignore contigs not present in bam files 
    # TODO
    #keepidx = sp.where(sp.in1d(sp.array([options.chrm_lookup[x.chr] for x in genes[chunk_idx]]), sp.array([x.chr_num for x in regions])))[0]
    #genes = genes[keepidx]

    c = 0 
    num_introns_added = 0
    num_introns = 0

    contigs = sp.array([x.chr for x in genes], dtype='str')
    gene_strands = sp.array([x.strand for x in genes])
    for contig in sp.unique(contigs):
        bam_cache = dict()
        for si, s in enumerate(strands):
            cidx = sp.where((contigs == contig) & (gene_strands == s))[0]

            for i in cidx:

                if options.verbose and (c+1) % 100 == 0:
                    print('\r %i(%i) genes done (found %i new retentions in %i tested eligible introns, %2.1f%%)' % (c+1, genes.shape[0], num_introns_added, num_introns, 100 * num_introns_added / float(max(1, num_introns))), file=sys.stdout)

                gg = genes[i]
                assert(gg.strand == s)
                assert(gg.chr == contig)

                if options.sparse_bam:
                    if isinstance(options.bam_fnames, str):
                        [tracks] = add_reads_from_sparse_bam(gg, options.bam_fnames, contig, options.confidence, types=['exon_track'], filter=options.read_filter, cache=bam_cache, unstranded=True)
                    else:
                        tracks = None
                        for fname in options.bam_fnames:
                            [tmp_] = add_reads_from_sparse_bam(gg, fname, contig, options.confidence, types=['exon_track'], filter=options.read_filter, cache=bam_cache, unstranded=True)
                            if tracks is None:
                                tracks = tmp_
                            else:
                                tracks += tmp_
                    tracks = sp.asarray(tracks)
                else:
                    tracks = add_reads_from_bam(sp.array([gg], dtype='object'), options.bam_fnames, ['exon_track'], options.read_filter, options.var_aware, options.primary_only, options.ignore_mismatches)

                exon_coverage = sp.zeros((gg.splicegraph.vertices.shape[1],), dtype='float')
                for k in range(gg.splicegraph.vertices.shape[1]):
                    idx = sp.arange(gg.splicegraph.vertices[0, k], gg.splicegraph.vertices[1, k]) - gg.start
                    exon_coverage[k] = sp.median(sp.sum(tracks[:, idx], axis=0).astype('float')) # median coverage for exon k

                ### check for all vertex-pairs, if respective intron can be retained
                new_retention = sp.zeros(gg.splicegraph.edges.shape, dtype='int') 
                #for k in range(gg.splicegraph.edges.shape[0]):
                #    for l in range(k + 1, gg.splicegraph.edges.shape[0]):
                #        if gg.splicegraph.edges[k, l] == 1:
                for k,l in sp.array(sp.where(sp.triu(gg.splicegraph.edges))).T:
                    ### ignore introns that include at least one complete exon
                    if sp.sum((gg.splicegraph.vertices[0, :] > gg.splicegraph.vertices[1, k]) & (gg.splicegraph.vertices[1, :] < gg.splicegraph.vertices[0, l])) > 0:
                        continue

                    num_introns += 1
                    idx = sp.arange(gg.splicegraph.vertices[1, k], gg.splicegraph.vertices[0, l]) - gg.start
                    icov = sp.sum(tracks[:, idx], axis=0) 
                    if sp.median(icov) > options.intron_retention['min_retention_cov'] and \
                        sp.mean(icov > (0.5 * sp.mean(icov))) > options.intron_retention['min_retention_region'] and  \
                        max(exon_coverage[k], exon_coverage[l]) / (1e-6 + min(exon_coverage[k], exon_coverage[l])) <= options.intron_retention['min_retention_max_exon_fold_diff'] and \
                        sp.mean(icov) >= options.intron_retention['min_retention_rel_cov'] * (exon_coverage[k] + exon_coverage[l]) / 2.0 and \
                        sp.mean(icov) <= options.intron_retention['max_retention_rel_cov'] * (exon_coverage[k] + exon_coverage[l]) / 2.0:

                        new_retention[k, l] = 1
                        inserted += 1
                    #	fprintf(log_fd, '%s\tintron_retention\t%c\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%2.1f\n', gg.chr, gg.strand, gg.splicegraph{1}(1,k), gg.splicegraph{1}(2,k), gg.splicegraph{1}(1,l), gg.splicegraph{1}(2,l), ...
                    #			floor(median(icov(1,:)+icov(2,:))), floor(gg.exon_coverage(k)), floor(gg.exon_coverage(l)), 100*mean(icov(1,:)+icov(2,:)>0)) ;
                any_added = False
                if False:
                    if sp.sum(new_retention.ravel()) > 0:
                        new_retention = scipy.sparse.linalg.expm(new_retention)
                        new_retention[new_retention == 0] = 2
                        sp.fill_diagonal(new_retention, 2)
                        while True:
                            any_added = False
                            k,l = sp.unravel_index(new_retention.argmin(), new_retention.shape)
                            if new_retention[k, l] == 2:
                                break
                            if new_retention[k, l] > 0:
                                gg.splicegraph.add_intron_retention(k, l)
                                new_retention = sp.c_[new_retention, sp.ones((new_retention.shape[0], 1), dtype='int') * 2]
                                new_retention = sp.r_[new_retention, sp.ones((1, new_retention.shape[1]), dtype='int') * 2]
                                ### unset all inbetween retentions
                                for u in range(k, l + 1):
                                    for v in range(u + 1, l + 1):
                                        new_retention[u, v] = 2
                                #new_retention[k, l] = 0
                                any_added = True
                                num_introns_added += 1
                                #fprintf(log_fd, '%s\tintron_retention\t%i\t%i\t%i\t%i\t%i\t%2.1f\n', gg.chr, gg.splicegraph{1}(2,k), gg.splicegraph{1}(1,l), floor(median(icov(1,:)+icov(2,:))), gg.exon_coverage(k), gg.exon_coverage(l), 100*mean(icov(1,:)+icov(2,:)>0)) ;
                            exon_order = sp.argsort(gg.splicegraph.vertices[0, :])
                            gg.splicegraph.reorder(exon_order)
                            new_retention = new_retention[exon_order, :][:, exon_order]
                            if not any_added:
                                break
                else:
                    if sp.sum(new_retention.ravel()) > 0:
                        new_retention = scipy.sparse.linalg.expm(new_retention)
                        while True:
                            any_added = False
                            for k in range(new_retention.shape[1]):
                                for l in range(k + 1, new_retention.shape[1]):
                                    if new_retention[k, l] > 0:
                                        gg.splicegraph.add_intron_retention(k, l)
                                        new_retention = sp.c_[new_retention, sp.zeros((new_retention.shape[0], 1), dtype='int')]
                                        new_retention = sp.r_[new_retention, sp.zeros((1, new_retention.shape[1]), dtype='int')]
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
                genes[i] = gg
                c += 1
    return (genes, inserted)


def insert_intron_edges(genes, options):

    if not hasattr(options, 'debug'):
        options.debug = False

    print_intermediates = False

    strands = ['+', '-']
    P = []
    both_missing = [0, 0]
    one_missing = [0, 0]
    multi = 0
    next = 0
    prev = 0

    exon_vicinity_cnt1 = [0, 0] 
    exon_vicinity_cnt2 = [0, 0] 
    merge_idx = sp.zeros((0,))
    intron_tol = 0 

    inserted = dict()
    inserted['intron_in_exon'] = 0 
    inserted['alt_53_prime'] = 0 
    inserted['exon_skip'] = 0 
    inserted['gene_merge'] = 0 
    inserted['new_terminal_exon'] = 0 

    num_unused_introns = sp.zeros((genes.shape[0],), dtype='int')

    last_chr = ''

    if options.logfile == '-':
        fd_log = sys.stdout
    else:
        fd_log = open(options.logfile, 'w')

    for i in range(genes.shape[0]):
    
        if options.verbose and (i+1) % 1000 == 0:
            print('%i of %i genes' % (i+1, genes.shape[0]), file=fd_log)

        if options.sparse_bam and genes[i].chr != last_chr:
            bam_cache = dict()
        last_chr = genes[i].chr

        s = strands.index(genes[i].strand)

        if genes[i].introns is None or genes[i].introns[s].shape[0] == 0:
            continue

        unused_introns = []
        if options.debug:
            print('processing gene %i; with %i introns' % (i, len(genes[i].introns[s])), file=fd_log)
            ### TODO timing

        for j in range(genes[i].introns[s].shape[0]):
            intron_used = False

            if (j > 0 and genes[i].splicegraph.vertices.shape[1] > 1):
                genes[i].splicegraph.uniquify()

            ### find exons within same gene whose end coincides with intron start
            idx1 = sp.where(sp.absolute(genes[i].splicegraph.vertices[1, :] - genes[i].introns[s][j, 0]) <= intron_tol)[0]
            ### find exons within same gene whose start coincides with intron end
            idx2 = sp.where(sp.absolute(genes[i].splicegraph.vertices[0, :] - genes[i].introns[s][j, 1]) <= intron_tol)[0]

            ### intron boundaries do not coincide with any exon boundaries
            if idx1.shape[0] == 0 and idx2.shape[0] == 0:
                both_missing[s] += 1

                if options.intron_edges['insert_intron_retention']:
                    ### find all exons that completely include added introns 
                    idx1__ = sp.where((genes[i].introns[s][j, 0] > genes[i].splicegraph.vertices[0, :]) & (genes[i].introns[s][j, 1] < genes[i].splicegraph.vertices[1, :]))[0]
                    for idx1_ in idx1__:
                        genes[i].splicegraph.vertices = sp.c_[genes[i].splicegraph.vertices, genes[i].splicegraph.vertices[:, idx1_]]
                        genes[i].splicegraph.vertices[1, -1] = genes[i].introns[s][j, 0]
                                
                        genes[i].splicegraph.vertices = sp.c_[genes[i].splicegraph.vertices, genes[i].splicegraph.vertices[:, idx1_]]
                        genes[i].splicegraph.vertices[0, -1] = genes[i].introns[s][j, 1]
                                
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
                                
                        genes[i].splicegraph.terminals = sp.c_[genes[i].splicegraph.terminals, genes[i].splicegraph.terminals[:, idx1_]]
                        genes[i].splicegraph.terminals[1, -1] = 0 # cannot be an end
                        genes[i].splicegraph.terminals = sp.c_[genes[i].splicegraph.terminals, genes[i].splicegraph.terminals[:, idx1_]]
                        genes[i].splicegraph.terminals[0, -1] = 0 # cannot be a start
                                
                        inserted['intron_in_exon'] += 1
                        assert(sp.all(genes[i].splicegraph.vertices[0, :] <= genes[i].splicegraph.vertices[1, :]))

                        if options.debug:
                            print('%s\tintron_retention_exon\t%c\t%i\t%i\t%i\t%i\n' % (genes[i].chr, genes[i].strand, genes[i].splicegraph.vertices[0, -2],
                                                                                                         genes[i].splicegraph.vertices[1, -2], genes[i].splicegraph.vertices[0, -1], 
                                                                                                         genes[i].splicegraph.vertices[1, -1]), file=fd_log)
                        intron_used = True

                if not intron_used:
                    unused_introns.append(j)
                continue # with next intron

            # did not find exons in same gene sharing boundaries with intron start
            # find first end in previous gene on same strand
            if idx1.shape[0] == 0 and i > 0 and len(genes) > 1  and (genes[i - 1].chr == genes[i].chr) and (genes[i - 1].strand == genes[i].strand): 
                ### find all exon ends in previuos gene that coincide with intron start j
                idx1_ = sp.where(sp.absolute(genes[i-1].splicegraph.vertices[1, :] - genes[i].introns[s][j, 0]) <= intron_tol)[0]
                if idx1_.shape[0] > 0:
                    prev += 1
                    # mark the two genes for merging
                    if options.intron_edges['gene_merges']:
                        merge_idx = sp.c_[merge_idx, sp.array([i-1, i])]
                        intron_used = True
                    if not intron_used:
                        unused_introns.append(j)
                    continue # with next intron

            # did not find exons in same gene sharing boundaries with intron end
            # find second end in next gene on same strand
            if idx2.shape[0] == 0 and len(genes) > 1 and i + 1 < len(genes) and genes[i + 1].chr == genes[i].chr and genes[i+1].strand == genes[i].strand:
                ### find all exon starts in following gene that coincide with intron end j
                idx2_ = sp.where(sp.absolute(genes[i+1].splicegraph.vertices[0, :] - genes[i].introns[s][j, 1]) <= intron_tol)[0]
                if idx2_.shape[0] > 0:
                    next += 1
                    # mark the two genes for merging
                    if options.intron_edges['gene_merges']: 
                        merge_idx = sp.c_[merge_idx, sp.array([i, i+1])]
                        intron_used = True
                    if not intron_used:
                        unused_introns.append(j)
                    continue # with next intron

            # did not find exons in same gene sharing boundaries with intron start
            # check whether the intron starts in the vicinity of an exon
            if idx1.shape[0] == 0: 
                ### find all exons that overlap intron-start j +/- options.intron_edges.vicinity_region
                idx1__ = sp.where((genes[i].splicegraph.vertices[0, :] - options.intron_edges['vicinity_region'] <= genes[i].introns[s][j, 0]) & 
                                  (genes[i].splicegraph.vertices[1, :] + options.intron_edges['vicinity_region'] > genes[i].introns[s][j, 0]))[0]

                ### check, if we can find an exon before the current intron and there is continuous coverage between intron end and exon
                if idx1__.shape[0] == 0:
                    #idx1__ = sp.argmax(genes[i].splicegraph.vertices[0, :] >= genes[i].introns[s][j, 1])
                    idx1__ = sp.where(genes[i].splicegraph.vertices[1, :] <= genes[i].introns[s][j, 0])[0]
                    if len(idx1__.shape) > 0 and idx1__.shape[0] > 0:
                        idx1__ = idx1__.max()[sp.newaxis]
                        gg = genes[i]
                        gg.strand = strands[s]
                        gg.strands = strands[s]
                        gg.start = genes[i].splicegraph.vertices[1, idx1__][0] ### stop of previous exon
                        gg.stop = genes[i].introns[s][j, 0]  ### end of presumable exon

                        if options.sparse_bam:
                            if isinstance(options.bam_fnames, str):
                                [tracks] = add_reads_from_sparse_bam(gg, options.bam_fnames, gg.chr, options.confidence, types=['exon_track'], filter=options.read_filter, cache=bam_cache, unstranded=True)
                            else:
                                tracks = None
                                for fname in options.bam_fnames:
                                    [tmp_] = add_reads_from_sparse_bam(gg, fname, gg.chr, options.confidence, types=['exon_track'], filter=options.read_filter, cache=bam_cache, unstranded=True)
                                    if tracks is None:
                                        tracks = tmp_
                                    else:
                                        tracks += tmp_
                            tracks = sp.asarray(tracks)
                        else:
                            tracks = add_reads_from_bam(sp.array([gg], dtype='object'), options.bam_fnames, ['exon_track'], options.read_filter, options.var_aware, options.primary_only, options.ignore_mismatches)

                        ### TODO: make the following a configurable
                        if sp.mean(sp.sum(tracks, axis=0) > 10) < 0.9:
                            idx1__ = sp.array([])
                
                # only take the case closest to an existing splice site
                if len(idx1__.shape) > 0 and idx1__.shape[0] > 0:
                    diff1 = sp.absolute(genes[i].splicegraph.vertices[0, idx1__] - genes[i].introns[s][j, 0])
                    diff2 = sp.absolute(genes[i].splicegraph.vertices[1, idx1__] - genes[i].introns[s][j, 0] - 1)
                    diff = sp.minimum(diff1, diff2)
                    if diff.shape[0] > 0:
                        idx1__ = sp.array([idx1__[sp.argmin(diff)]])
                        for idx1_ in idx1__:
                            if genes[i].introns[s][j, 0] - genes[i].splicegraph.vertices[0, idx1_] >= options.intron_edges['min_exon_len']:
                                exon_vicinity_cnt1[s] += 1
                                genes[i].splicegraph.vertices = sp.c_[genes[i].splicegraph.vertices, genes[i].splicegraph.vertices[:, idx1_]]
                                genes[i].splicegraph.vertices[1, -1] = genes[i].introns[s][j, 0]  # set exon end to intron start (half open)
                                genes[i].splicegraph.new_edge()
                                genes[i].splicegraph.edges[:, -1] = genes[i].splicegraph.edges[:, idx1_]
                                genes[i].splicegraph.edges[-1, :] = genes[i].splicegraph.edges[idx1_, :]
                                genes[i].splicegraph.terminals = sp.c_[genes[i].splicegraph.terminals, genes[i].splicegraph.terminals[:, idx1_]]  # copy from original exon
                                genes[i].splicegraph.terminals[1, -1] = 0 # cannot be an end
                                            
                                assert(sp.all(genes[i].splicegraph.vertices[0, :] <= genes[i].splicegraph.vertices[1, :]))

                                # check exons whose start coincides with intron end
                                genes[i].splicegraph.add_intron(sp.array([genes[i].splicegraph.edges.shape[0] - 1]), 0, idx2, 1)
                                            
                                inserted['alt_53_prime'] += 1
                                if options.debug:
                                    for idx2_ in idx2:
                                        print('%s\talternative_53_prime1\t%c\t%i\t%i\t%i\n' % (genes[i].chr, genes[i].strand, genes[i].splicegraph,vertices[1, idx1_], 
                                                                                                                 genes[i].splicegraph.vertices[1, -1], genes[i].splicegraph.vertices[0, idx2_]), file=fd_log)
                                intron_used = True

                ### if no proximal exon was found, insert new terminal exon, if wished
                if  not intron_used and options.intron_edges['append_new_terminal_exons']:
                    inserted['new_terminal_exon'] += 1

                    iregion = sp.array([[max(0, genes[i].introns[s][j, 0] - options.intron_edges['append_new_terminal_exons_len'])], [genes[i].introns[s][j, 0]]])
                    idx_iregion = sp.where((genes[i].introns[s][:, 1] >= iregion[0]) & (genes[i].introns[s][:, 1] < iregion[1]))[0]
                    if idx_iregion.shape[0] > 0:
                        if not idx_iregion.shape[0] == 1:
                            idx_iregion = idx_iregion[sp.argmax(genes[i].introns[s][idx_iregion, 1])]
                        iregion[0] = genes[i].introns[s][idx_iregion, 1]
                        assert(iregion[0] < iregion[1])

                    genes[i].splicegraph.vertices = sp.c_[genes[i].splicegraph.vertices, iregion]
                    genes[i].splicegraph.new_edge()
                    genes[i].splicegraph.terminals = sp.c_[genes[i].splicegraph.terminals, sp.array([1, 0])] # can be a start, but cannot be an end

                    ### correct terminal exon starts
                    for tmp_idx in idx2:
                        if genes[i].splicegraph.terminals[0, tmp_idx] == 1 and genes[i].introns[s][j, 1] < genes[i].splicegraph.vertices[1, tmp_idx]:
                            genes[i].splicegraph.vertices[0, tmp_idx] = genes[i].introns[s][j, 1]
                    assert(sp.all(genes[i].splicegraph.vertices[1, :] >= genes[i].splicegraph.vertices[0, :]))

                    genes[i].splicegraph.add_intron(idx2, 1, sp.array([genes[i].splicegraph.edges.shape[0] - 1]), 0)
                    genes[i].splicegraph.uniquify()
                    intron_used = True

                if intron_used:
                    continue

            # did not find exons in same gene sharing boundaries with intron end
            # check whether the intron ends in the vicinity of an exon
            if idx2.shape[0] == 0: 
                idx2__ = sp.where((genes[i].splicegraph.vertices[0, :] - options.intron_edges['vicinity_region'] < genes[i].introns[s][j, 1]) &
                                  (genes[i].splicegraph.vertices[1, :] + options.intron_edges['vicinity_region'] >= genes[i].introns[s][j, 1]))[0]

                ### check, if we can find an exon after the current intron and there is continuous coverage between intron end and exon
                if idx2__.shape[0] == 0:
                    #idx2__ = sp.argmax(genes[i].splicegraph.vertices[0, :] > genes[i].introns[s][j, 1])
                    idx2__ = sp.where(genes[i].splicegraph.vertices[0, :] > genes[i].introns[s][j, 1])[0]
                    if len(idx2__.shape) > 0 and idx2__.shape[0] > 0:
                        idx2__ = idx2__.min()[sp.newaxis]
                        gg = genes[i]
                        gg.strand = strands[s]
                        gg.strands = strands[s]
                        gg.start = genes[i].introns[s][j, 1]  ### start of presumable exon
                        gg.stop = genes[i].splicegraph.vertices[1, idx2__][0]  ### stop of next exon

                        if options.sparse_bam:
                            if isinstance(options.bam_fnames, str):
                                [tracks] = add_reads_from_sparse_bam(gg, options.bam_fnames, gg.chr, options.confidence, types=['exon_track'], filter=options.read_filter, cache=bam_cache, unstranded=True)
                            else:
                                tracks = None
                                for fname in options.bam_fnames:
                                    [tmp_] = add_reads_from_sparse_bam(gg, fname, gg.chr, options.confidence, types=['exon_track'], filter=options.read_filter, cache=bam_cache, unstranded=True)
                                    if tracks is None:
                                        tracks = tmp_
                                    else:
                                        tracks += tmp_
                            tracks = sp.asarray(tracks)
                        else:
                            tracks = add_reads_from_bam(sp.array([gg], dtype='object'), options.bam_fnames, ['exon_track'], options.read_filter, options.var_aware, options.primary_only, options.ignore_mismatches)

                        ### TODO: make configurable
                        if sp.mean(sp.sum(tracks, axis=0) > 10) < 0.9:
                            idx2__ = sp.array([])

                # only take the case closest to an existing splice site
                if len(idx2__.shape) > 0 and idx2__.shape[0] > 0:
                    diff1 = sp.absolute(genes[i].splicegraph.vertices[0, idx2__] - genes[i].introns[s][j, 1] + 1)
                    diff2 = sp.absolute(genes[i].splicegraph.vertices[1, idx2__] - genes[i].introns[s][j, 1])
                    diff = sp.minimum(diff1, diff2)
                    if diff.shape[0] > 0:
                        idx2__ = sp.array([idx2__[sp.argmin(diff)]])
                        for idx2_ in idx2__:
                            if genes[i].splicegraph.vertices[1, idx2_] - genes[i].introns[s][j, 1] >= options.intron_edges['min_exon_len']:
                                exon_vicinity_cnt2[s] += 1
                                genes[i].splicegraph.vertices = sp.c_[genes[i].splicegraph.vertices, genes[i].splicegraph.vertices[:, idx2_]]
                                genes[i].splicegraph.vertices[0, -1] = genes[i].introns[s][j, 1]
                                genes[i].splicegraph.new_edge()
                                genes[i].splicegraph.edges[:, -1] = genes[i].splicegraph.edges[:, idx2_]
                                genes[i].splicegraph.edges[-1, :] = genes[i].splicegraph.edges[idx2_, :]
                                genes[i].splicegraph.terminals = sp.c_[genes[i].splicegraph.terminals, genes[i].splicegraph.terminals[:, idx2_]]   # copy from original exon
                                genes[i].splicegraph.terminals[0, -1] = 0   # cannot be a start
                                    
                                assert(sp.all(genes[i].splicegraph.vertices[0, :] <= genes[i].splicegraph.vertices[1, :]))

                                genes[i].splicegraph.add_intron(idx1, 1, sp.array([genes[i].splicegraph.edges.shape[0] - 1]), 0)
                                
                                inserted['alt_53_prime'] += 1

                                if options.debug:
                                    for idx1_ in idx1:
                                        print('%s\talternative_53_prime2\t%c\t%i\t%i\t%i\n' % (genes[i].chr, genes[i].strand, genes[i].splicegraph.vertices[1, idx1_], 
                                                                                                                 genes[i].splicegraph.vertices[0, -1], genes[i].splicegraph.vertices[0, idx2_]), file=fd_log)
                                intron_used = True

                ### if no proximal exon was found, insert new terminal exon, if wished
                if not intron_used and options.intron_edges['append_new_terminal_exons']:

                    ### define range of new exon
                    iregion = sp.array([[genes[i].introns[s][j, 1]], [genes[i].introns[s][j, 1] + options.intron_edges['append_new_terminal_exons_len']]])
                    ### find introns starting within new exon
                    idx_iregion = sp.where((genes[i].introns[s][:, 0] > iregion[0]) & (genes[i].introns[s][:, 0] <= iregion[1]))[0]

                    if idx_iregion.shape[0] > 0:
                        if not idx_iregion.shape[0] == 1: 
                            idx_iregion = idx_iregion[sp.argmin(genes[i].introns[s][idx_iregion, 0])]
                        ### let new exon end at position before next intron starts
                        iregion[1] = genes[i].introns[s][idx_iregion, 0]
                        assert(iregion[0] < iregion[1])

                    inserted['new_terminal_exon'] += 1
                    genes[i].splicegraph.vertices = sp.c_[genes[i].splicegraph.vertices, iregion]
                    genes[i].splicegraph.new_edge()
                    genes[i].splicegraph.terminals = sp.c_[genes[i].splicegraph.terminals, sp.array([0, 1])]  # cannot be a start but can be an end

                    ### adapt terminal exon ends if the new intron starts within them
                    for tmp_idx in idx1:
                        if genes[i].splicegraph.terminals[1, tmp_idx] == 1 and genes[i].introns[s][j, 0] > genes[i].splicegraph.vertices[0, tmp_idx]:
                            genes[i].splicegraph.vertices[1, tmp_idx] = genes[i].introns[s][j, 0]
                    assert(sp.all(genes[i].splicegraph.vertices[1, :] >= genes[i].splicegraph.vertices[0, :]))
                            
                    genes[i].splicegraph.add_intron(idx1, 1, sp.array([genes[i].splicegraph.edges.shape[0] - 1]), 0)
                    intron_used = True

                if intron_used:
                    continue
            
            if idx1.shape[0] == 0 or idx2.shape[0] == 0:
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
            ### insert exon skips     TODO this is not always a skipped exon, can also reconnect alternative ends, etc.
            for idx1_ in idx1:
                for idx2_ in idx2:
                    if genes[i].splicegraph.edges[idx1_, idx2_] == 0:
                        inserted['exon_skip'] += 1
                                    
            genes[i].splicegraph.add_intron(idx1, 1, idx2, 1)

        unused_introns = sp.array(unused_introns, dtype='int')
        idx_unused = sp.where((genes[i].introns[s][unused_introns, 1] >= genes[i].start) & (genes[i].introns[s][unused_introns, 0] <= genes[i].stop))[0]
        unused_introns = unused_introns[idx_unused]
        if options.debug and unused_introns.shape[0] > 0:
            print('Warning: unused introns: %s' % str(unused_introns))
        num_unused_introns[i] += unused_introns.shape[0]

    if print_intermediates:
        print('one missing: %s' % str(one_missing), file=fd_log)
        print('multi: %s' % str(multi), file=fd_log)
        print('num_unused_introns: $i' % sp.sum(num_unused_introns), file=fd_log)

    if fd_log != sys.stdout:
        fd_log.close()

    merge_idx = unique_rows(merge_idx)
    rm_map = sp.zeros_like(genes, dtype='int')

    for i in range(merge_idx.shape[0]):
        
        while rm_map[merge_idx[i, 0]] == 1: 
            merge_idx[i, 0] -= 1

        # merge transcripts
        for j in range(genes[merge_idx[i, 1]].exons.shape[0]):
            genes[merge_idx[i, 0]].transcripts.append(genes[merge_idx[i, 1]].transcripts[j])
            genes[merge_idx[i, 0]].exons.append(genes[merge_idx[i, 1]].exons[j])
        # merge intron lists
        # TODO !!! check this
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
        genes[ix].splicegraph.sort()
    
    return (genes, inserted)


def insert_cassette_exons(genes, options):
    # written by Andre Kahles, Mpi Tuebingen, Germany, 2012

    inserted = 0

    strands = ['+', '-']

    ### form all possible combinations of contigs and strands --> regions
    (regions, options) = init_regions(options.bam_fnames, options.confidence, options, sparse_bam=options.sparse_bam)

    ### ignore contigs not present in bam files 
    # TODO TODO
    #keepidx = sp.where(sp.in1d(sp.array([CFG['chrm_lookup'][x.chr] for x in genes]), sp.array([x.chr_num for x in regions])))[0]
    #genes = genes[keepidx]

    c = 0
    num_exons_added = 0
    num_exons = 0

    contigs = sp.array([x.chr for x in genes], dtype='str')
    gene_strands = sp.array([x.strand for x in genes])
    for contig in sp.unique(contigs):
        bam_cache = dict()
        for si, s in enumerate(strands):
            cidx = sp.where((contigs == contig) & (gene_strands == s))[0]

            for i in cidx:
                if options.verbose and (c+1) % 100 == 0:
                    print('\r %i(%i) genes done (found %i new cassette exons in %i tested intron pairs, %2.1f%%)' % (c+1, genes.shape[0], num_exons_added, num_exons, 100*num_exons_added/float(max(1, num_exons))))

                gg = genes[i]
                assert(gg.strand == s)
                assert(gg.chr == contig)

                if options.sparse_bam:
                    if isinstance(options.bam_fnames, str):
                        [tracks] = add_reads_from_sparse_bam(gg, options.bam_fnames, contig, options.confidence, types=['exon_track'], filter=options.read_filter, cache=bam_cache, unstranded=True)
                    else:
                        tracks = None
                        for fname in options.bam_fnames:
                            [tmp_] = add_reads_from_sparse_bam(gg, fname, contig, options.confidence, types=['exon_track'], filter=options.read_filter, cache=bam_cache, unstranded=True)
                            if tracks is None:
                                tracks = tmp_
                            else:
                                tracks += tmp_
                    tracks = sp.asarray(tracks)
                else:
                    tracks = add_reads_from_bam(sp.array([gg], dtype='object'), options.bam_fnames, ['exon_track'], options.read_filter, options.var_aware, options.primary_only, options.ignore_mismatches)

                ### add introns implied by splicegraph to the list
                all_introns = gg.introns[si][:, :2]
                for k in range(gg.splicegraph.edges.shape[0]):
                    for l in range(k+1, gg.splicegraph.edges.shape[0]):
                        if gg.splicegraph.edges[k, l] == 1:
                            all_introns = sp.r_[all_introns, sp.array([[gg.splicegraph.vertices[1, k], gg.splicegraph.vertices[0, l]]])] # introns are half open
                all_introns = unique_rows(all_introns)
           
                ### use only relevant introns (inside gene boundaries)
                if all_introns.shape[0] > 0:
                    keep_idx = sp.where((all_introns[:, 1] > gg.start) & (all_introns[:, 0] < gg.stop))[0]
                    all_introns = all_introns[keep_idx, :]

                segment_starts = sp.sort(sp.unique(all_introns[:, 0]))
                segment_ends = sp.sort(sp.unique(all_introns[:, 1]))

                ### check for all intron-pairs, if exon could exist between them
                new_cassette = sp.zeros((all_introns.shape[0], all_introns.shape[0]), dtype='int') 
                for k in range(all_introns.shape[0]):
                    for l in range(k + 1, all_introns.shape[0]):
                        if all_introns[k, 1] >= all_introns[l, 0]:
                            continue
                        ### only take intron pair, if outer ends are supported by current exons
                        if (not all_introns[k, 0] in gg.splicegraph.vertices[1, :]) or (not all_introns[l, 1] in gg.splicegraph.vertices[0, :]): 
                            continue
                        curr_exon = [all_introns[k, 1], all_introns[l, 0]]
                        ### do not allow curr_exon to overlap existing exon
                        if sp.sum((gg.splicegraph.vertices[0, :] < curr_exon[1]) & (gg.splicegraph.vertices[1, :] > curr_exon[0])) > 0:
                            continue

                        num_exons += 1

                        if not ismember(curr_exon, gg.splicegraph.vertices.T, rows=True):
                            idx = sp.arange(curr_exon[0], curr_exon[1]) - gg.start
                            exon_cov = sp.sum(tracks[:, idx], axis=0)

                            pre_segment_end = sp.where(segment_ends < curr_exon[0])[0]
                            if pre_segment_end.shape[0] > 0:
                                pre_segment_cov = sp.sum(tracks[:, sp.arange(segment_ends[pre_segment_end.max()], curr_exon[0]) - gg.start], axis=0)
                            else:
                                pre_segment_cov = sp.sum(tracks[:, sp.arange(curr_exon[0] - gg.start)], axis=0)
                            min_len_pre = min(pre_segment_cov.shape[0], exon_cov.shape[0])

                            aft_segment_start = sp.where(segment_starts > curr_exon[1])[0]
                            if aft_segment_start.shape[0] > 0:
                                aft_segment_cov = sp.sum(tracks[:, sp.arange(curr_exon[1], segment_starts[aft_segment_start.min()]) - gg.start], axis=0)
                            else:
                                aft_segment_cov = sp.sum(tracks[:, (curr_exon[1] - gg.start):], axis=0)
                            min_len_aft = min(aft_segment_cov.shape[0], exon_cov.shape[0])

                            if sp.mean(exon_cov > (0.2 * sp.mean(exon_cov))) > options.cassette_exon['min_cassette_region'] and \
                               sp.median(exon_cov) > options.cassette_exon['min_cassette_cov']:
                                median_aft = max(sp.median(aft_segment_cov[:min_len_aft]), 1)
                                median_pre = max(sp.median(pre_segment_cov[-min_len_pre:]), 1)
                                if (sp.median(exon_cov[-min_len_aft:]) / median_aft) - 1 >= options.cassette_exon['min_cassette_rel_diff'] and \
                                   (sp.median(exon_cov[:min_len_pre]) / median_pre) - 1 >= options.cassette_exon['min_cassette_rel_diff']:
                                    new_cassette[k, l] = 1
                                    inserted += 1 
                any_added = False
                if any(new_cassette.ravel()):
                    curr_sg = gg.splicegraph.vertices
                    for k in range(new_cassette.shape[1]):
                        for l in range(k + 1, new_cassette.shape[1]):
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
                genes[i] = gg
                c += 1

    return (genes, inserted)


def infer_splice_graph_caller(genes):
    # genes = infer_splice_graph_caller(genes)
      
    MISALIGN = 5

    print('Number of genes:%d' % genes.shape[0])

    #####################################################################
    # Some simple merging to reduce exons
    #####################################################################
    print('Merging some genes')
    genes = reduce_splice_graph(genes)


    #####################################################################
    #
    # infer_splice_graph infers the relevant splice graph given a set of ESTs. 
    #
    # It is based on two assumptions, neither which are true in reality, but
    # are necessary for any reasonable processing.
    # 1. All splicing events are independent. That is a splice site
    #    at location A in ESTA does not depend on splice site B in ESTA, and
    #    neither does it depend on splice site C in ESTB.
    # 2. There exists a splicing event at a particular location if and only if
    #    there exists an EST that is spliced at this point.
    #
    #####################################################################

    print('performing inference')

    for gene_idx in range(genes.shape[0]):

        if gene_idx % 1000 ==0:
            print('%d ' % gene_idx)
      
        #####################################################################
        # For each gene, we consider each pair of exons in the splicegraph.
        # After sorting the exons by the location of the left side, we make
        # a local copy of the splicegraph.
        #
        # Note that 'left' means the 5 prime end of a positive strand, but
        # the 3 prime end of a negative strand.
        #
        # The two exons considered are indexed by exon_idx, and test_exon_idx.
        # The general idea of the algorithm is as follows:
        # For each pair:
        #   check which case it is in (more on this later)
        #   if it is possible to merge:
        #     form a new exon
        #   end merge.
        #   if new exon is not in current set
        #     add to splicegraph
        #   end new exon insert
        # end for
        # 
        # In words:
        # Each time two exons are considered, and they can be merged to form
        # a new exon, they are added to the list of exons. Then, each of the
        # original parent exons are considered for deletion (variable to_keep).
        # All the exons are kept till the very end, when all exons labelled
        # to_keep == 0 are removed.
        #
        # Hence the number of exons could potentially grow quadratically.
        # Since we have to also consider merging the pairs of exons which are
        # newly created, the number of pairs considered can be larger than
        # m^2 where m is the number of original exons. However, I believe that
        # the algorithm should terminate because there are only finitely many
        # splice sites and hence at worst, we have all possible left and right
        # exon ends used in our splicegraph.
        #
        # I am not sure what the best termination condition for the loops are.
        #
        #####################################################################
        if genes[gene_idx].splicegraph.vertices.shape[1] == 0:
            continue
      
        genes[gene_idx].splicegraph.sort()
        vertices = genes[gene_idx].splicegraph.vertices.cope()
        edges = genes[gene_idx].splicegraph.edges.copy()
        to_keep = sp.ones((vertices.shape[1],))
      
        changed = False
        exon_idx = 0
        first_merge_exon = 0
        test_exon_idx = -1
        while exon_idx <= vertices.shape[1]:
      
            if changed:
              changed = False
              exon_idx = first_merge_exon
              continue
            
            num_exons = vertices.shape[1]
            ### sort exons ascending if neccessary
            exon_order = sp.argsort(vertices[0, :])
            if ~sp.all(vertices == vertices[:, exon_order]):
                vertices = vertices[:, exon_order]
                edges = edges[exon_order, :][:, exon_order]
                to_keep = to_keep[exon_order]
            
            ### take all exons overlapping to exon_idx +- MISALIGN window
            exon_take_idx = sp.where(((vertices[0, exon_idx] - 2*MISALIGN) < vertices[1, :]) & ((vertices[1, exon_idx] + 2 * MISALIGN) > vertices[0, :]))[0]
            if exon_take_idx.shape[0] > 0:
                first_merge_exon = exon_take_idx[0]
            ### consider only exons downstream of exon_idx
            exon_take_idx = exon_take_idx[sp.where(exon_take_idx > exon_idx)[0]]
            internal_idx = 0
      
            while (internal_idx <= exon_take_idx.shape[0]) and not changed:
                test_exon_idx = exon_take_idx[internal_idx]
                assert (test_exon_idx < exon_idx)
                cur_edge_left = sp.sum(edges[:exon_idx + 1, exon_idx]) > 0
                test_edge_left = sp.sum(edges[:test_exon_idx + 1, test_exon_idx]) > 0
                cur_edge_right = sp.sum(edges[exon_idx:, exon_idx]) > 0
                test_edge_right = sp.sum(edges[test_exon_idx:, test_exon_idx]) > 0
              
                new_vertex = sp.zeros((2, 1), dtype='int')
          
                #####################################################################
                # In the following, the cases are labelled with a binary string
                # [cur_edge_left, cur_edge_right, test_edge_left, test_edge_right]
                #
                # cur refers to exon_idx
                # test refers to test_exon_idx
                # edge means introns.
                # e.g. cur_edge_left == 1 means that there exists an intron on
                # the left of exon_idx
                #
                # Hence there are 16 cases, and we treat each one individually.
                # For each case, the two logical tests are:
                # 1. should we merge?
                # 2. is the new exon known?
                #
                # Note that a splice site is defined by an intron-exon boundary.
                # Furthermore, we allow misalignment by MISALIGN nucleotides,
                # that is we shorten an exonAendA if it doesn't have a splice site
                # and there exists a splice site exonBendB within MISALIGN nucleotides.
                # 
                #####################################################################
                
                
                # 0000
                if not cur_edge_left and not cur_edge_right and not test_edge_left and not test_edge_right:
                  ### if exon_idx and test_exon_idx differ in at least one border and
                  ### at least one of both is not yet deleted and 
                  ###   --- exon ---X               X >=Y
                  ###         Y--- test_exon ---  
                  if (~sp.all(vertices[0, exon_idx] == vertices[0, test_exon_idx]) or ~sp.all(vertices[1, exon_idx] == vertices[1, test_exon_idx])) and \
                     (to_keep[test_exon_idx] != 0 or to_keep[exon_idx] != 0) and \
                     (vertices[1, exon_idx] >= vertices[0, test_exon_idx]):
                    
                      new_vertex[0] = min(vertices[0, exon_idx], vertices[0, test_exon_idx])
                      new_vertex[1] = max(vertices[1, exon_idx], vertices[1, test_exon_idx])
                    
                      known = 0
                      for ix in range(vertices.shape[1]):
                          if sp.all(new_vertex == vertices[:, ix]):
                              known = ix
                              break
                    
                      if known == 0:
                          vertices = sp.c_[vertices, new_vertex]
                          edges = sp.c_[edges, sp.zeros((edges.shape[1], 1), dtype='int')]
                          edges = sp.r_[edges, sp.zeros((1, edges.shape[0]), dtype='int')]
                          to_keep = sp.r_[to_keep, 1]
                          to_keep[exon_idx] = 0
                          to_keep[test_exon_idx] = 0
                          changed = True
                      else:
                          if known != exon_idx and to_keep[exon_idx] > 0:
                              to_keep[exon_idx] = 0
                              changed = True
                          if known != test_exon_idx and to_keep[test_exon_idx] > 0:
                              to_keep[test_exon_idx] = 0
                              changed = True
                  
                # 0001
                elif not cur_edge_left and not cur_edge_right and not test_edge_left and test_edge_right:
                    ###   --- exon ---M
                    ###         --- test_exon ---<
                    if (vertices[1, exon_idx] - MISALIGN <= vertices[1, test_exon_idx]) and \
                       (vertices[1, exon_idx] >= vertices[0, test_exon_idx]):
                  
                        new_vertex[0] = min(vertices[0, exon_idx], vertices[0, test_exon_idx])
                        new_vertex[1] = vertices[1, test_exon_idx]
                        new_edges = (edges[exon_idx, :] | edges[test_exon_idx, :])
                      
                        known = 0
                        for ix in range(vertices.shape[1]):
                            if sp.all(new_vertex == vertices[:, ix]) and issubset(new_edges, edges[ix, :]):
                                known = ix
                                break
                        assert(known != exon_idx)
                      
                        if to_keep[exon_idx] > 0:
                            to_keep[exon_idx] = 0
                            changed = True
                      
                        if known == 0:
                            vertices = sp.c_[vertices, new_vertex]
                            edges = sp.c_[edges, new_edges]
                            edges = sp.r_[edges, sp.r_[new_edges, 0]]
                            to_keep = sp.r_[to_keep, 1]
                            to_keep[exon_idx] = 0
                            to_keep[test_exon_idx] = 0
                            changed = True
                        elif known != test_exon_idx and to_keep[test_exon_idx] > 0:
                            to_keep[test_exon_idx] = 0
                            changed = True
      
                # 0100
                elif not cur_edge_left and cur_edge_right and not test_edge_left and not test_edge_right:
                    ###            --- exon ---<
                    ###  --- test_exon ---M
                    if (vertices[1, exon_idx] >= vertices[1, test_exon_idx] - MISALIGN) and \
                       (vertices[0, exon_idx] <= vertices[1, test_exon_idx]):
                  
                        new_vertex[0] = min(vertices[0, exon_idx], vertices[0, test_exon_idx])
                        new_vertex[1] = vertices[1, exon_idx]
                        new_edges = (edges[exon_idx, :] | edges[test_exon_idx, :])
                  
                        known = 0
                        for ix in range(vertices.shape[1]):
                            if sp.all(new_vertex == vertices[:, ix]) and issubset(new_edges, edges[ix, :]):
                                known = ix
                                break
                        assert(known != test_exon_idx)
                  
                        if to_keep[test_exon_idx] > 0:
                            to_keep[test_exon_idx] = 0
                            changed = True
                  
                        if known == 0:
                            vertices = sp.c_[vertices, new_vertex]
                            edges = sp.c_[edges, new_edges]
                            edges = sp.r_[edges, sp.r_[new_edges, 0]]
                            to_keep = sp.r_[to_keep, 1]
                            to_keep[exon_idx] = 0
                            to_keep[test_exon_idx] = 0
                            changed = True
                        elif known != exon_idx and to_keep[exon_idx] > 0:
                            to_keep[exon_idx] = 0
                            changed = True
          
                # 0010
                elif not cur_edge_left and not cur_edge_right and test_edge_left and not test_edge_right:
                    ###            M--- exon ---
                    ### >--- test_exon ---
                    if (vertices[0, exon_idx] + MISALIGN >= vertices[0, test_exon_idx]) and \
                       (vertices[0, exon_idx] <= vertices[1, test_exon_idx]):
                        
                        new_vertex[0] = vertices[0, test_exon_idx]
                        new_vertex[1] = max(vertices[1, exon_idx], vertices[1, test_exon_idx])
                        new_edges = (edges[exon_idx, :] | edges[test_exon_idx, :])
                      
                        known = 0
                        for ix in range(vertices.shape[1]):
                            if sp.all(new_vertex == vertices[:, ix]) and issubset(new_edges, edges[ix, :]):
                                known = ix
                                break
                        assert(known != exon_idx)
                      
                        if to_keep[exon_idx] > 0:
                            to_keep[exon_idx] = 0
                            changed = True
                      
                        if known == 0:
                            vertices = sp.c_[vertices, new_vertex]
                            edges = sp.c_[edges, new_edges]
                            edges = sp.r_[edges, sp.r_[new_edges, 0]]
                            to_keep = sp.r_[to_keep, 1]
                            to_keep[exon_idx] = 0
                            to_keep[test_exon_idx] = 0
                            changed = True
                        elif known != test_exon_idx and to_keep[test_exon_idx] > 0:
                            to_keep[test_exon_idx] = 0
                            changed = True
          
                # 1000
                elif cur_edge_left and not cur_edge_right and not test_edge_left and not test_edge_right:
                    ###  >--- exon ---
                    ###      M--- test_exon ---
                    if (vertices[0, exon_idx] <= vertices[0, test_exon_idx] + MISALIGN) and \
                       (vertices[1, exon_idx] >= vertices[0, test_exon_id]):
                  
                        new_vertex[0] = vertices[0, exon_idx]
                        new_vertex[1] = max(vertices[1, exon_idx], vertices[1, test_exon_idx])
                        new_edges = (edges[exon_idx, :] | edges[test_exon_idx, :])
                      
                        known = 0
                        for ix in range(vertices.shape[1]):
                            if sp.all(new_vertex == vertices[:, ix]) and issubset(new_edges, edges[ix, :]):
                                known = ix
                                break
                        assert(known != test_exon_idx)
                      
                        if to_keep[test_exon_idx] > 0:
                            to_keep[test_exon_idx] = 0
                            changed = True
                      
                        if known == 0:
                            vertices = sp.c_[vertices, new_vertex]
                            edges = sp.c_[edges, new_edges]
                            edges = sp.r_[edges, sp.r_[new_edges, 0]]
                            to_keep = sp.r_[to_keep, 1]
                            to_keep[exon_idx] = 0
                            to_keep[test_exon_idx] = 0
                            changed = Tru
                        elif known != exon_idx and to_keep[exon_idx]:
                            to_keep[exon_idx] = 0
                            changed = True
                
                # 0011
                elif not cur_edge_left and not cur_edge_right and test_edge_left and test_edge_right:
                    ###     M--- exon ---M
                    ###  >--- test_exon ---<
                    if (vertices[0, exon_idx] + MISALIGN >= vertices[0, test_exon_idx]) and \
                       (vertices[1, exon_idx] - MISALIGN <= vertices[1, test_exon_idx]):
                  
                        new_vertex = vertices[:, test_exon_idx]
                        new_edges = (edges[exon_idx, :] | edges[test_exon_idx, :])
                        
                        known = 0
                        for ix in range(vertices.shape[1]):
                            if sp.all(new_vertex == vertices[:, ix]) and issubset(new_edges, edges[ix, :]):
                                known = ix
                                break
                        assert(known != exon_idx)
                      
                        if to_keep[test_exon_idx] > 0:
                            to_keep[test_exon_idx] = 0
                            changed = True
                      
                        if known == 0:
                            vertices = sp.c_[vertices, new_vertex]
                            edges = sp.c_[edges, new_edges]
                            edges = sp.r_[edges, sp.r_[new_edges, 0]]
                            to_keep = sp.r_[to_keep, 1]
                            to_keep[exon_idx] = 0
                            to_keep[test_exon_idx] = 0
                            changed = True
                        elif known != test_exon_idx and to_keep[test_exon_idx]:
                            to_keep[test_exon_idx] = 0
                            changed = True
      
                # 1100
                elif cur_edge_left and cur_edge_right and not test_edge_left and not test_edge_right:
                    ###    >-------- exon ---------<
                    ###       M--- test_exon ---M
                    if (vertices[0, exon_idx] <= vertices[0, test_exon_idx] + MISALIGN) and \
                       (vertices[1, exon_idx] >= vertices[1, test_exon_idx] - MISALIGN):
                  
                        new_vertex = vertices[:, exon_idx]
                        new_edges = (edges[exon_idx, :] | edges[test_exon_idx, :])
                        
                        known = 0
                        for ix in range(vertices.shape[1]):
                            if sp.all(new_vertex == vertices[:, ix]) and issubset(new_edges, edges[ix, :]):
                                known = ix
                                break
                        assert(known != test_exon_idx)
                      
                        if to_keep[test_exon_idx] > 0:
                            to_keep[test_exon_idx] = 0
                            changed = True
                      
                        if known == 0:
                            vertices = sp.c_[vertices, new_vertex]
                            edges = sp.c_[edges, new_edges]
                            edges = sp.r_[edges, sp.r_[new_edges, 0]]
                            to_keep = sp.r_[to_keep, 1]
                            to_keep[exon_idx] = 0
                            to_keep[test_exon_idx] = 0
                            changed = True
                        elif known != exon_idx and to_keep[exon_idx]:
                            to_keep[exon_idx] = 0
                            changed = True
          
                # 0101
                elif not cur_edge_left and cur_edge_right and not test_edge_left and test_edge_right:
                    ###    -------- exon ------<
                    ###       --- test_exon ---<
                    if vertices[1, exon_idx] == vertices[1, test_exon_idx]:
                        ### same edges
                        if not sp.all(edges[exon_idx, :] == edges[test_exon_idx, :]):
                            edges[exon_idx, :] = (edges[exon_idx, :] | edges[test_exon_idx, :])
                            edges[:, exon_idx] = (edges[:, exon_idx] | edges[:, test_exon_idx])
                            edges[test_exon_idx, :] = (edges[exon_idx, :] | edges[test_exon_idx, :])
                            edges[:, test_exon_idx] = (edges[:, exon_idx] | edges[:, test_exon_idx])
                            ### keep longer exon
                            if vertices[0, exon_idx] > vertices[0, test_exon_idx]:
                                to_keep[exon_idx] = 0
                            else:
                                to_keep[test_exon_idx] = 0
                            changed = True
                        ### different edges
                        else:
                            ###        ---- exon ------<
                            ###   ------- test_exon ---<
                            if (vertices[0, exon_idx] > vertices[0, test_exon_idx]) and to_keep[exon_idx] > 0:
                                to_keep[exon_idx] = 0
                                changed = True
                            ###   --------- exon ------<
                            ###       --- test_exon ---<
                            if (vertices[0, exon_idx] <= vertices[0, test_exon_idx]) and to_keep[test_exon_idx] > 0:
                                to_keep[test_exon_idx] = 0
                                changed = True
          
                # 1010
                elif cur_edge_left and not cur_edge_right and test_edge_left and not test_edge_right:
                    ###    >-------- exon ------
                    ###    >--- test_exon ---
                    if vertices[0, exon_idx] == vertices[0, test_exon_idx]:
                        ### same edges
                        if not sp.all(edges[exon_idx, :] == edges[test_exon_idx, :]):
                            edges[exon_idx, :] = (edges[exon_idx, :] | edges[test_exon_idx, :])
                            edges[:, exon_idx] = (edges[:, exon_idx] | edges[:, test_exon_idx])
                            edges[test_exon_idx, :] = (edges[exon_idx, :] | edges[test_exon_idx, :])
                            edges[:, test_exon_idx] = (edges[:, exon_idx] | edges[:, test_exon_idx])
                            ### keep longer exon
                            if vertices[1, exon_idx] < vertices[1, test_exon_idx]:
                                to_keep[exon_idx] = 0
                            else:
                                to_keep[test_exon_idx] = 0
                            changed = True
                        ### different edges
                        else:
                            ###    >----- exon ------
                            ###    >------ test_exon ---
                            if (vertices[1, exon_idx] < vertices[1, test_exon_idx]) and to_keep[exon_idx] > 0:
                                to_keep[exon_idx] = 0
                                changed = 1
                            ###    >-------- exon ------
                            ###    >--- test_exon ---
                            if (vertices[1, exon_idx] >= vertices[1,test_exon_idx]) and to_keep[test_exon_idx] > 0:
                                to_keep[test_exon_idx] = 0
                                changed = 1
          
                # 0110  
                elif not cur_edge_left and cur_edge_right and test_edge_left and not test_edge_right:
                    ###          MX----- exon -----<    X<=Y
                    ###    >--- test_exon ---YM
                    if (vertices[0, exon_idx] + MISALIGN >= vertices[0, test_exon_idx]) and \
                       (vertices[1, exon_idx] >= vertices[1, test_exon_idx] - MISALIGN) and \
                       (vertices[0, exon_idx] <= vertices[1, test_exon_idx]):
      
                        new_vertex = sp.array([[vertices[0, test_exon_idx]], [vertices[1, exon_idx]]])
                        new_edges = (edges[exon_idx, :] | edges[test_exon_idx, :])
                  
                        known = 0
                        for ix in range(vertices.shape[1]):
                            # assert(isequal(find(new_edges),intersect(find(new_edges),find(edges(ix,:)))) == ~any(new_edges & ~edges(ix,:)))
                            if sp.all(new_vertex == vertices[:, ix]) and not sp.any(new_edges & (edges[ix, :] == 0)):
                                known = ix
                                break
                        assert((known != exon_idx) and (known != test_exon_idx))
                      
                        if known == 0:
                            vertices = sp.c_[vertices, new_vertex]
                            edges = sp.c_[edges, new_edges]
                            edges = sp.r_[edges, sp.r_[new_edges, 0]]
                            to_keep = sp.r_[to_keep, 1]
                            to_keep[exon_idx] = 0
                            to_keep[test_exon_idx] = 0
                            changed = True
                        elif to_keep[exon_idx] > 0 or to_keep[test_exon_idx] > 0:
                            to_keep[exon_idx] = 0
                            to_keep[test_exon_idx] = 0
                            changed = True
      
                # 1001
                elif cur_edge_left and not cur_edge_right and not test_edge_left and test_edge_right:
                    ###   >----- exon -----XM            X>=Y
                    ###          MY--- test_exon ---<
                    if (vertices[0, exon_idx] <= vertices[0, test_exon_idx] + MISALIGN) and \
                       (vertices[1, exon_idx] - MISALIGN <= vertices[1, test_exon_idx]) and \
                       (vertices[1, exon_idx] >= vertices[0, test_exon_idx]):
      
                        new_vertex = sp.array([[vertices[0, exon_idx]], [vertices[1, test_exon_idx]]])
                        new_edges = (edges[exon_idx, :] | edges[test_exon_idx, :])
                        known = 0
                        for ix in range(vertices.shape[1]):
                            # assert(isequal(find(new_edges),intersect(find(new_edges),find(edges(ix,:))))==~any(new_edges & ~edges(ix,:)))
                            if sp.all(new_vertex == vertices[:, ix]) and not sp.any(new_edges & (edges[ix, :] == 0)):
                                known = ix
                                break
                        assert((known != exon_idx) and (known != test_exon_idx))
                      
                        if known == 0:
                            vertices = sp.c_[vertices, new_vertex]
                            edges = sp.c_[edges, new_edges]
                            edges = sp.r_[edges, sp.r_[new_edges, 0]]
                            to_keep = sp.r_[to_keep, 1]
                            to_keep[exon_idx] = 0
                            to_keep[test_exon_idx] = 0
                            changed = True
                        elif to_keep[exon_idx] > 0 or to_keep[test_exon_idx] > 0:
                            to_keep[exon_idx] = 0
                            to_keep[test_exon_idx] = 0
                            changed = True
          
                # 0111 
                elif not cur_edge_left and cur_edge_right and test_edge_left and test_edge_right:
                    ###        ----- exon -----<           
                    ###      >--- test_exon ---<
                    if vertices[1, exon_idx] == vertices[1, test_exon_idx]:
                        #    (vertices(1,exon_idx)+MISALIGN>=vertices(1,test_exon_idx))
                        # index to the right
                        right_edge = max(vertices[1, exon_idx], vertices[1, test_exon_idx])
                        idx = sp.argmax(vertices[0, :] > right_edge)
                      
                        if not isequal(edges[exon_idx, idx:], edges[test_exon_idx, idx:]):
                          
                            edges[exon_idx, idx:] = (edges[exon_idx, idx:] | edges[test_exon_idx, idx:])
                            edges[test_exon_idx, idx:] = edges[exon_idx, idx:]
                            edges[idx:, exon_idx] = (edges[exon_idx, idx:] | edges[test_exon_idx, idx:])
                            edges[idx:, test_exon_idx] = edges[idx:, exon_idx]
                          
                            #to_keep[exon_idx] = 0
                            changed = True
                        ###        M---- exon -----<           
                        ###      >--- test_exon ---<
                        if (vertices[0, exon_idx] + MISALIGN >= vertices[0, test_exon_idx]) and to_keep[exon_idx] > 0:
                            to_keep[exon_idx] = 0
                            changed = True
          
                # 1101
                elif cur_edge_left and cur_edge_right and not test_edge_left and test_edge_right:
                    ###        >------ exon -----<           
                    ###         --- test_exon ---<
                    if vertices[1, exon_idx] == vertices[1, test_exon_idx]:
                        #    (vertices(1,exon_idx)-MISALIGN<=vertices(1,test_exon_idx))
                        # index to the right
                        right_edge = max(vertices[1, exon_idx], vertices[1, test_exon_idx])
                        idx = sp.argmax(vertices[0, :] > right_edge)
                  
                        if not isequal(edges[exon_idx, idx:], edges[test_exon_idx, idx:]):
                            
                            edges[exon_idx, idx:] = (edges[exon_idx, idx:] | edges[test_exon_idx, idx:])
                            edges[test_exon_idx, idx:] = edges[exon_idx, idx:]
                            edges[idx:, exon_idx] = (edges[exon_idx, idx:] | edges[test_exon_idx, idx:])
                            edges[idx:, test_exon_idx] = edges[idx:, exon_idx]
                          
                            #to_keep[test_exon_idx] = 0
                            changed = True
                        ###    M-------- exon -----<           
                        ###      >--- test_exon ---<
                        if (vertices[0, exon_idx] - MISALIGN <= vertices[0, test_exon_idx]) and to_keep[test_exon_idx] > 0:
                            to_keep[test_exon_idx] = 0
                            changed = True
      
                # 1011  
                elif cur_edge_left and not cur_edge_right and test_edge_left and test_edge_right:
                    ###      >------ exon ---           
                    ###      >--- test_exon ---<
                    if vertices[0, exon_idx] == vertices[0, test_exon_idx]:
                        #(vertices(2,exon_idx)-MISALIGN<=vertices(2,test_exon_idx))
                        # index to the left
                        left_edge = min(vertices[0, exon_idx], vertices[0, test_exon_idx])
                        idx = sp.argmax(vertices[1, :] < left_edge)
                      
                        if not isequal(edges[exon_idx, :idx], edges[test_exon_idx, :idx]):
      
                            edges[exon_idx, :idx] = (edges[exon_idx, :idx] | edges[test_exon_idx, :idx])
                            edges[test_exon_idx, :idx] = edges[exon_idx, :idx]
                            edges[:idx, exon_idx] = (edges[exon_idx, :idx] | edges[test_exon_idx, :idx])
                            edges[:idx, test_exon_idx] = edges[:idx, exon_idx]
                        
                            #to_keep[exon_idx] = 0
                            changed = True
                        ###     >------ exon ---M           
                        ###     >--- test_exon ----<
                        if (vertices[1, exon_idx] - MISALIGN <= vertices[1, test_exon_idx]) and to_keep[exon_idx] > 0:
                            to_keep[exon_idx] = 0
                            changed = True
      
                # 1110
                elif cur_edge_left and cur_edge_right and test_edge_left and not test_edge_right:
                    ###      >------ exon ------<           
                    ###      >--- test_exon ---
                    if vertices[0, exon_idx] == vertices[0, test_exon_idx]:
                        #    (vertices(2,exon_idx)+MISALIGN>=vertices(2,test_exon_idx))
                        # index to the left
                        left_edge = min(vertices[0, exon_idx], vertices[0, test_exon_idx])
                        idx = sp.argmax(vertices[1, :] < left_edge)
                  
                        if not isequal(edges[exon_idx, :idx], edges[test_exon_idx, :idx]):
                            edges[exon_idx, :idx] = (edges[exon_idx, :idx] | edges[test_exon_idx, :idx])
                            edges[test_exon_idx, :idx] = edges[exon_idx, :idx]
                            edges[:idx, exon_idx] = (edges[exon_idx, :idx] | edges[test_exon_idx, :idx])
                            edges[:idx, test_exon_idx] = edges[:idx, exon_idx]
                      
                            #to_keep[test_exon_idx] = 0
                            changed = True
                        ###      >------ exon -------<           
                        ###      >--- test_exon ---M
                        #if (vertices(2,exon_idx) >= vertices(2,test_exon_idx)-MISALIGN) && (to_keep[test_exon_idx])
                        if (vertices[1, exon_idx] + MISALIGN >= vertices[1, test_exon_idx]) and to_keep[test_exon_idx] > 0:
                            to_keep[test_exon_idx] = 0
                            changed = True
      
                # 1111
                elif cur_edge_left and cur_edge_right and test_edge_left and test_edge_right:
                  ###      >------ exon ------<           
                  ###      >--- test_exon ----<
                  if (vertices[0, exon_idx] == vertices[0, test_exon_idx]) and \
                     (vertices[1, exon_idx] == vertices[1, test_exon_idx]):
                      if not isequal(edges[exon_idx, :], edges[test_exon_idx, :]):
                          edges[exon_idx, :] = (edges[exon_idx, :] | edges[test_exon_idx, :])
                          edges[test_exon_idx, :] = edges[exon_idx, :]
                          edges[:, exon_idx] = (edges[exon_idx, :] | edges[test_exon_idx, :])
                          edges[:, test_exon_idx] = edges[:, exon_idx]
                    
                          to_keep[test_exon_idx] = 0
                          changed = True
                      elif to_keep[test_exon_idx] > 0:
                          to_keep[test_exon_idx] = 0
                          changed = True
                  ###      >------ exon ----<   OR   >------ exon ------<          
                  ###      >--- test_exon ----<      >--- test_exon --<
                  elif (vertices[0, exon_idx] == vertices[0, test_exon_idx]) and \
                       (vertices[1, exon_idx] != vertices[1, test_exon_idx]) and \
                       (to_keep[exon_idx] > 0 or to_keep[test_exon_idx] > 0):
                    # index to the left
                    #left_edge = min(vertices(1,exon_idx),vertices(1,test_exon_idx))
                    #idx = find(vertices(2,:)<left_edge, 1, 'last')
                    idx = sp.argmax(vertices[1, :] < vertices[0, exon_idx])
      
                    if not isequal(edges[exon_idx, :idx], edges[test_exon_idx, :idx]):
                        edges[exon_idx, :idx] = (edges[exon_idx, :idx] | edges[test_exon_idx, :idx])
                        edges[test_exon_idx, :idx] = edges[exon_idx, :idx]
                        edges[:idx, exon_idx] = (edges[exon_idx, :idx] | edges[test_exon_idx, :idx])
                        edges[:idx, test_exon_idx] = edges[:idx, exon_idx]
                            
                        changed = True
                  ###      >------ exon ----<   OR   >------ exon ------<          
                  ###    >--- test_exon ----<          >--- test_exon --<
                  elif (vertices[0, exon_idx] != vertices[0, test_exon_idx]) and \
                       (vertices[1, exon_idx] == vertices[1, test_exon_idx]) and \
                       (to_keep[exon_idx] > 0 or to_keep[test_exon_idx] > 0):
                    # index to the right
                    #right_edge = max(vertices(2,exon_idx),vertices(2,test_exon_idx))
                    #idx = find(vertices(1,:)>right_edge,1,'first')
                    idx = sp.argmax(vertices[0, :] > vertices[1, exon_idx])
      
                    if not isequal(edges[exon_idx, idx:], edges[test_exon_idx, idx:]):
      
                        edges[exon_idx, idx:] = (edges[exon_idx,idx:] | edges[test_exon_idx, idx:])
                        edges[test_exon_idx, idx:] = edges[exon_idx, idx:]
                        edges[idx:, exon_idx] = (edges[exon_idx,idx:] | edges[test_exon_idx, idx:])
                        edges[idx:, test_exon_idx] = edges[idx:, exon_idx]
                    
                        changed = True
                else:
                    raise Exception('Unknown case!')
                #test_exon_idx = test_exon_idx+1
                internal_idx = internal_idx + 1
            exon_idx += 1
      
        genes[gene_idx].splicegraph.vertices = vertices.copy()
        genes[gene_idx].splicegraph.edges = edges.copy()
        genes[gene_idx].splicegraph.subset(sp.where(to_keep)[0])
        genes[gene_idx].splicegraph.update_terminals()

    return genes
