def infer_splice_graph_caller(genes):
    # genes = infer_splice_graph_caller(genes)
      
    MISALIGN = 5

    print 'Number of genes:%d' % genes.shape[0]

    #####################################################################
    # Some simple merging to reduce exons
    #####################################################################
    print 'Merging some genes'
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

    print 'performing inference'

    for gene_idx in range(genes.shape[0]):

        if gene_idx % 1000 ==0:
            print '%d ' % gene_idx
      
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
                if (test_exon_idx < exon_idx):
                    pdb.set_trace()
                cur_edge_left = sp.sum(edges[:exon_idx + 1, exon_idx]) > 0
                test_edge_left = sp.sum(edges[:test_exon_idx + 1, test_exon_idx]) > 0
                cur_edge_right = sp.sum(edges[exon_idx:, exon_idx]) > 0
                test_edge_right = sp.sum(edges[test_exon_idx:, test_exon_idx]) > 0
              
                new_vertex = sp.zeros((2, 1))
          
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
                  if (~sp.all(vertices[0, exon_idx] == vertices[0, test_exon_idx]) or ~sp.all(vertices[1, exon_idx] == vertices[1, test_exon_idx])) and 
                     (to_keep[test_exon_idx] != 0 or to_keep[exon_idx] != 0) and
                     (vertices[1, exon_idx] >= vertices[0, test_exon_idx]):
                    
                      new_vertex[0] = min([vertices[0, exon_idx], vertices[0, test_exon_idx])
                      new_vertex[1] = max([vertices[1, exon_idx], vertices[1, test_exon_idx])
                    
                      known = 0
                      for ix in range(vertices.shape[1]):
                          if sp.all(new_vertex == vertices[:, ix]):
                              known = ix
                              break
                    
                      if known == 0:
                          vertices = sp.c_[vertices, new_vertex]
                          edges = sp.c_[edges, sp.zeros((edges.shape[1], 1))]
                          edges = sp.r_[edges, sp.zeros((1, edges.shape[0]))]
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
                    if (vertices[1, exon_idx] - MISALIGN <= vertices[1, test_exon_idx]) and
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
                    if (vertices[1, exon_idx] >= vertices[1, test_exon_idx] - MISALIGN) and 
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
                    if (vertices[0, exon_idx] + MISALIGN >= vertices[0, test_exon_idx]) and
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
                    if (vertices[0, exon_idx] <= vertices[0, test_exon_idx] + MISALIGN) and 
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
                    if (vertices[0, exon_idx] + MISALIGN >= vertices[0, test_exon_idx]) and
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
                    if (vertices[0, exon_idx] <= vertices[0, test_exon_idx] + MISALIGN) and 
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
                        else
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
                    if (vertices[0, exon_idx] + MISALIGN >= vertices[0, test_exon_idx]) and
                       (vertices[1, exon_idx] >= vertices[1, test_exon_idx] - MISALIGN) and 
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
                    if (vertices[0, exon_idx] <= vertices[0, test_exon_idx] + MISALIGN) and
                       (vertices[1, exon_idx] - MISALIGN <= vertices[1, test_exon_idx]) and
                       (vertices[1, exon_idx] >= vertices[0, test_exon_idx]):
      
                        new_vertex = sp.array([[vertices[0,exon_idx)]], [vertices[1, test_exon_idx]]])
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
                  if (vertices[0, exon_idx] == vertices[0, test_exon_idx]) and 
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
                  elif (vertices[0, exon_idx] == vertices[0, test_exon_idx]) and 
                       (vertices[1, exon_idx] != vertices[1, test_exon_idx]) and
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
                  elif (vertices[0, exon_idx] != vertices[0, test_exon_idx]) and
                       (vertices[1, exon_idx] == vertices[1, test_exon_idx]) and
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
                else
                    print >> sys.stderr, 'Unknown case!!'
                    pdb.set_trace()
                #test_exon_idx = test_exon_idx+1
                internal_idx = internal_idx + 1
            exon_idx += 1
      
        genes[gene_idx].splicegraph.vertices = vertices.copy()
        genes[gene_idx].splicegraph.edges = edges.copy()
        genes[gene_idx].splicegraph.subset(sp.where(to_keep)[0])
        genes[gene_idx].splicegraph.update_terminals()

    return genes
