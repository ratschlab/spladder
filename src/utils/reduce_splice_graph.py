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
