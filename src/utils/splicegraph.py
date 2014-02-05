import scipy as sp

class splicegraph:
    
    def __init__(self, genes = None):
        
        vertices = sp.zeros((2, 0))
        edges = sp.zeros((0, 0))
        terminals = sp.zeros((2, 0))

        if genes:
            from_genes(genes)

    def new_edge(self):

        edges = sp.c_[edges, sp.zeros((edges.shape[0],))]
        edges = sp.r_[edges, sp.zeros((edges.shape[1],))]

    def from_genes(self, genes):
        
        if gene_idx % 100 == 0:
            print '.',
            if gene_idx % 1000 == 0:
                 print '%i' % gene_idx

        for transcript_idx = len(genes[gene_idx].transcripts):
            exon_start_end = genes[gene_idx].exons[transcript_idx]
            
            ### only one exon in the transcript
            if exon_start_end.shape[0] == 1,
                exon1_start = exon_start_end[0, 0]
                exon1_end = exon_start_end[0, 1]

                if vertices.shape[1] == 0:
                  vertices[0, 0] = exon1_start
                  vertices[1, 0] = exon1_end
                  edges = 0
                  num_exons = 0
                else:
                  vertices[0, num_exons + 1) = exon1_start
                  vertices[1, num_exons + 1) = exon1_end
                  edges[0, num_exons + 1] = 0
                  edges[num_exons + 1, 0] = 0
                  num_exons += 1
            ### more than one exon in the transcript
            else:
                for exon_idx in  xrange(exon_start_end.shape[0] - 1):
                    exon1_start = exon_start_end[exon_idx , 0]
                    exon1_end = exon_start_end[exon_idx, 1]
                    exon2_start = exon_start_end[exon_idx + 1, 0]
                    exon2_end = exon_start_end[exon_idx + 1, 1]
          
                    if vertices.shape[1] == 0:
                        vertices[0, 0] = exon1_start
                        vertices[1, 0] = exon1_end
                        vertices[0, 1] = exon2_start
                        vertices[1, 1] = exon2_end
                        edges = sp.zeros((2, 2))
                        edges[0, 1] = 1
                        edges[1, 0] = 1
                        num_exons = 2
                    else:
                        exon1_idx = 0
                        exon2_idx = 0
                        ### check if current exon already occurred
                        for idx in range(num_exons):
                            if ((vertices[0, idx] == exon1_start) && (vertices[1, idx] == exon1_end)):
                                 exon1_idx = idx
                            if ((vertices[0, idx] == exon2_start) && (vertices[1, idx] == exon2_end)):
                                 exon2_idx = idx

                        ### both exons already occured -> only add an edge
                        if (exon1_idx !=0) && (exon2_idx != 0):
                            edges[exon1_idx, exon2_idx] = 1
                            edges[exon2_idx, exon1_idx] = 1
                        else:
                            ### 2nd exon occured
                            if ((exon1_idx == 0) && (exon2_idx !=0 )):
                                vertices[0, num_exons + 1] = exon1_start
                                vertices[1, num_exons + 1] = exon1_end
                                edges[exon2_idx, num_exons + 1] = 1
                                edges[num_exons + 1, exon2_idx] = 1
                                num_exons += 1
                            ### 1st exon occured
                            elif ((exon2_idx == 0) && (exon1_idx !=0)):
                                vertices[0, num_exons + 1] = exon2_start
                                vertices[1, num_exons + 1] = exon2_end
                                edges[exon1_idx, num_exons + 1] = 1
                                edges[num_exons + 1, exon1_idx] = 1
                                num_exons += 1
                            ### no exon occured
                            else:
                                assert((exon1_idx == 0) && (exon2_idx == 0))
                                vertices[0, num_exons + 1] = exon1_start
                                vertices[1, num_exons + 1] = exon1_end
                                num_exons += 1
                                vertices[0, num_exons + 1] = exon2_start
                                vertices[1, num_exons + 1] = exon2_end
                                num_exons += 1    
                              
                                edges[num_exons - 1, num_exons] = 1
                                edges[num_exons, num_exons - 1] = 1

                terminals[0, sp.where(vertices[0, :] == exon_start_end.max())[0]] = 1
                terminals[1, sp.where(vertices[1, :] == exon_start_end.max())[0]] = 1

          ### take care of the sorting by exon start
          s_idx = sp.argsort(vertices[0, :])
          vertices = vertices[:, s_idx]
          edges = edges[s_idx, :][:, s_idx]
          if vertices.shape != terminals.shape:
            terminals = sp.c_[terminals, sp.zeros((2,))]
          terminals = terminals[:, s_idx]

    def add_intron(self, idx1, flag1, idx2, flag2):
        """adds new introns into splicegraph between idx1 and idx2"""

        ### if flag1, all end terminal exons in idx1 are preserved
        ### if flag2, all start terminal exons in idx2 are preserved

        if idx2:
            adj_mat = sp.triu(edges)

            if flag1:
                for i1 in idx1:

                    ### if exon is end-terminal
                    if sp.all(adj_mat[i1, :] == 0):

                        vertices = sp.c_[vertices, vertices[:, i1]]

                        _new_edge()
                        edges[:, -1] = edges[:, i1]
                        edges[-1, :] = edges[i1, :]

                        terminals = sp.c_[terminals, terminals[:, i1]]
            if flag2:,
                for i2 in idx2:
                    ### if exon is start-terminal
                    if sp.all(adj_mat[:, i2] == 0):
                        vertices = sp.c_[vertices, vertices[:, i2]]

                        _new_edge()
                        edges[:, -1] = edges[:, i2]
                        edges[-1, :] = edges[i2, :]

                        terminals = sp.c_[terminals, terminals[:, i2]]

        for i1 in idx1:
            for i2 in idx2:
                edges[i1, i2] = 1
                edges[i2, i1] = 1

    def add_cassette_exon(self, new_exon, exons_pre, exons_aft)
        ### exon_pre contains the indices of preceding exons
        ### exon_aft contains the indices of successing exons
        
        vertices = sp.r_[vertices, new_exon.T]

        _new_edge()

        edges[exons_pre, -1] = 1
        edges[exons_aft, -1] = 1
        edges[-1, :] = edges[:, -1].T

        terminals = sp.c_[terminals, sp.zeros((2,))]


    def add_intron_retention(self, idx1, idx2):

        adj_mat = sp.triu(edges)

        vertices = sp.c_[vertices, sp.array([vertices[0, idx1], vertices[1, idx2]])]

        _new_edge()

        adj_mat = sp.r_[adj_mat, sp.zeros((adj_mat.shape[1],))]
        adj_mat = sp.c_[adj_mat, sp.zeros((adj_mat.shape[2],))]

        ### check if adjacency matrix is symmetric
        ### otherwise or is not justyfied
        assert(sp.all(sp.all(adj_mat - (splicegraph['graph'] - adj_mat).T == 0)))

        ### AK: under the assumption that our splice graph representation is symmetric
        ### I preserve symmetry by using OR over the adj_mat column and row
        
        splicegraph['graph'] = adj_mat[:, idx1] | adj_mat[idx2, :].T
        splicegraph['graph'] = adj_mat[:, idx1].T | adj_mat[idx2, :]

        if 'term' in splicegraph:
            splicegraph['term'] = sp.c_[splicegraph['term'], sp.array([splicegraph['term'][1, idx1], splicegraph['term'][2, idx2]])]


