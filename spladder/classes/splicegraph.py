import scipy as sp

from ..utils import *
from ..init import get_tags_gtf

class Splicegraph:

    def __init__(self, gene = None):
        
        self.vertices = sp.zeros((2, 0), dtype='int')
        self.edges = sp.zeros((0, 0), dtype='int')
        self.terminals = sp.zeros((2, 0), dtype='int')

        if gene:
            self.from_gene(gene)

    def get_len(self):
        
        return self.vertices.shape[1]

    def new_edge(self):

        self.edges = sp.c_[self.edges, sp.zeros((self.edges.shape[0], 1), dtype='int')]
        self.edges = sp.r_[self.edges, sp.zeros((1, self.edges.shape[1]), dtype='int')]
    
    def subset(self, keep_idx):
        
        self.vertices = self.vertices[:, keep_idx]
        self.edges = self.edges[keep_idx, :][:, keep_idx]
        self.terminals = self.terminals[:, keep_idx]

    def update_terminals(self):
        
        self.terminals = sp.zeros(self.vertices.shape, dtype='int')
        self.terminals[0, sp.where(sp.sum(sp.tril(self.edges), axis=1) == 0)[0]] = 1
        self.terminals[1, sp.where(sp.sum(sp.triu(self.edges), axis=1) == 0)[0]] = 1

    def reorder(self, idx):
        
        self.vertices = self.vertices[:, idx]
        self.edges = self.edges[idx, :][:, idx]
        self.terminals = self.terminals[:, idx]

    def sort(self):
        
        s_idx = sp.lexsort([self.vertices[1, :], self.vertices[0, :]])
        self.reorder(s_idx)

    def from_gene(self, gene):
        
        for transcript_idx in range(len(gene.transcripts)):
            exon_start_end = gene.exons[transcript_idx]
            
            ### only one exon in the transcript
            if exon_start_end.shape[0] == 1:
                exon1_start = exon_start_end[0, 0]
                exon1_end = exon_start_end[0, 1]

                if self.vertices.shape[1] == 0:
                    self.vertices = sp.array([[exon1_start], [exon1_end]], dtype='int')
                    self.edges = sp.array([[0]], dtype='int')
                else:
                    self.vertices = sp.c_[self.vertices, [exon1_start, exon1_end]]
                    self.new_edge()
            ### more than one exon in the transcript
            else:
                for exon_idx in range(exon_start_end.shape[0] - 1):
                    exon1_start = exon_start_end[exon_idx , 0]
                    exon1_end = exon_start_end[exon_idx, 1]
                    exon2_start = exon_start_end[exon_idx + 1, 0]
                    exon2_end = exon_start_end[exon_idx + 1, 1]
          
                    if self.vertices.shape[1] == 0:
                        self.vertices = sp.array([[exon1_start, exon2_start], [exon1_end, exon2_end]], dtype='int')
                        self.edges = sp.array([[0, 1], [1, 0]], dtype='int')
                    else:
                        exon1_idx = -1
                        exon2_idx = -1
                        ### check if current exon already occurred
                        for idx in range(self.vertices.shape[1]):
                            if ((self.vertices[0, idx] == exon1_start) and (self.vertices[1, idx] == exon1_end)):
                                 exon1_idx = idx
                            if ((self.vertices[0, idx] == exon2_start) and (self.vertices[1, idx] == exon2_end)):
                                 exon2_idx = idx

                        ### both exons already occured -> only add an edge
                        if (exon1_idx != -1) and (exon2_idx != -1):
                            self.edges[exon1_idx, exon2_idx] = 1
                            self.edges[exon2_idx, exon1_idx] = 1
                        else:
                            ### 2nd exon occured
                            if ((exon1_idx == -1) and (exon2_idx != -1)):
                                self.vertices = sp.c_[self.vertices, [exon1_start, exon1_end]]
                                self.new_edge()
                                self.edges[exon2_idx, -1] = 1
                                self.edges[-1, exon2_idx] = 1
                            ### 1st exon occured
                            elif ((exon2_idx == -1) and (exon1_idx != -1)):
                                self.vertices = sp.c_[self.vertices, [exon2_start, exon2_end]]
                                self.new_edge()
                                self.edges[exon1_idx, -1] = 1
                                self.edges[-1, exon1_idx] = 1
                            ### no exon occured
                            else:
                                assert((exon1_idx == -1) and (exon2_idx == -1))
                                self.vertices = sp.c_[self.vertices, [exon1_start, exon1_end]]
                                self.vertices = sp.c_[self.vertices, [exon2_start, exon2_end]]
                                self.new_edge()
                                self.new_edge()
                                self.edges[-2, -1] = 1
                                self.edges[-1, -2] = 1

        ### take care of the sorting by exon start
        s_idx = sp.argsort(self.vertices[0, :])
        self.vertices = self.vertices[:, s_idx]
        self.edges = self.edges[s_idx, :][:, s_idx]
        self.terminals = sp.zeros(self.vertices.shape, dtype='int')
        self.terminals[0, sp.where(sp.tril(self.edges).sum(axis=1) == 0)[0]] = 1
        self.terminals[1, sp.where(sp.triu(self.edges).sum(axis=1) == 0)[0]] = 1

    def from_sg_gtf(self, fname):
        """generate a splicing graph from a Splicegrapher GTF"""

        ### prepare containers to collect graph structure from file
        exon_ids = dict()
        exon_coords = []
        parent_pairs = []
        child_pairs = []
        exon_cnt = 0
        start_terminals = []
        end_terminals = []

        for line in open(fname, 'r'):
            sl = line.strip().split('\t')
            if not sl[2] in ['parent', 'child']:
                continue

            ### SplAdder coords are 0-based and half open
            start = int(sl[3]) - 1
            stop = int(sl[4])

            tags = get_tags_gtf(sl[8])

            ### collect node information
            exon_ids[tags['ID']] = exon_cnt
            exon_cnt += 1
            exon_coords.append([start, stop])

            ### collect edge information
            if len(tags['putative_children']) == 0:
                end_terminals.append(c)
            else:
                for c in tags['putative_children'].strip(',').split(','):
                    child_pairs.append(tags['ID'], c)
            if len(tags['putative_parents']) == 0:
                start_terminals.append(c)
            else:
                for c in tags['putative_parents'].strip(',').split(','):
                    parent_pairs.append(tags['ID'], c)

        ### use information to build up graph data structure
        self.vertices = sp.array(exon_coords, dtype='int')
        self.edges = sp.zeros((self.vertices.shape[1], self.vertices.shape[1]), dtype='int')
        self.terminals = sp.zeros_like(self.vertices, dtype='int')
        for pair in child_pairs:
            self.edges[exon_ids[pair[0]], exon_ids[pair[1]]] = 1
            self.edges[exon_ids[pair[1]], exon_ids[pair[0]]] = 1
        for t in start_terminals:
            self.terminals[0, t] = 1
        for t in end_terminals:
            self.terminals[1, t] = 1
            
        

    def from_matfile(self, mat_struct):
        """generates a splicing graph structure from a matfile structure"""

        self.vertices = mat_struct['splicegraph'][0, 0].astype('int')
        self.edges = mat_struct['splicegraph'][0, 1].astype('int')
        self.terminals = mat_struct['splicegraph'][0, 2].astype('int')
        

    def add_intron(self, idx1, flag1, idx2, flag2):
        """adds new introns into splicegraph between idx1 and idx2"""

        ### if flag1, all end terminal exons in idx1 are preserved
        ### if flag2, all start terminal exons in idx2 are preserved

        if idx2.shape[0] > 0:
            adj_mat = sp.triu(self.edges)

            if flag1:
                for i1 in idx1:

                    ### if exon is end-terminal
                    if sp.all(adj_mat[i1, :] == 0):

                        self.vertices = sp.c_[self.vertices, self.vertices[:, i1]]

                        self.new_edge()
                        self.edges[:, -1] = self.edges[:, i1]
                        self.edges[-1, :] = self.edges[i1, :]

                        self.terminals = sp.c_[self.terminals, self.terminals[:, i1]]
            if flag2:
                for i2 in idx2:
                    ### if exon is start-terminal
                    if sp.all(adj_mat[:, i2] == 0):
                        self.vertices = sp.c_[self.vertices, self.vertices[:, i2]]

                        self.new_edge()
                        self.edges[:, -1] = self.edges[:, i2]
                        self.edges[-1, :] = self.edges[i2, :]

                        self.terminals = sp.c_[self.terminals, self.terminals[:, i2]]

        for i1 in idx1:
            for i2 in idx2:
                self.edges[i1, i2] = 1
                self.edges[i2, i1] = 1
        
        self.uniquify()

    def add_cassette_exon(self, new_exon, exons_pre, exons_aft):
        ### exon_pre contains the indices of preceding exons
        ### exon_aft contains the indices of successing exons
        
        self.vertices = sp.c_[self.vertices, new_exon]

        self.new_edge()

        self.edges[exons_pre, -1] = 1
        self.edges[exons_aft, -1] = 1
        self.edges[-1, :] = self.edges[:, -1].T

        self.terminals = sp.c_[self.terminals, sp.zeros((2,), dtype='int')]


    def add_intron_retention(self, idx1, idx2):
        
        adj_mat = sp.triu(self.edges)

        self.vertices = sp.c_[self.vertices, sp.array([self.vertices[0, idx1], self.vertices[1, idx2]], dtype='int')]

        self.new_edge()

        adj_mat = sp.r_[adj_mat, sp.zeros((1, adj_mat.shape[1]), dtype='int')]
        adj_mat = sp.c_[adj_mat, sp.zeros((adj_mat.shape[0], 1), dtype='int')]

        ### check if adjacency matrix is symmetric
        ### otherwise or is not justyfied
        assert(sp.all(sp.all(adj_mat - (self.edges - adj_mat).T == 0)))

        ### AK: under the assumption that our splice graph representation is symmetric
        ### I preserve symmetry by using OR over the adj_mat column and row
        
        self.edges[:, -1] = adj_mat[:, idx1] | adj_mat[idx2, :].T
        self.edges[-1, :] = adj_mat[:, idx1].T | adj_mat[idx2, :]

        self.terminals = sp.c_[self.terminals, sp.array([self.terminals[0, idx1], self.terminals[1, idx2]], dtype='int')]

    def uniquify(self):
        # OUTPUT: splice graph that has been made unique on exons for each gene

        self.sort()
        (s_tmp, s_idx) = sort_rows(self.vertices.T, index=True)
        self.vertices = s_tmp.T
        self.edges = self.edges[s_idx, :][:, s_idx]
        self.terminals = self.terminals[:, s_idx]

        rm_idx = []
        for j in range(1, self.vertices.shape[1]):
            if sp.all(self.vertices[:, j-1] == self.vertices[:, j]):
                self.edges[:, j] = self.edges[:, j-1] | self.edges[:, j]
                self.edges[j, :] = self.edges[j-1, :] | self.edges[j, :]
                rm_idx.append(j - 1)

        keep_idx = sp.where(~sp.in1d(sp.array(sp.arange(self.vertices.shape[1])), rm_idx))[0]
        self.vertices = self.vertices[:, keep_idx]
        self.edges = self.edges[keep_idx, :][:, keep_idx]
        self.terminals = self.terminals[:, keep_idx]

