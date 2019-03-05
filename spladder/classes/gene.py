import scipy as sp
import scipy.sparse as spsp
import warnings

if __package__ is None:
    __package__ = 'modules.classes'

from ..utils import *
from .segmentgraph import Segmentgraph
from .splicegraph import Splicegraph

class Gene:
    
    def __init__(self, name=None, start=None, stop=None, chr=None, strand=None, source=None, gene_type=None):
        self.name = name
        self.start = start
        self.stop = stop
        self.exons = []
        self.chr = chr
        if strand in ['+', '-']:
            self.strand = strand
        else:
            warnings.warn('WARNING: strand of gene was provided as %s - automatically set to +') % strand
            self.strand = '+'
        self.transcripts = []
        self.source = source
        self.splicegraph = Splicegraph()
        self.segmentgraph = Segmentgraph()
        self.gene_type = gene_type
        self.is_alt = None
        self.is_alt_spliced = None

    def add_exon(self, exon, idx):
        if idx > (len(self.exons) - 1): 
            self.exons.append(sp.zeros((0, 2), dtype='int'))
        self.exons[idx] = sp.r_[self.exons[idx], sp.expand_dims(exon, axis=0)]

    def label_alt(self):

        # This script labels a gene as 'is_alt_spliced' that
        # has any transcript position that can be intronic or
        # exonic.
        #
        # This script labels each gene 'is_alt', that has a
        # simple alternative transcript start or end. Simple meaning
        # that the remaining transcript structure is the same and no 
        # other exons overlap to the alternative start or end.

        num_exons = self.splicegraph.vertices.shape[1]
      
        vertices = self.splicegraph.vertices
        edges = self.splicegraph.edges
       
        ### no edges in graph --> return
        if num_exons < 2:
            self.is_alt_spliced = 0
            self.is_alt = 0
            return
        
        init_alt_idx = []
        term_alt_idx = []
        for i in range(num_exons):
            ### handle start terminal exons -> no incoming edges
            if not sp.any(edges[:i, i]):
                ### find other exons with same end
                idx = sp.where(~sp.in1d(sp.where(vertices[1, i] == vertices[1, :])[0], i))[0]
                if idx.shape[0] > 0:
                    is_simple = True
                    for j in idx:
                        ### not simple, if different introns follow --> break
                        if not isequal(sp.where(edges[j, j+1:])[0] + j, sp.where(edges[i, i+1:])[0] + i):
                            is_simple = False
                            break
                        ### not simple, if exons previous to j overlap the start of i --> break
                        idx_prev = sp.where(edges[:j, j])[0]
                        for k in idx_prev:
                            if vertices[1, k] > vertices[0, i]:
                                is_simple = False
                                break
                    if not is_simple:
                        continue
                    ### if simple, report alternative initial exon
                    init_alt_idx.append(i)
            ### handle end terminal exons -> no outgoing edges
            if not sp.any(edges[i, i+1:]):
                ### find other exons with the same start
                idx = sp.where(~sp.in1d(sp.where(vertices[0, i] == vertices[0, :])[0], i))[0]
                if idx.shape[0] > 0: 
                    is_simple = True
                    for j in idx:
                        ### not simple, if different introns precede --> break
                        if not isequal(sp.where(edges[:j, j])[0], sp.where(edges[:i, i])[0]):
                            is_simple = False
                            break
                        ### not simple, if exons following to j start before i ends --> break
                        idx_next = sp.where(edges[j, j+1:])[0] + j
                        for k in idx_next:
                            if vertices[0, k] < vertices[1, i]:
                                is_simple = False
                                break
                    if not is_simple:
                        continue
                    ### if simple, report alternative terminal exon
                    term_alt_idx.append(i)
      
        ### further only consider exons that are neither init_alt nor term_alt
        take_idx = sp.where(~sp.in1d(sp.arange(num_exons), [init_alt_idx + term_alt_idx]))[0]
        if take_idx.shape[0] > 0:
            vertices = self.splicegraph.vertices[:, take_idx]
            edges = self.splicegraph.edges[take_idx, :][:, take_idx]
          
            start = vertices[0, :].min()
            stop = vertices[1, :].max()
            exon_loc = sp.zeros((stop - start,))
            
            ### map all edges (introns) to genomic coordinates
            for i in range(edges.shape[0]):
                for j in range(i + 1, edges.shape[0]):
                    if edges[i, j] == 1:
                        cur_edge = sp.array([vertices[1, i], vertices[0, j]]) - start
                        exon_loc[cur_edge[0]:cur_edge[1]] += 1
          
            ### map all vertices (exons) to genomic coordinates 
            for i in range(vertices.shape[1]):
                cur_vertex = vertices[:, i] - start
                exon_loc[cur_vertex[0]:cur_vertex[1]] += 1
        else:
            exon_loc = sp.zeros((1,))
      
        ### if at any position more than one exon or intron -> is_alt_spliced
        if max(exon_loc) > 1:
            self.is_alt_spliced = 1
            self.is_alt = 1
        else:
            self.is_alt_spliced = 0
            ### if not alt_spliced but term_alt or init_alt --> is_alt
            if len(init_alt_idx) == 0 and len(term_alt_idx) == 0:
                self.is_alt = 0 
            else:
                self.is_alt = 1 

    def get_non_alt_seg_ids(self):
        """Returns a list of segment indices that are not alternatively used in the graph.
           
           Currently, the definition of non-alternative means the segment is not spanned
           by any intron. In the narrow sense of non-alternative, one should also exclude
           alternative starts and ends but this could lead to cases where no segment is 
           considered as non-alternative, so we use this approximation for now.
        """

        tmp = sp.ones((self.segmentgraph.seg_edges.shape[0],), dtype='bool')
        for i in range(self.segmentgraph.seg_edges.shape[0] - 1):
            ### get index of last acceptor
            idx = sp.where(self.segmentgraph.seg_edges[i, i + 1:])[0]
            ### mask all segments between current segment and acceptor
            if idx.shape[0] > 0:
                tmp[i + 1:idx[-1] + i + 1] = 0

        return sp.where(tmp)[0]

    def to_sparse(self):
        
        self.splicegraph_edges_shape = self.splicegraph.edges.shape
        self.splicegraph.edges = spsp.coo_matrix(self.splicegraph.edges) 
        self.splicegraph_edges_data = self.splicegraph.edges.data
        self.splicegraph_edges_row = self.splicegraph.edges.row
        self.splicegraph_edges_col = self.splicegraph.edges.col
        self.splicegraph.edges = None

        if hasattr(self, 'edge_count'):
            self.edge_count_shape = self.edge_count.shape
            self.edge_count = spsp.coo_matrix(self.edge_count)
            self.edge_count_data = self.edge_count.data
            self.edge_count_row = self.edge_count.row
            self.edge_count_col = self.edge_count.col
            self.edge_count = None 
        

    def from_sparse(self):
        
        if hasattr(self, 'splicegraph_edges_data'):
            self.splicegraph.edges = spsp.coo_matrix((self.splicegraph_edges_data, (self.splicegraph_edges_row, self.splicegraph_edges_col)), shape=self.splicegraph_edges_shape).toarray()
            del self.splicegraph_edges_data
            del self.splicegraph_edges_shape
            del self.splicegraph_edges_row
            del self.splicegraph_edges_col

        if hasattr(self, 'edge_count_data'):
            self.edge_count = spsp.coo_matrix((self.edge_count_data, (self.edge_count_row, self.edge_count_col)), shape=self.edge_count_shape).toarray()
            del self.edge_count_data
            del self.edge_count_shape
            del self.edge_count_row
            del self.edge_count_col

