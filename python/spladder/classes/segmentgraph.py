import scipy as sp

class Segmentgraph:

    def __init__(self, gene = None):
    
        self.segments = sp.zeros((2, 0), dtype='int')
        self.seg_match = sp.zeros((0, 0), dtype='bool')
        self.seg_edges = sp.zeros((0, 0), dtype='bool')

        if gene is not None:
            self.from_gene(gene)

    def is_empty(self):
        
        return (self.segments.shape[1] == 0) and (self.seg_match.shape[0] == 0) and (self.seg_edges.shape[0] == 0)

    def from_gene(self, gene): 

        sg = gene.splicegraph.vertices
        breakpoints = sp.unique(sg.ravel())
        self.segments = sp.zeros((2, 0), dtype='int')
        for j in range(1, breakpoints.shape[0]):
            s = sp.sum(sg[0, :] < breakpoints[j])
            e = sp.sum(sg[1, :] < breakpoints[j])
            if s > e:
                self.segments = sp.c_[self.segments, [breakpoints[j-1], breakpoints[j]]]

        ### match nodes to segments
        self.seg_match = sp.zeros((0, sg.shape[1]), dtype='bool')
        for j in range(sg.shape[1]):
            tmp = ((sg[0, j] <= self.segments[0, :]) & (sg[1, j] >= self.segments[1, :]))
            if self.seg_match.shape[0] == 0:
                self.seg_match = tmp.copy().reshape((1, tmp.shape[0]))
            else:
                self.seg_match = sp.r_[self.seg_match, tmp.reshape((1, tmp.shape[0]))]

        ### create edge graph between segments
        self.seg_edges = sp.zeros((self.segments.shape[1], self.segments.shape[1]), dtype='bool')
        k, l = sp.where(sp.triu(gene.splicegraph.edges))

        for m in range(k.shape[0]):
            ### donor segment
            d = sp.where(self.seg_match[k[m], :])[0][-1]
            ### acceptor segment
            a = sp.where(self.seg_match[l[m], :])[0][0]
            self.seg_edges[d, a] = True


