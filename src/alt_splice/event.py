import scipy as sp

class Event:

    def __init__(self, event_type, chr=None, chr_num=None, strand=None):
        
        self.event_type = event_type
        self.chr = chr
        self.chr_num = chr_num
        self.strand = strand
        self.strain = ''
        self.exons1 = sp.zeros((0, 2), dtype = 'int')
        self.exons2 = sp.zeros((0, 2), dtype = 'int')
        self.exons1_col = sp.zeros((2, 0), dtype = 'int')
        self.exons2_col = sp.zeros((2, 0), dtype = 'int')
        self.p_values = None
        self.gene_name = None
        self.transcript_type = None
        self.num_detected = None
        self.id = None

    def get_len(self, trafo=False):

        if trafo:
            return max(exons1_col[-1, -1], exons1_col[-1, -1]) - min(exons1_col[0, 0], exons2_col[0, 0])
        else:
            return max(exons1[-1, -1], exons1[-1, -1]) - min(exons1[0, 0], exons2[0, 0])

    def get_inner_coords(self, trafo=False):
        
        if trafo:
            tmp = sp.sort(sp.unique(sp.c_[exons1_col.ravel(), exons2_col.ravel()]))
        else:
            tmp = sp.sort(sp.unique(sp.c_[exons1.ravel(), exons2.ravel()]))
        
        if tmp.shape[0] > 2:
            return tmp[1:-1]
        else:
            return []

    def get_coords(self, trafo=False):
        
        if trafo:
            return sp.sort(sp.unique(sp.c_[exons1_col.ravel(), exons2_col.ravel()]))
        else:
            return sp.sort(sp.unique(sp.c_[exons1.ravel(), exons2.ravel()]))
        
