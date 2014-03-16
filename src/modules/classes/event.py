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
        self.detected = None

    def get_len(self, trafo=False):

        if trafo:
            return max(self.exons1_col[-1, -1], self.exons1_col[-1, -1]) - min(self.exons1_col[0, 0], self.exons2_col[0, 0])
        else:
            return max(self.exons1[-1, -1], self.exons1[-1, -1]) - min(self.exons1[0, 0], self.exons2[0, 0])

    def get_inner_coords(self, trafo=False):
        
        if trafo:
            return sp.sort(sp.unique(sp.c_[sp.sort(self.exons1_col.ravel())[1:-1], sp.sort(self.exons2_col.ravel())[1:-1]]))
        else:
            return sp.sort(sp.unique(sp.c_[sp.sort(self.exons1.ravel())[1:-1], sp.sort(self.exons2.ravel())[1:-1]]))
        

    def get_coords(self, trafo=False):
        
        if trafo:
            #return sp.sort(sp.unique(sp.c_[self.exons1_col.ravel(), self.exons2_col.ravel()]))
            return sp.sort(sp.r_[self.exons1_col.ravel(), self.exons2_col.ravel()])
        else:
            #return sp.sort(sp.unique(sp.c_[self.exons1.ravel(), self.exons2.ravel()]))
            return sp.sort(sp.r_[self.exons1.ravel(), self.exons2.ravel()])
        
