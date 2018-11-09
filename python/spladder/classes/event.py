import scipy as sp

class Event:

    def __init__(self, event_type, chr=None, strand=None):
        
        self.event_type = event_type
        self.chr = chr
        self.strand = strand
        self.strain = ''
        self.exons1 = sp.zeros((0, 2), dtype = 'int')
        self.exons2 = sp.zeros((0, 2), dtype = 'int')
        self.exons1_col = sp.zeros((2, 0), dtype = 'int')
        self.exons2_col = sp.zeros((2, 0), dtype = 'int')
        self.gene_name = None
        self.transcript_type = None
        self.num_detected = None
        self.id = None
        self.detected = None

    def get_len(self, trafo=False):

        if trafo:
            return max(self.exons1_col.max(), self.exons2_col.max()) - min(self.exons1_col.min(), self.exons2_col.min())
        else:
            return max(self.exons1.max(), self.exons2.max()) - min(self.exons1.min(), self.exons2.min())

    def get_inner_coords(self, trafo=False):
        
        if self.event_type == 'mult_exon_skip':
            if trafo:
                return sp.sort(sp.unique(sp.r_[sp.sort(self.exons2_col.ravel())[1:4], sp.sort(self.exons2_col.ravel())[-4:-1]]))
                #return sp.unique(self.exons2_col.ravel())[1:-1]
            else:
                return sp.sort(sp.unique(sp.r_[sp.sort(self.exons2.ravel())[1:4], sp.sort(self.exons2.ravel())[-4:-1]]))
                #return sp.unique(self.exons2.ravel())[1:-1]
        elif self.event_type == 'mutex_exons':
            if trafo:
                return sp.sort(sp.r_[self.exons1_col.ravel()[1:4], self.exons2_col[1, :], self.exons1_col[2, 0]])
            else:
                return sp.sort(sp.r_[self.exons1.ravel()[1:4], self.exons2[1, :], self.exons1[2, 0]])
        else:
            if trafo:
                return sp.sort(sp.unique(sp.r_[sp.sort(self.exons1_col.ravel())[1:-1], sp.sort(self.exons2_col.ravel())[1:-1]]))
            else:
                return sp.sort(sp.unique(sp.r_[sp.sort(self.exons1.ravel())[1:-1], sp.sort(self.exons2.ravel())[1:-1]]))
            
        

    def get_coords(self, trafo=False):
        
        if self.event_type != 'mult_exon_skip':
            if trafo:
                #return sp.sort(sp.unique(sp.c_[self.exons1_col.ravel(), self.exons2_col.ravel()]))
                return sp.sort(sp.r_[self.exons1_col.ravel(), self.exons2_col.ravel()])
            else:
                #return sp.sort(sp.unique(sp.c_[self.exons1.ravel(), self.exons2.ravel()]))
                return sp.sort(sp.r_[self.exons1.ravel(), self.exons2.ravel()])
        else:
            if trafo:
                return sp.sort(sp.r_[self.exons1_col.ravel()[:4], self.exons2_col.ravel()[-4:]])
            else:
                return sp.sort(sp.r_[self.exons1.ravel()[:4], self.exons2.ravel()[-4:]])
            
    def get_introns(self):
        
        _introns = sp.reshape(self.exons1.ravel()[1:-1], (self.exons1.shape[0] - 1, 2))
        if len(self.exons2.shape) > 1:
            _introns = sp.r_[_introns, sp.reshape(self.exons2.ravel()[1:-1], (self.exons2.shape[0] - 1, 2))]

        return _introns


    def get_intron_lens(self):

        _introns = self.get_introns()
        return _introns[:, 1] - _introns[:, 0]

