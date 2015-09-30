import scipy as sp

class Counts:

    def __init__(self, seg_num):
        self.segments = sp.zeros((seg_num,), dtype='float') 
        self.seg_pos = sp.zeros((seg_num,), dtype='float')
        self.edges = sp.zeros((0, 2), dtype='float')
