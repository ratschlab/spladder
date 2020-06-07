import numpy as np

class Counts:

    def __init__(self, seg_num):
        self.segments = np.zeros((seg_num,), dtype='float') 
        self.seg_pos = np.zeros((seg_num,), dtype='float')
        self.edges = np.zeros((0, 2), dtype='float')
