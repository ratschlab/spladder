class Counts:

    def __init__(self, seg_num):
        self.segments = sp.zeros((seg_num,), dtype='int') 
        self.seg_pos = sp.zeros((seg_num,), dtype='int')
        self.edges = sp.zeros((0, 2), dtype='int')
