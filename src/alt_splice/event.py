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
        self.introns1 = None
        self.introns2 = None
        self.introns1_col = None
        self.introns2_col = None
        self.p_values = None
        self.gene_name = None
        self.transcript_type = None

