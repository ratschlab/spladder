import scipy as sp

class Gene:
    
    def __init__(self, name=None, start=None, stop=None, chr=None, chr_num=None, strand=None, source=None, gene_type=None):
        self.name = name
        self.start = start
        self.stop = stop
        self.exons = []
        self.chr = chr
        self.chr_num = chr_num
        self.strand = strand
        self.transcripts = []
        self.source = source
        self.splicegraph = None
        self.gene_type=gene_type

    def add_exon(self, exon, idx):
        if idx > (len(self.exons) - 1): 
            self.exons.append(sp.zeros((0, 2), dtype='int'))
        self.exons[idx] = sp.r_[self.exons[idx], sp.expand_dims(exon, axis=0)]
