import numpy as np

from ..utils import sort_rows

class Event:

    def __init__(self, event_type, chr=None, strand=None):
        
        self.event_type = event_type
        self.chr = chr
        self.strand = strand
        self.strain = ''
        self.exons1 = np.zeros((0, 2), dtype = 'int')
        self.exons2 = np.zeros((0, 2), dtype = 'int')
        self.gene_name = None
        self.transcript_type = None
        self.num_detected = None
        self.id = None
        self.detected = None
        self.annotated = None  ### 0 - both novel; 1 - iso 1 annotated and iso 2 novel; 2 - iso 2 annotated and iso 1 novel; 3 - both annotated 


    def get_len(self):

        return max(self.exons1.max(), self.exons2.max()) - min(self.exons1.min(), self.exons2.min())


    def get_inner_coords(self):
        
        if self.event_type == 'mult_exon_skip':
            return np.sort(np.unique(np.r_[np.sort(self.exons2.ravel())[1:4], np.sort(self.exons2.ravel())[-4:-1]]))
        elif self.event_type == 'mutex_exons':
            return np.sort(np.r_[self.exons1.ravel()[1:4], self.exons2[1, :], self.exons1[2, 0]])
        else:
            return np.sort(np.unique(np.r_[np.sort(self.exons1.ravel())[1:-1], np.sort(self.exons2.ravel())[1:-1]]))
            
        
    def get_coords(self):
        
        if self.event_type != 'mult_exon_skip':
            return np.sort(np.r_[self.exons1.ravel(), self.exons2.ravel()])
        else:
            return np.sort(np.r_[self.exons1.ravel()[:4], self.exons2.ravel()[-4:]])


    def get_exon_coordinate_strings(self):
        
        _iso1 = np.array(['%i-%i' % (_[0], _[1]) for _ in self.exons1])
        _iso2 = np.array(['%i-%i' % (_[0], _[1]) for _ in self.exons2])
        _iso_both = np.unique(np.r_[_iso1, _iso2])
        sidx = np.argsort([int(_.split('-')[0]) for _ in _iso_both])
        _iso_both = _iso_both[sidx]

        _usage = (~(np.in1d(_iso_both, _iso1) & np.in1d(_iso_both, _iso2))).astype('int')

        return ':'.join(_iso_both), ':'.join(_usage.astype('str'))
            
    def get_introns(self):
        
        _introns = np.reshape(self.exons1.ravel()[1:-1], (self.exons1.shape[0] - 1, 2))
        if self.exons2.shape[0] > 1:
            _introns = np.r_[_introns, np.reshape(self.exons2.ravel()[1:-1], (self.exons2.shape[0] - 1, 2))]

        return _introns


    def get_intron_lens(self):

        _introns = self.get_introns()
        return _introns[:, 1] - _introns[:, 0]

    def set_annotation_flag(self, anno_introns):

        self.sort_exons()
        
        ### check annotation status of isoform 1
        self.annotated = 3
        for i in range(self.exons1.shape[0] - 1):
            if not (self.exons1[i, 1], self.exons1[i + 1, 0]) in anno_introns:
                self.annotated -= 1
                break

        ### check annotation status of isoform 2
        if self.exons2.shape[0] < 2:
            return
        for i in range(self.exons2.shape[0] - 1):
            if not (self.exons2[i, 1], self.exons2[i + 1, 0]) in anno_introns:
                self.annotated -= 2
                break

    def sort_exons(self):
        
        self.exons1 = sort_rows(self.exons1)
        self.exons2 = sort_rows(self.exons2)
