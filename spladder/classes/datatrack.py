import sys
import os
import re
import scipy as sp

class DataTrack:

    TYPES = ['coverage', 'segments', 'splicegraph', 'event', 'transcript', 'gene']

    def __init__(self, track_type, data):
        self.type = track_type
        if not self.type in DataTrack.TYPES:
            sys.stderr.write('ERROR: Track type %s was given. Allowed track types are: %s\n' % (self.type, ','.join(DataTrack.TYPES)))
            sys.exit(1)
        self.strains = []
        self.bam_fnames = []
        self.group_labels = []
        self.event_info = []
        ### handle events
        if track_type == 'event':
            self.event_info.extend(data)
        else:
            for data_items in data:
                ### handle any other type of data
                label = ''
                if ':' in data_items:
                    assert len(data_items.split(':')) < 3, 'ERROR: At most one label can be given per data group!\n'
                    label, data_items = data_items.split(':')
                if data_items.endswith('.txt'):
                    self.bam_fnames.extend([str(x) for x in sp.atleast_1d(sp.loadtxt(data_items, dtype='str'))])
                else:
                    self.bam_fnames.append(sp.array(data_items.split(',')))
                if track_type == 'coverage':
                    for group in self.bam_fnames:
                        for fname in group:
                            if not os.path.isfile(fname):
                                print('ERROR: Input file %s can not be found\n\n' % fname, file=sys.stderr)
                                sys.exit(2)
                    
                self.strains.append(sp.array([re.sub(r'(.[bB][aA][mM]|.[hH][dD][fF]5)$', '', os.path.basename(_)) for _ in self.bam_fnames[-1]]))
                self.group_labels.append(label)
                
