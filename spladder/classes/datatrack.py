import sys
import os
import re
import scipy as sp

class DataTrack:

    TYPES = ['coverage', 'segments', 'splicegraph', 'events', 'transcripts']
    NODATA = ['splicegraph', 'events', 'transcripts']

    def __init__(self, track_type, data):
        self.type = track_type
        if not self.type in DataTrack.TYPES:
            sys.stderr.write('ERROR: Track type %s was given. Allowed track types are: %s\n' % (self.type, ','.join(DataTrack.TYPES)))
            sys.exit(1)
        self.strains = []
        self.bam_fnames = []
        self.group_labels = []
        for data_items in data:
            label = ''
            if ':' in data_items:
                assert len(data_items.split(':')) < 3, 'ERROR: At most one label can be given per data group!\n'
                label, data_items = data_items.split(':')
            self.bam_fnames.append(sp.array(data_items.split(',')))
            self.strains.append(sp.array([re.sub(r'(.[bB][aA][mM]|.[hH][dD][fF]5)$', '', os.path.basename(_)) for _ in self.bam_fnames[-1]]))
            self.group_labels.append(label)
        if self.type in DataTrack.NODATA and len(self.strains) > 0:
            sys.stderr.write('WARNING: Track types %s do not accept any further data items. Given samples will be ignored!\n' % ','.join(DataTrack.NODATA))
            
        
            
