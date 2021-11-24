import sys
import os
import numpy as np
import pickle

from .gen_graphs import gen_graphs 

def spladder_core(bam_fnames, out_fname, options):

    print('Augmenting splice graphs.', file=sys.stdout)
    print('=========================', file=sys.stdout)
    if not hasattr(options, 'genes'):
        genes = pickle.load(open(options.annotation, 'rb'), encoding='latin1')
    else:
        genes = options.genes

    ### mark which introns have been annotated
    for gene in genes:
        gene.populate_annotated_introns()

    ### augment
    genes = gen_graphs(genes, bam_fnames, options)

    print('Saving genes to %s' % (out_fname))
    pickle.dump(genes, open(out_fname, 'wb'), -1)

