""" This module is concerned with functions that are specific to our identity (matlab vs python)"""

import copy
import cPickle
import os

### return identity
def identity():
    return "python"

### return a gene object
def get_gene(gene):
    """This function returns a deep copy of the gene object"""

    return copy.deepcopy(gene)

### load the gene list in the right format
def load_genes(options):
    """This is a helper function to load the gene data from file"""

    if options.validate_sg:
        (genes, events) = cPickle.load(open(os.path.join(options.outdir, 'spladder', 'genes_graph_conf%s.merge_graphs.validated.pickle' % options.confidence), 'r'))
    else:
        (genes, events) = cPickle.load(open(os.path.join(options.outdir, 'spladder', 'genes_graph_conf%s.merge_graphs.pickle' % options.confidence), 'r'))

    return genes


### load a single event from file
def load_event(options, event_info):

    (events, _) = cPickle.load(open(os.path.join(options.outdir, 'merge_graphs_%s_C%s.pickle' % (event_info[0], options.confidence)),'r'))

    return events[int(event_info[1]) - 1]
