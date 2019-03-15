import scipy as sp
import pdb
import sys
import copy

if __package__ is None:
    __package__ = 'modules.core'

from ..reads import *
from ..helpers import *
from ..editgraph import *
from ..merge import *

def gen_graphs(genes, options=None):

    if options is None and isinstance(genes, dict):
        PAR = genes
        genes  = PAR['genes']
        options = PAR['options']

    if options.logfile == '-':
        fd_log = sys.stdout
    else:
        fd_log = open(options.logfile, 'w')

    ### init the stats for inserted elements
    inserted = dict()
    inserted['cassette_exon'] = 0
    inserted['intron_retention'] = 0
    inserted['intron_in_exon'] = 0
    inserted['alt_53_prime'] = 0
    inserted['exon_skip'] = 0
    inserted['gene_merge'] = 0
    inserted['new_terminal_exon'] = 0

    # build splice graph for all genes 
    ##############################################################################%%
    print('Generating splice graph ...', file=fd_log)
    ### merge exons if possible / reduce graph
    ### originially implemented for ESTs, reduces complexity, 
    ### but removes alternative transcript starts and ends !
    if options.infer_sg:
        genes = infer_splice_graph(genes)

    ### sort exons by start position in ascending order
    for ix in range(genes.shape[0]):
        genes[ix].splicegraph.sort()

    ### label alternative and constitutive genes
    for ix in range(genes.shape[0]):
        genes[ix].label_alt()
    if options.verbose:
        print('\nTotal genes:\t\t\t\t\t\t\t%d' % genes.shape[0], file=fd_log)
        print('Total genes with alternative isoforms:\t\t\t\t%d' % sp.sum([x.is_alt for x in genes]), file=fd_log)
        print('Total genes alternatively spliced:\t\t\t\t%d' % sp.sum([x.is_alt_spliced for x in genes]), file=fd_log)
        print('Total constitutively spliced:\t\t\t\t\t%d' % (genes.shape[0] - sp.sum([x.is_alt_spliced for x in genes])), file=fd_log)

    ### update terminals, start terminals in row 1, end terminals in row 2
    for i in range(genes.shape[0]):
        genes[i].splicegraph.update_terminals()

    ### reset gene start and gene stop according to exons and splice graph
    for i in range(genes.shape[0]):
        genes[i].start = min([x.min() for x in genes[i].exons])
        genes[i].stop = max([x.max() for x in genes[i].exons])
    print('...done.\n', file=fd_log)

    ### sort genes by positions
    s_idx = sp.argsort([x.stop for x in genes])
    genes = genes[s_idx]
    s_idx = sp.argsort([x.start for x in genes], kind='mergesort')
    genes = genes[s_idx]
    s_idx = sp.argsort([x.chr for x in genes], kind='mergesort')
    genes = genes[s_idx]

    # append list of introns supported by RNA-seq data to 
    # the genes structure (only in case we want to augment the graph)
    ##############################################################################
    if (options.insert_es or options.insert_ir or options.insert_ni): 
        print('Loading introns from file ...', file=fd_log)
        introns = get_intron_list(genes, options)
        print('...done.\n', file=fd_log)

        ### clean intron list
        ### remove all introns that overlap more than one gene on the same strand
        print('Filtering introns for ambiguity ...', file=fd_log)
        introns = filter_introns(introns, genes, options)
        print('...done.\n', file=fd_log)

        ### check feasibility
        ### TODO when working exclusively with sparse bam, we need to skip this ...
        print('Testing for infeasible genes ...', file=fd_log)
        introns = make_introns_feasible(introns, genes, options)
        print('...done.\n', fd_log)

        for i in range(genes.shape[0]):
            genes[i].introns = introns[i, :]

    ### insert casssette exons
    if options.insert_es:
        print('Inserting cassette exons ...', file=fd_log)
        CFG_ = dict()
        if hasattr(options, 'cassette_exon') and 'read_filter' in options.cassette_exon:
            CFG_['read_filter'] = options.read_filter.copy()
            options.read_filter = options.cassette_exon['read_filter']
        genes, inserted_ = insert_cassette_exons(genes, options)
        inserted['cassette_exon'] = inserted_
        for key in CFG_:
            setattr(options, key, CFG_[key].copy())
        print('\n... inserted %i casette exons ....\n... done.\n' % inserted['cassette_exon'], file=fd_log)

    if options.insert_ir:
        print('Inserting intron retentions ...', file=fd_log)
        CFG_ = dict()
        if 'read_filter' in options.intron_retention:
            CFG_['read_filter'] = options.read_filter.copy()
            options.read_filter = options.intron_retention['read_filter']
        genes, inserted_ = insert_intron_retentions(genes, options)
        inserted['intron_retention'] = inserted_
        for key in CFG_:
            setattr(options, key, CFG_[key].copy())
        print('\n... inserted %i new intron retentions ...\n...done.\n' % inserted['intron_retention'], file=fd_log)

    if options.remove_se:
        print('Removing short exons ...', file=fd_log)
        genes = remove_short_exons(genes, options)
        for i in range(genes.shape[0]):
            if sp.any(genes[i].splicegraph.vertices[:, 1] - genes[i].splicegraph.vertices[:, 0] < options.remove_exons['min_exon_len_remove']):
                print('WARNING: could not remove all short exons', file=sys.stderr)
        print('... done.\n', file=fd_log)


    # test all exons if the reading frame is larger if exon is skipped
    ##############################################################################%%
    #print('find exons to skip to elongate reading frame', file=fd_log)
    #genes = insert_cds_exon_skips(genes, genome_info)
    #genes = splice_graph(genes)


    # sanity checking
    for g in genes:
        assert(all(g.splicegraph.vertices[0, :] <= g.splicegraph.vertices[1, :]))

    ### insert new intron edges
    if options.insert_ni:
        # re-set list of introns supported by RNA-seq data to 
        # the genes structure
        ##############################################################################%%
        for i in range(genes.shape[0]):
            genes[i].introns = introns[i, :]

        print('Inserting new intron edges ...', file=fd_log)
        chrms = sp.array([x.chr for x in genes], dtype='str')
        for chr_idx in sp.unique(chrms):
            genes_before = genes[sp.where(chrms == chr_idx)[0]]
            tmp_genes = copy.deepcopy(genes_before)
            #
            ##############################################################################%%
            if not hasattr(options, 'insert_intron_iterations'):
                options.insert_intron_iterations = 5
            for iter in range(1, options.insert_intron_iterations + 1):
                print('... chr %s - iteration %i/%i\n' % (chr_idx, iter, options.insert_intron_iterations), file=fd_log)

                genes_mod, inserted_ = insert_intron_edges(tmp_genes, options)

                inserted['intron_in_exon'] += inserted_['intron_in_exon']
                inserted['alt_53_prime'] += inserted_['alt_53_prime']
                inserted['exon_skip'] += inserted_['exon_skip']
                inserted['gene_merge'] += inserted_['gene_merge']
                inserted['new_terminal_exon'] += inserted_['new_terminal_exon']

                # in case any exon was inserted that already existed, we merge them into one exon 
                print('... removing duplicate exons ...', file=fd_log)
                genes_mod = merge_duplicate_exons(genes_mod, options)

                # inserted
                if isequal(genes_mod, genes_before):
                    break
                tmp_genes = genes_mod
            chrms = sp.array([x.chr for x in genes], dtype='str')
            genes[sp.where(chrms == chr_idx)[0]] = copy.deepcopy(genes_mod)
        print('... done.\n', file=fd_log)

    print('Re-labeleling new alternative genes ...', file=fd_log)
    for ix in range(genes.shape[0]):
        genes[ix].start = genes[ix].splicegraph.vertices.min()
        genes[ix].stop = genes[ix].splicegraph.vertices.max()
        genes[ix].label_alt()
    print('... done.\n', file=fd_log)

    ### print summary to log file
    print('Inserted:', file=fd_log)
    for fn in inserted:
        print('\t%s:\t%i' % (fn, inserted[fn]), file=fd_log)

    if fd_log != sys.stdout:
        fd_log.close()

    return (genes, inserted)

