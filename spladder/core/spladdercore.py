import sys
import os
import scipy as sp
import pickle

from .gen_graphs import gen_graphs 

def spladder_core(CFG):

    genes_loaded = False

    ### check if result file exists and start gen graph step if necessary
    if not os.path.exists(CFG['out_fname']):
        print('Augmenting splice graphs.', file=sys.stdout)
        print('=========================', file=sys.stdout)
        if not 'genes' in CFG:
            genes = pickle.load(open(CFG['anno_fname'], 'rb'), encoding='latin1')
        else:
            genes = CFG['genes']

        genes = gen_graphs(genes, CFG)

        print('Saving genes to %s' % (CFG['out_fname']))
        pickle.dump(genes, open(CFG['out_fname'], 'wb'), -1)

        genes_loaded = True
    else:
        print('Augmenting splice graphs already completed.')

    ### prune splice graph if necessary
    if CFG['do_prune']:
        load_fn = CFG['out_fname']
        CFG['out_fname'] = re.sub(r'.pickle$', '_pruned.pickle', CFG['out_fname'])
        if not os.path.exists(CFG['out_fname']):
            ### load genes if not present yet
            if not genes_loaded:
                genes = pickle.load(open(load_fn), 'rb')

            ### make splice graphs unique
            genes = uniquify_splicegraph(genes)

            ### prune graphs
            num_paths_before = count_all_paths(genes)
            genes = prune_graph(genes, bam_fnames)
            num_paths_after = count_all_paths(genes)

            ### save pruned genes
            print('saving genes to %s' % (CFG['out_fname']))
            pickle.dump(genes, open(CFG['out_fname'], 'wb'), -1)

            genes_loaded = True;
        else:
            print('Pruning of splice graphs already done')
    else:
        print('No pruning requested!')

    ### generate isoforms if necessary
    if CFG['do_gen_isoforms']:
        load_fn = CFG['out_fname']
        CFG['out_fame'] = re.sub(r'.pickle$', '_with_isoforms.pickle', CFG['out_fname']);
        if not os.path.exists(CFG['out_fname']):
            ### load genes if not present yet
            if not genes_loaded:
                genes = pickle.load(open(load_fn, 'rb'))

            ### generate isoforms
            print('Generating all isoforms')
            genes = generate_isoforms(genes, PAR['conf'])

            ### re-constitute splicing graph
            print('\tRe-constituting simplified splice graph from generated transcripts')
            genes_unsimplified = genes
            genes = splice_graph(genes, CFG['do_infer_splice_graph'])

            ### save splicing graph with isoforms
            print('\tSaving genes to %s' % CFG['out_fname'])
            pickle.dump((genes, genes_unsimplified), open(CFG['out_fname'], 'wb'), -1)
        else:
            print('Generating all isoforms already done')
    else:
        print('Generating all isoforms not requested')
