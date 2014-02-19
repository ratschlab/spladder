def run_merge(CFG):

    merge_all = (CFG['merge_strategy'] == 'merge_all')
    merge_all_tag = ''
    if merge_all:
        merge_all_tag = '_merged_bams'

    prune_tag = ''
    if CFG['do_prune']:
        prune_tag = '_pruned'

    chunksize = 50

    fn_out = '%s/spladder/genes_graph_conf%i.%s%s.mat' % (CFG['out_dirname'] , CFG['confidence_level'], CFG['merge_strategy'], prune_tag)
    fn_out_val = '%s/spladder/genes_graph_conf%i.%s%s.validated.mat' % (CFG['out_dirname'], CFG['confidence_level'], CFG['merge_strategy'], prune_tag)

    if not os.path.exists(fn_out):
        if not CFG['rproc']:
            merge_genes_by_splicegraph(CFG)
        else:
            jobinfo =[]
            PAR = dict()
            PAR['CFG'] = CFG
            if chunksize > 0:
                merge_list_len = len(CFG['samples'])
                if merge_all:
                    merge_list_len += 1
                for c_idx in range(0, merge_list_len, chunksize):
                    fn = '%s/spladder/genes_graph_conf%i.%s%s_chunk%i_%i.mat' % (CFG['out_dirname'], CFG['confidence_level'], CFG['merge_strategy'], prune_tag, c_idx, min(merge_list_len, c_idx + chunksize - 1))
                    if os.path.exists(fn):
                        continue
                    else:
                        print 'submitting chunk %i to %i' % (c_idx, min(merge_list_len, c_idx + chunksize - 1))
                        PAR['chunk_idx'] = range(c_idx, min(merge_list_len, c_idx + chunksize - 1))
                        jobinfo.append(rproc('merge_genes_by_splicegraph', PAR, 50000, CFG['options_rproc'], 40*60))
            else:
                jobinfo.append(rproc('merge_genes_by_splicegraph', PAR, 10000, CFG['options_rproc'], 40*60))

            jobinfo = rproc_wait(jobinfo, 30, 1, 1)
            ### merge chunks
            if chunksize > 0:
                PAR['chunksize'] = chunksize
                merge_chunks_by_splicegraph(PAR)
    else:
        print 'File %s already exists!' % fn_out

    ### generate validated version of splice graph
    if CFG['validate_splicegraphs'] and not os.path.exists(fn_out_val):
        filter_by_edgecount(CFG, fn_out, fn_out_val)

    if CFG['do_gen_isoforms']:
        fn_out = '%s/spladder/genes_graph_conf%i.%s%s_isoforms.mat' % (CFG['out_dirname'], CFG['confidence_level'], CFG['merge_strategy'], prune_tag)
        if not os.path.exists(fn_out):
            if not CFG['rproc']:
                merge_genes_by_isoform(CFG['out_dirname'], CFG['confidence_level'], merge_all, experiment)
            else:
                jobinfo = rproc('merge_genes_by_isoform', PAR, 10000, CFG['options_rproc'], 40*60)
                jobinfo = rproc_wait(jobinfo, 30, 1, 1)
        else:
            print 'File %s already exists!' % fn_out
