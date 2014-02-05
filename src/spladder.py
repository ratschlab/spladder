def spladder(ARGS):
    ### return genes

    ### parse parameters from ARGS string
    if type(ARGS) == 'dict':
        CFG = ARGS
        CFG = parse_args('', CFG)
    else:
        CFG = parse_args(ARGS)

    ### add dependencies provided in config section
    if 'paths' in CFG:
        for i in CFG['paths']:
            eval('import %s'% CFG['paths'][i])

    ### load confidence level settings
    if not CFG['no_reset_conf']:
        CFG = set_confidence_level(CFG)

    if not 'spladder_infile' in CFG:
        ### iterate over files, if merge strategy is single
        if CFG['merge_strategy'] in ['single', 'merge_graphs']:
            idxs = range(CFG.samples)
        else:
            idxs = [0]
        
        ### set parallelization
        if CFG.rproc:
            ### TODO init pygrid

        ### create out-directory
        if not os.path.exists(CFG['out_dirname']):
            os.makedirs(CFG['out_dirname'])

        ### create spladder sub-directory
        if not os.path.exists(os.path.join(CFG['out_dirname'], 'spladder')):
            os.makedirs(os.path.join(CFG['out_dirname'], 'spladder'))

        for idx in idxs:
            CFG_ = CFG
            if CFG['merge_strategy'] != 'merge_bams':
                CFG['bam_fnames'] = CFG['bam_fnames'][idx]
                CFG['samples'] = CFG['samples'][idx]
                CFG['out_fname'] = '%s/spladder/genes_graph_conf%i.%s.hdf5' % (CFG['out_dirname'], CFG['confidence_level'], CFG['samples'][0])
            else:
                CFG['out_fname'] = '%s/spladder/genes_graph_conf%i.%s.hdf5' % (CFG['out_dirname'], CFG['confidence_level'], CFG['merge_strategy'])

            ### assemble out filename to check if we are already done
            fn_out = CFG['out_fname']
            if CFG['do_prune']:
                fn_out = re.sub('.mat$', '_pruned.mat', fn_out)
            if CFG['do_gen_isoforms']:
                fn_out = re.sub('.mat$', '_with_isoforms.mat', fn_out)
    
            if os.path.exists(fn_out):
                print >> sys.stdout, 'All result files already exist.'
            else:
                if CFG['rproc']:
                    # TODO: start pygrid
                   # jobinfo(job_nr) = rproc('spladder_core', CFG, 10000, CFG.options_rproc, 40*60) ;
                   # job_nr = job_nr + 1;
                else:
                    spladder_core(CFG)

            CFG = CFG_

        ### collect results after parallelization
        if CFG['rproc']:
            # TODO: pygrid heartbeat
            # jobinfo = rproc_wait(jobinfo, 30, 1, 1) ;

        ### merge parts if necessary
        if CFG['merge_strategy'] == 'merge_graphs':
            run_merge(CFG)

    ### handle alternative splicing part
    if CFG['run_as_analysis']:
        alt_genes_collect(CFG)

        for idx in range(len(CFG['event_types'])):
            alt_genes_analyze(CFG, CFG['event_types'][idx])
