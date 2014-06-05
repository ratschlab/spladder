function run_merge(CFG)

    merge_all = strcmp(CFG.merge_strategy, 'merge_all');
    merge_all_tag = '';
    if merge_all,
        merge_all_tag = '_merged_bams';
    end;

    chunksize = 50;

    if CFG.do_prune,
        prune_tag = '_pruned';
    else
        prune_tag = '';
    end;

    fn_out = sprintf('%s/spladder/genes_graph_conf%i.%s%s.mat', CFG.out_dirname, CFG.confidence_level, CFG.merge_strategy, prune_tag);
    fn_out_val = sprintf('%s/spladder/genes_graph_conf%i.%s%s.validated.mat', CFG.out_dirname, CFG.confidence_level, CFG.merge_strategy, prune_tag);

    if ~exist(fn_out, 'file'),
        if ~CFG.rproc,
            merge_genes_by_splicegraph(CFG);
        else
            jobinfo = rproc_empty(0) ;
            PAR = struct();
            PAR.CFG = CFG;
            if chunksize > 0,
                merge_list_len = length(CFG.samples);
                if strcmp(CFG.merge_strategy, 'merge_all'),
                    merge_list_len = merge_list_len + 1;
                end;
                for c_idx = 1:chunksize:merge_list_len,
                    fn = sprintf('%s/spladder/genes_graph_conf%i.%s%s_chunk%i_%i.mat', CFG.out_dirname, CFG.confidence_level, CFG.merge_strategy, prune_tag, c_idx, min(merge_list_len, c_idx + chunksize - 1));
                    if exist(fn, 'file'),
                        continue;
                    else
                        fprintf(1, 'submitting chunk %i to %i\n', c_idx, min(merge_list_len, c_idx + chunksize - 1));
                        PAR.chunk_idx = c_idx:min(merge_list_len, c_idx + chunksize - 1);
                        jobinfo(end + 1) = rproc('merge_genes_by_splicegraph', PAR, 50000, CFG.options_rproc, 40*60);
                    end;
                end;
            else
                jobinfo(end + 1) = rproc('merge_genes_by_splicegraph', PAR, 10000, CFG.options_rproc, 40*60);
            end;
            jobinfo = rproc_wait(jobinfo, 30, 1, 1) ;
            %%% merge chunks
            if chunksize > 0,
                PAR.chunksize = chunksize;
                merge_chunks_by_splicegraph(PAR);
            end;
        end;
    else
        fprintf(1, 'File %s already exists!\n', fn_out);
    end;

    %%% generate validated version of splice graph
    if CFG.validate_splicegraphs && ~exist(fn_out_val, 'file'),
        filter_by_edgecount(CFG, fn_out, fn_out_val);
    end;

    if CFG.do_gen_isoforms,
        fn_out = sprintf('%s/spladder/genes_graph_conf%i.%s%s_isoforms.mat', CFG.out_dirname, CFG.confidence_level, CFG.merge_strategy, prune_tag);
        if ~exist(fn_out, 'file');
            if ~CFG.rproc,
                merge_genes_by_isoform(CFG.out_dirname, CFG.confidence_level, merge_all, experiment);
            else
                jobinfo = rproc_empty(0) ;
                jobinfo(end + 1) = rproc('merge_genes_by_isoform', PAR, 10000, CFG.options_rproc, 40*60);
                jobinfo = rproc_wait(jobinfo, 30, 1, 1) ;
            end;
        else
            fprintf(1, 'File %s already exists! \n', fn_out);
        end;
    end;
