function run_merge(CFG)

    merge_all = strcmp(CFG.merge_strategy, 'merge_all');
    merge_all_tag = '';
    if merge_all,
        merge_all_tag = '_merged_bams';
    end;

    %options.addpaths = {'/cbio/grlab/home/akahles/git/projects/2013/THCA/alternative_splicing', '/cbio/grlab/home/akahles/git/projects/2013/THCA/alternative_splicing/settings', '/cbio/grlab/home/akahles/git/projects/2013/THCA/alternative_splicing/writer', '/cbio/grlab/home/akahles/git/projects/2013/THCA/alternative_splicing/tools'}; 

    chunksize = 50;

    if CFG.do_prune,
        prune_tag = '_pruned';
    else
        prune_tag = '';
    end;

    fn_out = sprintf('%s/spladder/genes_graph_conf%i.merged%s_graphs%s.mat', CFG.out_dirname, CFG.confidence_level, prune_tag, merge_all_tag);
    if ~exist(fn_out, 'file'),
        if ~CFG.rproc,
            merge_genes_by_splicegraph(PAR);
        else
            jobinfo = rproc_empty(0) ;
            PAR.CFG = CFG;
            if chunksize > 0,
                merge_list_len = length(CFG.samples);
                if merge_all,
                    merge_list_len = merge_list_len + 1;
                end;
                for c_idx = 1:chunksize:merge_list_len,
                    if merge_all,
                        fn = sprintf('%s/spladder/genes_graph_conf%i.merged%s_graphs_merged_bams_chunk%i_%i.mat', CFG.out_dirname, CFG.confidence_level, prune_tag, c_idx, min(merge_list_len, c_idx + chunksize - 1));
                    else
                        fn = sprintf('%s/spladder/genes_graph_conf%i.merged%s_graphs_chunk%i_%i.mat', CFG.out_dirname, CFG.confidence_level, prune_tag, c_idx, min(merge_list_len, c_idx + chunksize - 1));
                    end
                    if exist(fn, 'file'),
                        continue;
                    else
                        fprintf(1, 'submitting chunk %i to %i\n', c_idx, min(merge_list_len, c_idx + chunksize - 1));
                        PAR.chunk_idx = c_idx:min(merge_list_len, c_idx + chunksize - 1);
                        jobinfo(end + 1) = rproc('merge_genes_by_splicegraph', PAR, 50000, options, 40*60);
                    end;
                end;
            else
                jobinfo(end + 1) = rproc('merge_genes_by_splicegraph', PAR, 10000, CFG.options_rproc, 40*60);
            end;
            jobinfo = rproc_wait(jobinfo, 30, 1, 1) ;
            if chunksize > 0,
                PAR.chunksize = chunksize;
                merge_chunks_by_splicegraph(PAR);
            end;
        end;
    else
        fprintf(1, 'File %s already exists!\n', fn_out);
    end;
    fn_merged = sprintf('%s/spladder/genes_graph_conf%i.merged%s_graphs%s.mat', CFG.out_dirname, CFG.confidence_level, prune_tag, merge_all_tag);
    fn_val = sprintf('%s/spladder/genes_graph_conf%i.merged%s_graphs%s.validated.mat', CFG.out_dirname, CFG.confidence_level, prune_tag, merge_all_tag);
    if validate_splicegraphs && ~exist(fn_val, 'file'),
        filter_by_edgecount(fn_merged, fn_val);
    end;
    if CFG.do_gen_isoforms,
        fn_out = sprintf('%s/spladder/genes_graph_conf%i.merged%s_isoforms%s.mat', CFG.out_dirname, CFG.confidence_level, prune_tag, merge_all_tag);
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
