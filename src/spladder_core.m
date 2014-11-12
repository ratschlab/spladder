function spladder_core(CFG)
% function spladder_core(CFG)

    %%% add paths
    if isfield(CFG, 'paths'),
        addpath(CFG.paths{:});
    end;

    genes_loaded = 0;

    %%% check if result file exists and start gen graph step if necessary
    if ~exist(CFG.out_fname, 'file'),
        fprintf('Augmenting splice graphs.\n');
        fprintf('=========================\n\n');
        if ~isfield(CFG, 'genes'),
            load(CFG.anno_fname, 'genes');
        else
            genes = CFG.genes ;
        end;

        %%% force closed intervals
        if CFG.is_half_open,
            genes = half_open_to_closed(genes);
            CFG.is_half_open = 0;
        end;

        genes = gen_graphs(genes, CFG) ;

        fprintf('Saving genes to %s\n', CFG.out_fname) ;
        save(CFG.out_fname, 'genes') ;

        genes_loaded = 1;
    else
        fprintf('Augmenting splice graphs already completed.\n') ;
    end ;

    %%% prune splice graph if necessary
    if CFG.do_prune,
        load_fn = CFG.out_fname;
        CFG.out_fname = strrep(CFG.out_fname, '.mat', '_pruned.mat');
        if ~exist(CFG.out_fname, 'file'),
            %%% load genes if not present yet
            if genes_loaded == 0,
                load(load_fn);
            end;

            %%% make splice graphs unique
            genes = uniquify_splicegraph(genes);

            %%% prune graphs
            num_paths_before = count_all_paths(genes);
            genes = prune_graph(genes, bam_fnames);
            num_paths_after = count_all_paths(genes);

            %%% save pruned genes
            fprintf('saving genes to %s\n', CFG.out_fname) ;
            save(CFG.out_fname, 'genes') ;

            genes_loaded = 1;
        else
            fprintf('Pruning of splice graphs already done\n');
        end;
    else
        fprintf('No pruning requested!\n');
    end;

    %%% generate isoforms if necessary
    if CFG.do_gen_isoforms,
        load_fn = CFG.out_fname;
        CFG.out_fame = strrep(CFG.out_fname, '.mat', '_with_isoforms.mat');
        if ~exist(CFG.out_fname, 'file'),
            %%% load genes if not present yet
            if genes_loaded == 0,
                load(load_fn);
            end;

            %%% generate isoforms
            fprintf('Generating all isoforms\n') ;
            genes = generate_isoforms(genes, PAR.conf);

            %%% re-constitute splicing graph
            fprintf('\tRe-constituting simplified splice graph from generated transcripts\n') ;
            genes_unsimplified = genes ;
            genes = splice_graph(genes, conf.do_infer_splice_graph);

            %%% save splicing graph with isoforms
            fprintf('\tSaving genes to %s\n', CFG.out_fname) ;
            save(CFG.out_fname, 'genes', 'genes_unsimplified') ;
        else
            fprintf('Generating all isoforms already done\n');
        end;
    else
        fprintf('Generating all isoforms not requested\n');
    end;
    fprintf('\n');
