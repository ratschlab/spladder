function spladder_core(CFG)
% function spladder_core(CFG)

    %%% add paths
    if isfield(CFG, 'paths'),
        addpath(CFG.paths{:});
    end;

    %%% check if result file exists and start gen graph step if necessary
    if ~exist(CFG.out_fname, 'file'),
        fprintf('Augmenting splice graphs.\n');
        fprintf('=========================\n\n');
        if ~isfield(CFG, 'genes'),
            load(CFG.anno_fname, 'genes');
        else
            genes = CFG.genes ;
        end;
        genes = gen_graphs(genes, CFG) ;

        fprintf('Saving genes to %s\n', CFG.out_fname) ;
        save(CFG.out_fname, 'genes') ;
    else
        fprintf('Augmenting splice graphs already completed.\n%s exists\n ... loading genes from %s ...\n', CFG.out_fname, CFG.out_fname);
        load(CFG.out_fname, 'genes') ;
        fprintf('...done.\n\n');
    end ;

    if CFG.do_prune,
        CFG.out_fname = strrep(CFG.out_fname, '.mat', '_pruned.mat');
    end;
    if CFG.do_prune && ~exist(CFG.out_fname, 'file'),
        genes = uniquify_splicegraph(genes);

        %%% prune graphs
        num_paths_before = count_all_paths(genes);
        genes = prune_graph(genes, bam_fnames);
        num_paths_after = count_all_paths(genes);

        %%% save pruned genes
        fprintf('saving genes to %s\n', CFG.out_fname) ;
        save(CFG.out_fname, 'genes') ;
    else
        if CFG.do_prune,
            fprintf('Pruning of splice graphs already done\n\tloading genes\n');
            load(CFG.out_fname, 'genes') ;
        else
            fprintf('No pruning requested!\n');
        end;
    end;

    if CFG.do_gen_isoforms,
        CFG.out_fame = strrep(CFG.out_fname, '.mat', '_with_isoforms.mat');
    end;
    if CFG.do_gen_isoforms && ~exist(CFG.out_fname, 'file')
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
        if CFG.do_gen_isoforms,
            fprintf('Generating all isoforms already done\n\tloading genes\n');
        else
            fprintf('Generating all isoforms not requested\n');
        end;
    end;
    fprintf('\n');
