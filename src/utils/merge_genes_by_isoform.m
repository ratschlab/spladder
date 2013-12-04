function merge_genes_by_isoform(CFG) 
    %%% This script takes several gene structures and merges them into one. If the number of isoforms of an annotated transcript should 
    %%% exceed max_num_isoforms, max_num_isoforms many are subsampled
    
    appended = 0;
    max_num_isoforms = 10;

    %%% generate merge list
    merge_list = {};
    if CFG.do_prune,
        prune_tag = '_pruned';
    else
        prune_tag = '';
    end;

    %%% add all single bam file isoforms
    for i = 1:length(CFG.samples),
        merge_list{end+1} = sprintf('%s/spladder/genes_graph_conf%i.%s%s.mat', CFG.out_dirname, CFG.confidence_level, CFG.samples{i}, prune_tag);
    end;
    %%% add also isoforms of all bam files combined
    fn_bam_merged = sprintf('%s/spladder/genes_graph_conf%i.merge_bams%s.mat', CFG.out_dirname, CFG.confidence_level, prune_tag);
    if strcmp(CFG.merge_strategy, 'merge_all') && exist(fn_bam_merged, 'file')
        merge_list{end+1} = fn_bam_merged;
    end;

    for i = 1:length(merge_list),
        %%% load gene structure from sample i
        fprintf(1, 'Loading %s ...\n', merge_list{i});
        load(merge_list{i});
        fprintf(1, '... done\n');

        %%% sort
        name_list = {genes(:).name};
        [name_list, s_idx] = sort(name_list);
        genes = genes(s_idx);

        %%% jump over first sample - nothig to add yet 
        if (i == 1)
            genes2 = genes;
            clear genes;
            continue;
        end

        %%% did we append genes in the last round? --> re-sort
        if appended,
            name_list = {genes2(:).name};
            [name_list, s_idx] = sort(name_list);
            genes2 = genes2(s_idx);
            appended = 0;
        end

        %%% iterate over current genes
        g_idx = 1;
        fprintf(1, 'Processing ...\n');
        for j = 1:length(genes);
            if mod(j, 100) == 0,
                fprintf(1, '.');
                if mod(j, 1000) == 0,
                    fprintf(1, '%i/%i\n', j, length(genes));
                end;
            end;
            g_idx_ = g_idx;
            while (strlexcmp(genes2(g_idx).name, genes(j).name) == -1),
                g_idx = g_idx + 1;
            end;
            % same gene
            if strcmp(genes2(g_idx).name, genes(j).name),
                %%% pairwise comparison of all exons
                for k = 1:length(genes(j).exons),
                    broken = 0;
                    for l = 1:length(genes2(g_idx).exons),
                        %%% transcripts are not identical in size --> continue
                        if ~all(size(genes(j).exons{k}) == size(genes2(g_idx).exons{l}))
                            continue ;
                        end
                        %%% found transcript in genes2 that is identical to current transcript --> break
                        if (min(min(genes(j).exons{k} == genes2(g_idx).exons{l})) == 1)
                            broken = 1;
                            break ;
                        end;
                    end;
                    %%% we did not find identical transcript in genes2 --> append transcript
                    if ~broken,
                        genes2(g_idx).exons{end+1} = genes(j).exons{k};
                        genes2(g_idx).transcripts{end+1} = sprintf('%s_%i', genes(j).transcripts{k}, i);
                        genes2(g_idx).start = min(genes2(g_idx).start, min(genes(j).exons{k}(:,1)));
                        genes2(g_idx).stop = max(genes2(g_idx).stop, max(genes(j).exons{k}(:,2)));
                        genes2(g_idx) = splice_graph(genes2(g_idx), 0);
                    end;
                end;
            %%% we did non find the gene name --> append new gene to genes2
            elseif strlexcmp(genes2(g_idx).name, genes(j).name) > 0,
                g_idx = g_idx_;
                genes2(end+1) = genes(j);
                appended = 1;
            end;
        end;
        fprintf(1, '... done\n\n');
        clear genes;
    end;

    genes = genes2;
    clear genes2;

    fn = sprintf('%s/spladder/genes_graph_conf%i.%s%s_merge_isoforms.mat', CFG.out_dirname, CFG.confidence_level, CFG.merge_strategy, prune_tag);
    fprintf(1, 'Store genes at: %s\n', fn);
    save(fn, 'genes') ;

    %%% subsample transcripts if neccessary 
    fprintf(1, 'Subsample genes ...\n');
    for i = 1:length(genes),
        if length(genes(i).transcripts) > max_num_isoforms,
            fprintf(1, 'Subsample for gene %i from %i transcripts.\n', i, length(genes(i).transcripts));
            r_idx = randperm(length(genes(i).transcripts));
            genes(i).transcripts = genes(i).transcripts(r_idx(1:max_num_isoforms));
            genes(i).exons = genes(i).exons(r_idx(1:max_num_isoforms));
            genes(i).start = min(genes(i).start, min(genes(i).exons{1}(:,1)));
            genes(i).stop = max(genes(i).stop, max(genes(i).exons{end}(:,2)));
            genes(i) = splice_graph(genes(i), 0);
        end;
    end;
    fprintf(1, '... done\n\n');

    fn = sprintf('%s/spladder/genes_graph_conf%i.%s%s_merge_isoforms_subsampled.mat', CFG.out_dirname, CFG.confidence_level, CFG.merge_strategy, prune_tag);
    fprintf(1, 'Store subsampled genes at: %s\n', fn);
    save(fn, 'genes') ;
