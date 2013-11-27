function merge_chunks_by_splicegraph(CFG, chunksize) 
%function merge_chunks_by_splicegraph(CFG, chunksize) 
%
%   This script takes several gene structures and merges them into one. 
%   the merge is based on the splicegraphs within the genes struct
    
    %%% if we are running with rproc we only get one parameter struct
    if isfield(CFG, 'CFG'),
        chunksize = CFG.chunksize;
        CFG = CFG.CFG;
    end;

    appended = 0;

    %%% generate merge list
    merge_list = {};
    if CFG.gen_graph_do_prune,
        prune_tag = '_pruned';
    else
        prune_tag = '';
    end;

    merge_list_len = length(CFG.samples);
    if CFG.do_merge_all,
        merge_list_len = merge_list_len + 1;
    end;

    %%% add all single bam file graphs
    for c_idx = 1:chunksize:merge_list_len,
        if CFG.do_merge_all,
            merge_list{end + 1} = sprintf('%s/spladder/genes_graph_conf%i.merged%s_graphs_merged_bams_chunk%i_%i.mat', CFG.out_dirname, CFG.confidence_level, prune_tag, c_idx, min(c_idx + chunksize - 1, merge_list_len));
        else
            merge_list{end + 1} = sprintf('%s/spladder/genes_graph_conf%i.merged%s_graphs_chunk%i_%i.mat', CFG.out_dirname, CFG.confidence_level, prune_tag, c_idx, min(c_idx + chunksize - 1, merge_list_len));
        end;
    end;
    %%% iterate over merge list
    for i = 1:length(merge_list),
        %%% load gene structure from sample i
        fprintf(1, 'Loading %s ...\n', merge_list{i});
        load(merge_list{i});
        fprintf(1, '... done (%i / %i)\n', i, length(merge_list));
        assert(isfield(genes, 'splicegraph'), 'Cannot merge by splicegraph - Field not available!');

        %%% sort genes by name
        name_list = {genes(:).name};
        [name_list, s_idx] = sort(name_list);
        genes = genes(s_idx);

        %%% make sure, that splicegraph is unique
        genes = uniquify_splicegraph(genes);

        %%% jump over first sample - nothig to add yet 
        if (i == 1)
            genes2 = genes;
            for j = 1:length(genes),
                genes2(j).edge_count = genes2(j).splicegraph{2};
            end;
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
            if mod(j, 100) == 0
                fprintf(1, '.');
                if mod(j, 1000) == 0
                    fprintf(1, '%i/%i\n', j, length(genes));
                end;
            end;
            g_idx_ = g_idx;
            while (g_idx <= length(genes2) && strlexcmp(genes2(g_idx).name, genes(j).name) == -1),
                g_idx = g_idx + 1;
            end
            % same gene
            if g_idx <= length(genes2) && strcmp(genes2(g_idx).name, genes(j).name),
                
                [genes(j).splicegraph{1} s_idx] = sortrows(genes(j).splicegraph{1}');
                genes(j).splicegraph{1} = genes(j).splicegraph{1}';

                splice1 = genes(j).splicegraph{2}(s_idx, s_idx);
                splice2 = genes2(g_idx).splicegraph{2};
                edgecnt = genes2(g_idx).edge_count;

                s1_len = size(genes(j).splicegraph{1}, 2);
                s2_len = size(genes2(g_idx).splicegraph{1}, 2);
          
                if s2_len > 10000,
                    fprintf('Do not further merge into gene %i, has more than 10000 vertices!\n', g_idx);
                    %%% still count edges that can be confirmed
                    [~, c_idx, a_idx] = intersect(genes2(g_idx).splicegraph{1}', genes(j).splicegraph{1}',  'rows');
                    if ~isempty(c_idx),
                        genes2(g_idx).edge_count(c_idx, c_idx) = genes2(g_idx).edge_count(c_idx, c_idx) + genes(j).splicegraph{2}(a_idx, a_idx);
                    end;
                else
                    m_graph = [genes(j).splicegraph{1}' ones(s1_len, 1); genes2(g_idx).splicegraph{1}' 2*ones(s2_len, 1)];
                    [tmp, s_idx] = sortrows(m_graph(:,1:3));
                    m_graph = m_graph(s_idx, :);

                    [um_graph, u_l] = unique(m_graph(:,1:2), 'last', 'rows');
                    [tmp, u_f] = unique(m_graph(:,1:2), 'first', 'rows');
                    u_graph = u_l - u_f;
                    
                    f_idx = find(u_graph == 0);

                    if ~isempty(f_idx),
                        splice1_ = zeros(size(u_graph, 1), size(u_graph, 1));
                        splice2_ = zeros(size(u_graph, 1), size(u_graph, 1));
                        edgecnt = zeros(size(u_graph, 1), size(u_graph, 1));
                        idx1_ = find(m_graph(u_f, 3) == 1);
                        idx2_ = find(m_graph(u_l, 3) == 2);
                        splice1_(idx1_, idx1_) = splice1;
                        splice2_(idx2_, idx2_) = splice2;
                        edgecnt(idx2_, idx2_) = genes2(g_idx).edge_count;
                    else
                        splice1_ = splice1;
                        splice2_ = splice2;
                        edgecnt = genes2(g_idx).edge_count;
                    end;
                    if ~all(size(splice1_) == size(splice2_))
                        fprintf(1, 'splice1_ and splice2_ differ in size!\n');
                        keyboard;
                    end;
                    genes2(g_idx).splicegraph{2} = splice1_ | splice2_;
                    genes2(g_idx).splicegraph{1} = um_graph';
                    genes2(g_idx).splicegraph{3} = [[sum(tril(genes2(g_idx).splicegraph{2}), 2) == 0]'; [sum(triu(genes2(g_idx).splicegraph{2}), 2) == 0]'];
                    genes2(g_idx).edge_count = edgecnt + splice1_;
                end;
            %%% we did not find the gene name --> append new gene to genes2
            elseif g_idx > length(genes2) || strlexcmp(genes2(g_idx).name, genes(j).name) > 0,
                g_idx = g_idx_;
                genes2(end+1) = genes(j);
                genes2(end).edge_count = genes(j).splicegraph{2};
                appended = 1;
            end;
        end;
        fprintf(1, '... done\n\n');
        clear genes;
    end;

    genes = genes2;
    clear genes2;

    if merge_all,
        fn = sprintf('%s/spladder/genes_graph_conf%i.merged%s_graphs_merged_bams.mat', CFG.out_dirname, CFG.confidence_level, prune_tag);
    else
        fn = sprintf('%s/spladder/genes_graph_conf%i.merged%s_graphs.mat', CFG.out_dirname, CFG.confidence_level, prune_tag);
    end
    fprintf(1, 'Store genes at: %s\n', fn);
    s = warning('error', 'MATLAB:save:sizeTooBigForMATFile');
    try
        save(fn, 'genes') ;
    catch
        save(fn, 'genes', '-v7.3') ;
    end;
    warning(s);

