function merge_genes_by_splicegraph(CFG, chunk_idx) 
%function merge_genes_by_splicegraph(CFG, chunk_idx) 
%
%   This script takes several gene structures and merges them into one. 
%   the merge is based on the splicegraphs within the genes struct
    
    if nargin < 2,
        chunk_idx = [];
    end;

    %%% if we are running with rproc we only get one parameter struct
    if isfield(CFG, 'CFG'),
        if isfield(CFG, 'chunk_idx'),
            chunk_idx = CFG.chunk_idx;
        else
            chunk_idx = [];
        end;
        CFG = CFG.CFG
    elseif nargin < 2,
        chunk_idx = [];
    end;

    %%% generate merge list
    merge_list = {};
    if CFG.do_prune,
        prune_tag = '_pruned';
    else
        prune_tag = '';
    end;

    %%% subset samples in case of chunked computation
    if ~isempty(chunk_idx),
        samples = CFG.samples(chunk_idx);
    else
        samples = CFG.samples;
    end;

    %%% add all single bam file graphs
    for i = 1:length(samples),
        merge_list{end+1} = sprintf('%s/spladder/genes_graph_conf%i.%s%s.mat', CFG.out_dirname, CFG.confidence_level, samples{i}, prune_tag);
    end;
    %%% add also graph of all bam files combined
    if CFG.do_merge_all && exist(sprintf('%s/spladder/genes_graph_conf%i.merge_bams%s.mat', CFG.out_dirname, CFG.confidence_level, prune_tag), 'file')
        merge_list{end+1} = sprintf('%s/spladder/genes_graph_conf%i.merge_bams%s.mat', CFG.out_dirname, CFG.confidence_level, prune_tag);
    end

    %%% iterate over merge list 
    appended = 0;
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

                s1_len = size(genes(j).splicegraph{1}, 2);
                s2_len = size(genes2(g_idx).splicegraph{1}, 2);

                if s2_len > 10000,
                    fprintf('Do not further merge into gene %i, has more than 10000 vertices!\n', g_idx);
                    %%% still count edges that can be confirmed
                    [~, c_idx, a_idx] = intersect(genes2(g_idx).splicegraph{1}', genes(j).splicegraph{1}',  'rows');
                    if ~isempty(c_idx),
                        genes2(g_idx).edge_count(c_idx, c_idx) = genes2(g_idx).edge_count(c_idx, c_idx) + splice1(a_idx, a_idx);
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

    fn_out = sprintf('%s/spladder/genes_graph_conf%i.%s%s.mat', CFG.out_dirname, CFG.confidence_level, CFG.merge_strategy, prune_tag);
    if ~isempty(chunk_idx),
        chunk_tag = sprintf('_chunk%i_%i', chunk_idx(1), chunk_idx(end));
        fn_out = strrep(fn_out, '.mat', sprintf('%s.mat', chunk_tag));
    end;
    
    fprintf(1, 'Store genes at: %s\n', fn_out);
    s = warning('error', 'MATLAB:save:sizeTooBigForMATFile');
    try
        save(fn_out, 'genes') ;
    catch
        save(fn_out, 'genes', '-v7.3');
    end;
    warning(s);
