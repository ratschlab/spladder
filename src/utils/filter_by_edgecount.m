function filter_by_edgecount(CFG, fn_genes, fn_out)
% function filter_by_edgecount(CFG, fn_genes, fn_out)

    %%% load gene structure
    fprintf(1, 'Loading merged genes from %s ...\n', fn_genes);
    load(fn_genes);
    fprintf(1, '... done.\n');

    %%% filter splicegraphs by support count over samples
    for i = 1:length(genes),
        u_exons = unique(vertcat(genes(i).exons{:}), 'rows');
        [tmp, tmp, k_idx] = intersect(u_exons, genes(i).splicegraph{1}', 'rows');
        %k_idx = find(sum(genes(i).splicegraph{2}, 2) == 0);
        genes(i).splicegraph{2} = genes(i).edge_count >= CFG.sg_min_edge_count;
        %%% remove all exons that have no incoming or outgoing edges (caused by validation, keep single exon transcripts that occured before)
        rm_idx = setdiff(find(sum(genes(i).splicegraph{2}, 2) == 0), k_idx);
        if ~isempty(rm_idx),
            genes(i).splicegraph{1}(:, rm_idx) = [];
            genes(i).splicegraph{2}(rm_idx, :) = [];
            genes(i).splicegraph{2}(:, rm_idx) = [];
            if length(genes(i).splicegraph) > 2,
                genes(i).splicegraph{3}(:, rm_idx) = [];
            end;
            genes(i).edge_count(:, rm_idx) = [];
            genes(i).edge_count(rm_idx, :) = [];
        end;
    end;

    %genes = rmfield(genes, 'edge_count');
    %%% saving validated splicegraph
    fprintf(1, 'Saving validated genes to: %s\n', fn_out);
    save(fn_out, 'genes');
