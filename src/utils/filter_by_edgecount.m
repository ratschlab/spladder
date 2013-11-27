function filter_by_edgecount(fn_genes, fn_out)

    %%% load gene structure
    fprintf(1, 'Loading merged genes from %s ...\n', fn_genes);
    load(fn_genes);
    fprintf(1, '... done.\n');

    %%% filter splicegraphs by support count over samples
    for i = 1:length(genes),
        genes(i).splicegraph{2} = genes(i).edge_count >= 5;
        %%% remove all exons that have no incoming or outgoing edges
        rm_idx = find(sum(genes(i).splicegraph{2}, 2) == 0);
        if ~isempty(rm_idx),
            genes(i).splicegraph{1}(:, rm_idx) = [];
            genes(i).splicegraph{2}(rm_idx, :) = [];
            genes(i).splicegraph{2}(:, rm_idx) = [];
            if length(genes(i).splicegraph) > 2,
                genes(i).splicegraph{3}(:, rm_idx) = [];
            end;
        end;
    end;

    genes = rmfield(genes, 'edge_count');
    %%% saving validated splicegraph
    fprintf(1, 'Saving validated genes to: %s\n', fn_out);
    save(fn_out, 'genes');
