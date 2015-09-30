function genes = splice_graph(genes, CFG)
% genes = splice_graph(genes, CFG)

	%%% use annotated exons in genes to build up a new splice graph
	genes = build_splice_graph(genes);

	%%% merge exons if possible / reduce graph
    %%% originially implemented for ESTs, reduces complexity, 
    %%% but removes alternative transcript starts and ends !
	if CFG.do_infer_splice_graph,
		infer_splice_graph;
	end

	%%% sort exons by start position in ascending order
	for ix = 1:length(genes)
		[dummy,exon_order] = sortrows(genes(ix).splicegraph{1}', [1, 2]);
		genes(ix).splicegraph{1} = genes(ix).splicegraph{1}(:,exon_order);
		genes(ix).splicegraph{2} = genes(ix).splicegraph{2}(exon_order,exon_order);
	end

	%%% label alternative and constitutive genes
	genes = label_alt_genes(genes, CFG) ;

    %%% update terminals, start terminals in row 1, end terminals in row 2
    %%% TODO: keep track of terminals in previous steps 
	for i = 1:length(genes)
		adj_mat = triu(genes(i).splicegraph{2}) ;
		genes(i).splicegraph{3} = zeros(2, size(adj_mat,1)) ;
		for j = 1:size(adj_mat, 1),
			genes(i).splicegraph{3}(1, j) = all(adj_mat(:,j) == 0) ;
			genes(i).splicegraph{3}(2, j) = all(adj_mat(j,:) == 0) ;
		end ;
	end ;

	%%% reset gene start and gene stop according to exons and splice graph
	for j = 1:length(genes),
		genes(j).start = inf ;
		genes(j).stop = 0 ;
		for q = 1:size(genes(j).exons),
			genes(j).start = min(genes(j).start, min(genes(j).exons{q}(:))) ;
			genes(j).stop = max(genes(j).stop, max(genes(j).exons{q}(:))) ;
		end ;
		genes(j).start = min(genes(j).start, min(genes(j).splicegraph{1}(1,:))) ;
		genes(j).stop = max(genes(j).stop, max(genes(j).splicegraph{1}(2,:))) ;
		assert(~isinf(genes(j).start) && genes(j).stop ~= 0) ;
	end ;

return

