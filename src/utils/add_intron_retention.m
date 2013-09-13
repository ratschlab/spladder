function splicegraph = add_intron_retention(splicegraph, idx1, idx2)

	adj_mat=triu(splicegraph{2}) ;

	splicegraph{1}(:,end+1) = [splicegraph{1}(1, idx1) splicegraph{1}(2, idx2)];
	splicegraph{2}(end+1, end+1) = 0 ;
	adj_mat(end+1, end+1) = 0 ;
	%splicegraph{2}(:,end) = adj_mat(:, idx1) ;
	%splicegraph{2}(end,:) = adj_mat(idx2, :) ;

	%%% check if adjacency matrix is symmetric
	%%% otherwise or is not justyfied
	assert(all(all(adj_mat - (splicegraph{2} - adj_mat)' == 0))) ;

	%%% AK: under the assumption that our splice graph representation is symmetric
	%%% I preserve symmetry by using OR over the adj_mat column and row


	splicegraph{2}(:,end) = or(adj_mat(:, idx1), adj_mat(idx2, :)') ;
	splicegraph{2}(end,:) = or(adj_mat(:, idx1)', adj_mat(idx2, :)) ;

	if length(splicegraph)>2
		splicegraph{3}(1, end+1) = splicegraph{3}(1, idx1);
		splicegraph{3}(2, end) = splicegraph{3}(2, idx2);
	end
