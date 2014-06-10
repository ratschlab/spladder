function genes = half_open_to_closed(genes)
% genes = half_open_to_closed(genes)

for j = 1:length(genes)
	for k = 1:length(genes(j).exons)
		if genes(j).strand=='+'
			genes(j).exons{k}(:, 2) = genes(j).exons{k}(:, 2) + 1;
		else
			genes(j).exons{k}(:, 1) = genes(j).exons{k}(:, 1) - 1;
		end
	end
    if genes(j).strand=='+'
        if isfield(genes(j), 'splicegraph'),
            genes(j).splicegraph{1}(2, :) = genes(j).splicegraph{1}(2, :) + 1;
        end;
        if isfield(genes(j), 'segmentgraph'),
            genes(j).segmentgraph{1}(2, :) = genes(j).segmentgraph{1}(2, :) + 1;
        end;
    else
        if isfield(genes(j), 'splicegraph'),
            genes(j).splicegraph{1}(1, :) = genes(j).splicegraph{1}(1, :) - 1;
        end;
        if isfield(genes(j), 'segmentgraph'),
            genes(j).segmentgraph{1}(1, :) = genes(j).segmentgraph{1}(1, :) - 1;
        end;
    end;
end
