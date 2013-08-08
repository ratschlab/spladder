function genes = uniquify_splicegraph(genes);
% function genes = uniquify_splicegraph(genes);
% 
% INPUT: genes      gene structure containing all gene information
% OUTPUT: genes     gene structure where splice graph has been made unique on exons for each gene

for i = 1:length(genes)

    [s_tmp, s_idx] = sortrows(genes(i).splicegraph{1}');
    genes(i).splicegraph{1} = s_tmp';
    genes(i).splicegraph{2} = genes(i).splicegraph{2}(s_idx,s_idx);
    genes(i).splicegraph{3} = genes(i).splicegraph{3}(:,s_idx);

    rm_idx = [];
    for j = 2:size(genes(i).splicegraph{1}, 2),
        if isequal(genes(i).splicegraph{1}(:,j-1), genes(i).splicegraph{1}(:,j))
            genes(i).splicegraph{2}(:, j) = genes(i).splicegraph{2}(:, j-1) | genes(i).splicegraph{2}(:, j);
            genes(i).splicegraph{2}(j, :) = genes(i).splicegraph{2}(j-1, :) | genes(i).splicegraph{2}(j, :);
            rm_idx = [rm_idx, j - 1];
        end
    end
    genes(i).splicegraph{1}(:,rm_idx) = [];
    genes(i).splicegraph{2}(:,rm_idx) = [];
    genes(i).splicegraph{2}(rm_idx,:) = [];
    genes(i).splicegraph{3}(:,rm_idx) = [];
end
