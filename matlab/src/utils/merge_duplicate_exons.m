function genes = merge_duplicate_exons(genes, CFG)
% genes = merge_duplicate_exons(genes, CFG)

if nargin < 2,
    CFG.verbose = 0;
end;

num_removed = 0 ;
for i = 1:length(genes),
    if CFG.verbose && mod(i, 100) == 0, 
        fprintf('%i\r', i);
    end ;

    %%% check if there are non unique exons
    if size(unique(genes(i).splicegraph{1}', 'rows'), 1) == size(genes(i).splicegraph{1}, 2),
        continue;
    end;

    [dummy, exon_order] = sort(genes(i).splicegraph{1}(1, :), 2, 'ascend');
    genes(i).splicegraph{1} = genes(i).splicegraph{1}(:, exon_order);
    genes(i).splicegraph{2} = genes(i).splicegraph{2}(exon_order, exon_order);
    if length(genes(i).splicegraph) > 2,
        genes(i).splicegraph{3} = genes(i).splicegraph{3}(:, exon_order);
    end ;
    exons = genes(i).splicegraph{1} ;
    admat = genes(i).splicegraph{2} ;
    if length(genes(i).splicegraph) > 2,
        initial = genes(i).splicegraph{3}(1, :) ;
        terminal = genes(i).splicegraph{3}(2, :) ;
    end ;
    remove = zeros(1, size(exons, 2)) ;
    for j = 1:size(exons, 2),
        if remove(j), continue ; end ;
        idx = find(exons(1, j) == exons(1, :) & exons(2, j) == exons(2, :)) ;
        idx = setdiff(idx, j) ;
        remove(idx) = 1 ;
        for k=idx,
            if length(genes(i).splicegraph) > 2,
                initial(j) = initial(j) | initial(k) ;
                terminal(j) = terminal(j) | terminal(k) ;
            end;
            admat(j,:) = admat(j,:) | admat(k,:) ;
            admat(:,j) = admat(:,j) | admat(:,k) ;
        end ;
    end ;
    if sum(remove) == 0
        continue ;
    end ;
    genes(i).splicegraph{2} = admat ;
    if length(genes(i).splicegraph) > 2,
        genes(i).splicegraph{3}(1,:) = initial ;
        genes(i).splicegraph{3}(2,:) = terminal ;
    end ;

    genes(i).splicegraph{1} = genes(i).splicegraph{1}(:,find(~remove));
    genes(i).splicegraph{2} = genes(i).splicegraph{2}(find(~remove),find(~remove));
    if length(genes(i).splicegraph) > 2,
        genes(i).splicegraph{3} = genes(i).splicegraph{3}(:, find(~remove));
    end ;
    num_removed = num_removed + sum(remove) ;
end ;
if CFG.verbose,
    fprintf('\n... removed %i duplicate exons ...\n', num_removed);
end;
