function splicegraph = add_cassette_exon(splicegraph, new_exon, exons_pre, exons_aft)
    %%% exon_pre contains the indices of preceding exons
    %%% exon_aft contains the indices of successing exons

	splicegraph{1}(:,end+1) = new_exon';
	splicegraph{2}(end+1, end+1) = 0 ;

    splicegraph{2}(exons_pre, end) = 1;
    splicegraph{2}(exons_aft, end) = 1;
    splicegraph{2}(end, :) = splicegraph{2}(:, end)';

    if length(splicegraph) > 2,
        splicegraph{3}(:, end + 1) = [0; 0];
    end;
