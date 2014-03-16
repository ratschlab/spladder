function [verified, info] = verify_mult_exon_skip(event, genes, counts, CFG)
% [verified, info] = verify_mult_exon_skip(event, genes, counts, CFG) 
%

verified = [0 0 0 0 0] ;

info.exons_cov = 0;
info.exon_pre_cov = 0 ;
info.exon_aft_cov = 0 ;
info.exon_pre_exon_conf = 0 ;
info.exon_exon_aft_conf = 0 ;
info.exon_pre_exon_aft_conf = 0 ;
info.sum_inner_exon_conf = 0 ;
info.num_inner_exon = 0 ;
info.valid = 1 ;

ho_offset = 0;
if isfield(CFG, 'is_half_open') && CFG.is_half_open,
    ho_offset = 1;
end;

%%% check validity of exon coordinates (>=0)
if any([[event.exons(:)]' event.exon_pre event.exon_aft] <= 0),
    info.valid = 0 ;
    return ;
%%% check validity of exon coordinates (start < stop && non-overlapping)
elseif event.exon_pre(1) >= event.exon_pre(2) || event.exon_aft(1) >= event.exon_aft(2) || ...
        any(event.exons(1:2:end) >= event.exons(2:2:end)) || event.exon_pre(2) >= min(event.exons(:)) || max(event.exons(:) >= event.exon_aft(1)),
    info.valid = 0 ;
    return ;
end ;
assert(event.exon_pre(2)<event.exon_aft(1)) ;

%%% find exons corresponding to event
sg = genes(event.gene_idx).splicegraph;
segs = genes(event.gene_idx).segmentgraph;
idx_exon_pre  = find(sg{1}(1, :) == event.exon_pre(1) & sg{1}(2, :) == event.exon_pre(2));
idx_exon_aft  = find(sg{1}(1, :) == event.exon_aft(1) & sg{1}(2, :) == event.exon_aft(2));
seg_exons = {};
for i = 1:2:(size(event.exons, 2) - 1),
    tmp = find(sg{1}(1, :) == event.exons(i) & sg{1}(2, :) == event.exons(i+1));
    seg_exons{end + 1} = sort(find(segs{2}(tmp, :)));
end;

%%% find segments corresponding to exons
seg_exon_pre = sort(find(segs{2}(idx_exon_pre, :)));
seg_exon_aft = sort(find(segs{2}(idx_exon_aft, :)));
seg_exons_u = sort(unique([seg_exons{:}]));

info.exon_pre_cov = sum(counts(event.gene_idx).segments(seg_exon_pre) .* (segs{1}(2, seg_exon_pre) - segs{1}(1, seg_exon_pre) + 1 - ho_offset)) / sum(segs{1}(2, seg_exon_pre) - segs{1}(1, seg_exon_pre) + 1 - ho_offset);
info.exon_aft_cov = sum(counts(event.gene_idx).segments(seg_exon_aft) .* (segs{1}(2, seg_exon_aft) - segs{1}(1, seg_exon_aft) + 1 - ho_offset)) /sum (segs{1}(2, seg_exon_aft) - segs{1}(1, seg_exon_aft) + 1 - ho_offset);
info.exons_cov = sum(counts(event.gene_idx).segments(seg_exons_u) .* (segs{1}(2, seg_exons_u) - segs{1}(1, seg_exons_u) + 1 - ho_offset)) / sum(segs{1}(2, seg_exons_u) - segs{1}(1, seg_exons_u) + 1 - ho_offset);

%%% check if coverage of skipped exon is >= than FACTOR times average of pre and after
if info.exons_cov >= CFG.mult_exon_skip.min_skip_rel_cov * (info.exon_pre_cov + info.exon_aft_cov)/2 
    verified(1) = 1 ;
end ;

%%% check intron confirmation as sum of valid intron scores
%%% intron score is the number of reads confirming this intron
idx = find(counts(event.gene_idx).edges(:, 1) == sub2ind(size(segs{3}), seg_exon_pre(end), seg_exons{1}(1)));
info.exon_pre_exon_conf = counts(event.gene_idx).edges(idx, 2);
idx = find(counts(event.gene_idx).edges(:, 1) == sub2ind(size(segs{3}), seg_exons{end}(end), seg_exon_aft(1)));
info.exon_exon_aft_conf = counts(event.gene_idx).edges(idx, 2);
idx = find(counts(event.gene_idx).edges(:, 1) == sub2ind(size(segs{3}), seg_exon_pre(end), seg_exon_aft(1)));
info.exon_pre_exon_aft_conf = counts(event.gene_idx).edges(idx, 2);
for i = 1:(size(seg_exons, 2)-1),
    idx = find(counts(event.gene_idx).edges(:, 1) == sub2ind(size(segs{3}), seg_exons{i}(end), seg_exons{i + 1}(1)));
    info.sum_inner_exon_conf = info.sum_inner_exon_conf + counts(event.gene_idx).edges(idx, 2);
end;
info.num_inner_exon = size(event.exons, 2)/2;
if info.exon_pre_exon_conf >= CFG.mult_exon_skip.min_non_skip_count,
    verified(2) = 1 ;
end ;
if info.exon_exon_aft_conf >= CFG.mult_exon_skip.min_non_skip_count,
    verified(3) = 1 ;
end ;
if (info.sum_inner_exon_conf / info.num_inner_exon) >= CFG.mult_exon_skip.min_non_skip_count,
    verified(4) = 1 ;
end ;
if info.exon_pre_exon_aft_conf >= CFG.mult_exon_skip.min_skip_count,
    verified(5) = 1 ;
end ;
