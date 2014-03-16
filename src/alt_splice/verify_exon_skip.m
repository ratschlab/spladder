function [verified, info] = verify_exon_skip(event, genes, counts, CFG)
% [verified, info] = verify_exon_skip(event, genes, counts, CFG)
%

verified = [0 0 0 0] ;

info.exon_cov = 0;
info.exon_pre_cov = 0 ;
info.exon_aft_cov = 0 ;
info.exon_pre_exon_conf = 0 ;
info.exon_exon_aft_conf = 0 ;
info.exon_pre_exon_aft_conf = 0 ;
info.valid = 1 ;

ho_offset = 0;
if isfield(CFG, 'is_half_open') && CFG.is_half_open,
    ho_offset = 1;
end;

%%% check validity of exon coordinates (>=0)
if any([event.exon event.exon_pre event.exon_aft] <= 0),
    info.valid = 0 ;
    return ;
%%% check validity of exon coordinates (start < stop && non-overlapping)
elseif event.exon_pre(1) >= event.exon_pre(2) || event.exon_aft(1) >= event.exon_aft(2) || ...
        event.exon(1) >= event.exon(2) || event.exon_pre(2) >= event.exon(1) || event.exon(2) >= event.exon_aft(1),
    info.valid = 0 ;
    return ;
end ;

assert(event.exon_pre(2)<event.exon_aft(1)) ;

sg = genes(event.gene_idx).splicegraph;
segs = genes(event.gene_idx).segmentgraph;

%%% find exons corresponding to event
idx_exon_pre = find(sg{1}(1, :) == event.exon_pre(1) & sg{1}(2, :) == event.exon_pre(2));
idx_exon_aft = find(sg{1}(1, :) == event.exon_aft(1) & sg{1}(2, :) == event.exon_aft(2));
idx_exon = find(sg{1}(1, :) == event.exon(1) & sg{1}(2, :) == event.exon(2));

%%% find segments corresponding to exons
seg_exon_pre = sort(find(segs{2}(idx_exon_pre, :)));
seg_exon_aft = sort(find(segs{2}(idx_exon_aft, :)));
seg_exon = sort(find(segs{2}(idx_exon, :)));

info.exon_pre_cov = sum(counts(event.gene_idx).segments(seg_exon_pre) .* (segs{1}(2, seg_exon_pre) - segs{1}(1, seg_exon_pre) + 1 - ho_offset)) / (segs{1}(2, seg_exon_pre(end)) - segs{1}(1, seg_exon_pre(1)) + 1 - ho_offset);
info.exon_aft_cov = sum(counts(event.gene_idx).segments(seg_exon_aft) .* (segs{1}(2, seg_exon_aft) - segs{1}(1, seg_exon_aft) + 1 - ho_offset)) / (segs{1}(2, seg_exon_aft(end)) - segs{1}(1, seg_exon_aft(1)) + 1 - ho_offset);
info.exon_cov = sum(counts(event.gene_idx).segments(seg_exon) .* (segs{1}(2, seg_exon) - segs{1}(1, seg_exon) + 1 - ho_offset)) / (segs{1}(2, seg_exon(end)) - segs{1}(1, seg_exon(1)) + 1 - ho_offset);

%%% check if coverage of skipped exon is >= than FACTOR times average of pre and after
if info.exon_cov >= CFG.exon_skip.min_skip_rel_cov * (info.exon_pre_cov + info.exon_aft_cov)/2 
    verified(1) = 1 ;
end ;

%%% check intron confirmation as sum of valid intron scores
%%% intron score is the number of reads confirming this intron
idx = find(counts(event.gene_idx).edges(:, 1) == sub2ind(size(segs{3}), seg_exon_pre(end), seg_exon(1)));
info.exon_pre_exon_conf = counts(event.gene_idx).edges(idx, 2);
idx = find(counts(event.gene_idx).edges(:, 1) == sub2ind(size(segs{3}), seg_exon(end), seg_exon_aft(1)));
info.exon_exon_aft_conf = counts(event.gene_idx).edges(idx, 2);
idx = find(counts(event.gene_idx).edges(:, 1) == sub2ind(size(segs{3}), seg_exon_pre(end), seg_exon_aft(1)));
info.exon_pre_exon_aft_conf = counts(event.gene_idx).edges(idx, 2);

if info.exon_pre_exon_conf >= CFG.exon_skip.min_non_skip_count,
    verified(2) = 1 ;
end ;
if info.exon_exon_aft_conf >= CFG.exon_skip.min_non_skip_count,
    verified(3) = 1 ;
end ;
if info.exon_pre_exon_aft_conf >= CFG.exon_skip.min_skip_count,
    verified(4) = 1 ;
end ;

