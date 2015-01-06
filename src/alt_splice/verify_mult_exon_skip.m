function [verified, info] = verify_mult_exon_skip(event, gene, counts_segments, counts_edges, CFG)
% [verified, info] = verify_mult_exon_skip(event, gene, counts_segments, counts_edges, CFG) 
%

verified = [0 0 0 0 0] ;
info = [1, 0, 0, 0, 0, 0, 0, 0, 0];
% (1) valid, (2) exon_pre_cov, (3) exons_cov, (4) exon_aft_cov
% (5) exon_pre_exon_conf, (6) exon_exon_aft_conf, (7) exon_pre_exon_aft_conf
% (8) sum_inner_exon_conf, (9) num_inner_exon

%%% check validity of exon coordinates (>=0)
if any([[event.exons(:)]' event.exon_pre event.exon_aft] <= 0),
    info(1) = 0 ;
    return ;
%%% check validity of exon coordinates (start < stop && non-overlapping)
elseif event.exon_pre(1) >= event.exon_pre(2) || event.exon_aft(1) >= event.exon_aft(2) || ...
        any(event.exons(1:2:end) >= event.exons(2:2:end)) || event.exon_pre(2) >= min(event.exons(:)) || max(event.exons(:) >= event.exon_aft(1)),
    info(1) = 0 ;
    return ;
end ;
assert(event.exon_pre(2)<event.exon_aft(1)) ;

%%% find exons corresponding to event
sg = gene.splicegraph;
segs = gene.segmentgraph;
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

seg_lens = segs{1}(2, :) - segs{1}(1, :) + 1;

% exon_pre_cov
info(2) = sum(counts_segments(seg_exon_pre) .* seg_lens(seg_exon_pre)) / sum(seg_lens(seg_exon_pre));
% exon_aft_cov
info(4) = sum(counts_segments(seg_exon_aft) .* seg_lens(seg_exon_aft)) / sum(seg_lens(seg_exon_aft));
% exons_cov
info(3) = sum(counts_segments(seg_exons_u) .* seg_lens(seg_exons_u)) / sum(seg_lens(seg_exons_u));

%%% check if coverage of skipped exon is >= than FACTOR times average of pre and after
if info(3) >= CFG.mult_exon_skip.min_skip_rel_cov * (info(2) + info(4))/2 
    verified(1) = 1 ;
end ;

%%% check intron confirmation as sum of valid intron scores
%%% intron score is the number of reads confirming this intron
% exon_pre_exon_conf
idx = find(counts_edges(:, 1) == sub2ind(size(segs{3}), seg_exon_pre(end), seg_exons{1}(1)));
info(5) = counts_edges(idx, 2);
% exon_exon_aft_conf
idx = find(counts_edges(:, 1) == sub2ind(size(segs{3}), seg_exons{end}(end), seg_exon_aft(1)));
info(6) = counts_edges(idx, 2);
% exon_pre_exon_aft_conf
idx = find(counts_edges(:, 1) == sub2ind(size(segs{3}), seg_exon_pre(end), seg_exon_aft(1)));
info(7) = counts_edges(idx, 2);
for i = 1:(size(seg_exons, 2)-1),
    % sum_inner_exon_conf
    idx = find(counts_edges(:, 1) == sub2ind(size(segs{3}), seg_exons{i}(end), seg_exons{i + 1}(1)));
    info(8) = info(8) + counts_edges(idx, 2);
end;
% num_inner_exon
info(9) = size(event.exons, 2)/2;
if info(5) >= CFG.mult_exon_skip.min_non_skip_count,
    verified(2) = 1 ;
end ;
if info(6) >= CFG.mult_exon_skip.min_non_skip_count,
    verified(3) = 1 ;
end ;
if (info(8) / info(9)) >= CFG.mult_exon_skip.min_non_skip_count,
    verified(4) = 1 ;
end ;
if info(7) >= CFG.mult_exon_skip.min_skip_count,
    verified(5) = 1 ;
end ;
