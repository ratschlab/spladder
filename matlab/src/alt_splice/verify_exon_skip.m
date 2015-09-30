function [verified, info] = verify_exon_skip(event, gene, counts_segments, counts_edges, CFG)
% [verified, info] = verify_exon_skip(event, gene, counts_segments, counts_edges, CFG)
%

verified = [0 0 0 0] ;

% (1) valid, (2) exon_cov, (3) exon_pre_cov, (4) exon_aft_cov, 
% (5) exon_pre_exon_conf, (6) exon_exon_aft_conf, (7) exon_pre_exon_aft_conf
info = [1, 0, 0, 0, 0, 0, 0];

%%% check validity of exon coordinates (>=0)
if any([event.exon event.exon_pre event.exon_aft] <= 0),
    info(1) = 0;
    return ;
%%% check validity of exon coordinates (start < stop && non-overlapping)
elseif event.exon_pre(1) >= event.exon_pre(2) || event.exon_aft(1) >= event.exon_aft(2) || ...
        event.exon(1) > event.exon(2) || event.exon_pre(2) >= event.exon(1) || event.exon(2) >= event.exon_aft(1),
    info(1) = 0;
    return ;
end ;

assert(event.exon_pre(2)<event.exon_aft(1)) ;

sg = gene.splicegraph;
segs = gene.segmentgraph;

%%% find exons corresponding to event
idx_exon_pre = find(sg{1}(1, :) == event.exon_pre(1) & sg{1}(2, :) == event.exon_pre(2));
idx_exon_aft = find(sg{1}(1, :) == event.exon_aft(1) & sg{1}(2, :) == event.exon_aft(2));
idx_exon = find(sg{1}(1, :) == event.exon(1) & sg{1}(2, :) == event.exon(2));

%%% find segments corresponding to exons
seg_exon_pre = sort(find(segs{2}(idx_exon_pre, :)));
seg_exon_aft = sort(find(segs{2}(idx_exon_aft, :)));
seg_exon = sort(find(segs{2}(idx_exon, :)));

seg_lens = segs{1}(2, :) - segs{1}(1, :) + 1;

% exon pre cov
info(3) = sum(counts_segments(seg_exon_pre) .* seg_lens(seg_exon_pre)) / sum(seg_lens(seg_exon_pre));
% exon aft cov
info(4) = sum(counts_segments(seg_exon_aft) .* seg_lens(seg_exon_aft)) / sum(seg_lens(seg_exon_aft));
% exon cov
info(2) = sum(counts_segments(seg_exon) .* seg_lens(seg_exon)) / sum(seg_lens(seg_exon));

%%% check if coverage of skipped exon is >= than FACTOR times average of pre and after
if info(2) >= CFG.exon_skip.min_skip_rel_cov * (info(3) + info(4))/2 
    verified(1) = 1 ;
end ;

%%% check intron confirmation as sum of valid intron scores
%%% intron score is the number of reads confirming this intron
idx = find(counts_edges(:, 1) == sub2ind(size(segs{3}), seg_exon_pre(end), seg_exon(1)));
% exon_pre_exon_conf
info(5) = counts_edges(idx, 2);
if info(5) >= CFG.exon_skip.min_non_skip_count,
    verified(2) = 1 ;
end ;
idx = find(counts_edges(:, 1) == sub2ind(size(segs{3}), seg_exon(end), seg_exon_aft(1)));
% exon_exon_aft_conf
info(6) = counts_edges(idx, 2);
if info(6) >= CFG.exon_skip.min_non_skip_count,
    verified(3) = 1 ;
end ;
idx = find(counts_edges(:, 1) == sub2ind(size(segs{3}), seg_exon_pre(end), seg_exon_aft(1)));
% exon_pre_exon_aft_conf
info(7) = counts_edges(idx, 2);
if info(7) >= CFG.exon_skip.min_skip_count,
    verified(4) = 1 ;
end ;
