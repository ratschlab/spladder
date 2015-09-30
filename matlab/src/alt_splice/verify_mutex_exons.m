function [verified, info] = verify_mutex_exons(event, gene, counts_segments, counts_edges, CFG)
% [verified, info] = verify_mutex_exons(event, gene, counts_segments, counts_edges, CFG)
%

verified = [0 0 0 0] ;

% (1) valid, (2) exon_pre_cov, (3) exon1_cov, (4) exon1_cov, (5) exon_aft_cov, 
% (6) exon_pre_exon1_conf, (7) exon_pre_exon2_conf, (8) exon1_exon_aft_conf, (9) exon2_exon_aft_conf
info = [1, 0, 0, 0, 0, 0, 0, 0, 0];

%%% check validity of exon coordinates (>=0)
if any([event.exon1 event.exon2 event.exon_pre event.exon_aft] <= 0),
    info(1) = 0;
    return ;
%%% check validity of exon coordinates (start < stop && non-overlapping)
elseif event.exon_pre(1) >= event.exon_pre(2) || event.exon_aft(1) >= event.exon_aft(2) || ...
        event.exon1(1) >= event.exon1(2) || event.exon2(1) >= event.exon2(2) || ...
        event.exon_pre(2) >= event.exon1(1) || event.exon1(2) >= event.exon2(1) || event.exon2(2) >= event.exon_aft(1),
    info(1) = 0;
    return ;
end ;

assert(event.exon_pre(2)<event.exon_aft(1)) ;

sg = gene.splicegraph;
segs = gene.segmentgraph;

%%% find exons corresponding to event
idx_exon_pre = find(sg{1}(1, :) == event.exon_pre(1) & sg{1}(2, :) == event.exon_pre(2));
idx_exon_aft = find(sg{1}(1, :) == event.exon_aft(1) & sg{1}(2, :) == event.exon_aft(2));
idx_exon1 = find(sg{1}(1, :) == event.exon1(1) & sg{1}(2, :) == event.exon1(2));
idx_exon2 = find(sg{1}(1, :) == event.exon2(1) & sg{1}(2, :) == event.exon2(2));

%%% find segments corresponding to exons
seg_exon_pre = sort(find(segs{2}(idx_exon_pre, :)));
seg_exon_aft = sort(find(segs{2}(idx_exon_aft, :)));
seg_exon1 = sort(find(segs{2}(idx_exon1, :)));
seg_exon2 = sort(find(segs{2}(idx_exon2, :)));

seg_lens = segs{1}(2, :) - segs{1}(1, :) + 1;

% exon pre cov
info(2) = sum(counts_segments(seg_exon_pre) .* seg_lens(seg_exon_pre)) / sum(seg_lens(seg_exon_pre));
% exon1 cov
info(3) = sum(counts_segments(seg_exon1) .* seg_lens(seg_exon1)) / sum(seg_lens(seg_exon1));
% exon2 cov
info(4) = sum(counts_segments(seg_exon2) .* seg_lens(seg_exon2)) / sum(seg_lens(seg_exon2));
% exon aft cov
info(5) = sum(counts_segments(seg_exon_aft) .* seg_lens(seg_exon_aft)) / sum(seg_lens(seg_exon_aft));

%%% check if coverage of first exon is >= than FACTOR times average of pre and after
if info(3) >= CFG.mutex_exons.min_skip_rel_cov * (info(2) + info(5))/2 
    verified(1) = 1 ;
end ;
if info(4) >= CFG.mutex_exons.min_skip_rel_cov * (info(2) + info(5))/2 
    verified(2) = 1 ;
end ;

%%% check intron confirmation as sum of valid intron scores
%%% intron score is the number of reads confirming this intron
% exon_pre_exon1_conf
idx = find(counts_edges(:, 1) == sub2ind(size(segs{3}), seg_exon_pre(end), seg_exon1(1)));
info(6) = counts_edges(idx, 2);
% exon_pre_exon1_conf
idx = find(counts_edges(:, 1) == sub2ind(size(segs{3}), seg_exon_pre(end), seg_exon2(1)));
info(7) = counts_edges(idx, 2);
% exon1_exon_aft_conf
idx = find(counts_edges(:, 1) == sub2ind(size(segs{3}), seg_exon1(end), seg_exon_aft(1)));
info(8) = counts_edges(idx, 2);
% exon2_exon_aft_conf
idx = find(counts_edges(:, 1) == sub2ind(size(segs{3}), seg_exon2(end), seg_exon_aft(1)));
info(9) = counts_edges(idx, 2);

% set verification flags for intron confirmation
if min(info(6), info(7)) >= CFG.mutex_exons.min_conf_count,
    verified(3) = 1 ;
end ;
if min(info(8), info(9)) >= CFG.mutex_exons.min_conf_count,
    verified(4) = 1 ;
end ;
