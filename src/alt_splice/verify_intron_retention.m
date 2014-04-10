function [verified, info] = verify_intron_retention(event, gene, counts_segments, counts_edges, counts_seg_pos, CFG)
% function [verified, info] = verify_intron_retention(event, gene, counts_segments, counts_edges, counts_seg_pos, CFG)
%

verified = [0 0] ;

% (1) valid, (2) intron_cov, (3) exon1_cov, (4), exon2_cov
% (5) intron_conf, (6) intron_cov_region
info = [1, 0, 0, 0, 0, 0];

ho_offset = 0;
if isfield(CFG, 'is_half_open') && CFG.is_half_open,
    ho_offset = 1;
end;

%%% check validity of exon coordinates (>=0)
if any([event.exon1 event.exon2 event.intron] <= 0),
    info(1) = 0 ;
    return ;
%%% check validity of exon coordinates (start < stop && non-overlapping)
elseif event.exon1(1)>=event.exon1(2) || event.exon2(1)>=event.exon2(2) || ...
        event.intron(1)>=event.intron(2) || event.intron(1)<event.exon1(1) || event.intron(2)>event.exon2(2),
    info(1) = 0 ;
    return ;
end ;

%%% find exons corresponding to event
sg = gene.splicegraph;
segs = gene.segmentgraph;
idx_exon1  = find(sg{1}(1, :) == event.exon1(1) & sg{1}(2, :) == event.exon1(2));
idx_exon2 = find(sg{1}(1, :) == event.exon2(1) & sg{1}(2, :) == event.exon2(2));

%%% find segments corresponding to exons
seg_exon1 = sort(find(segs{2}(idx_exon1, :)));
seg_exon2 = sort(find(segs{2}(idx_exon2, :)));
seg_all = [seg_exon1(1):seg_exon2(end)];
seg_intron = setdiff(seg_all, seg_exon1);
seg_intron = setdiff(seg_intron, seg_exon2);
assert(~isempty(seg_intron));

seg_lens = segs{1}(2, :) - segs{1}(1, :) + 1 - ho_offset;

%%% compute exon coverages as mean of position wise coverage
% exon1_cov
info(3) = sum(counts_segments(seg_exon1) .* seg_lens(seg_exon1)) / sum(seg_lens(seg_exon1));
% exon2_cov
info(4) = sum(counts_segments(seg_exon2) .* seg_lens(seg_exon2)) / sum(seg_lens(seg_exon2));
% intron_cov
info(2) = sum(counts_segments(seg_intron) .* seg_lens(seg_intron)) / sum(seg_lens(seg_intron));
% intron_cov_region
info(6) = sum(counts_seg_pos(seg_intron)) / sum(seg_lens(seg_intron));

%%% check if counts match verification criteria
if info(2) > CFG.intron_retention.min_retention_cov && ...
   info(6) > CFG.intron_retention.min_retention_region && ...
   info(2) >= CFG.intron_retention.min_retention_rel_cov*(info(3) + info(4))/2 
    verified(1) = 1 ;
end ;

%%% check intron confirmation as sum of valid intron scores
%%% intron score is the number of reads confirming this intron
% intron conf
idx = find(counts_edges(:, 1) == sub2ind(size(segs{3}), seg_exon1(end), seg_exon2(1)));
info(5) = counts_edges(idx, 2);

if info(5) >= CFG.intron_retention.min_non_retention_count,
    verified(2) = 1 ;
end;
