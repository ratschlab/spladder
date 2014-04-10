function [verified, info] = verify_alt_prime(event, gene, counts_segments, counts_edges, CFG)
% [verified, info] = verify_exon_skip(event, gene, counts_segments, counts_edges, CFG)
%

% (1) valid, (2) exon_diff_cov, (3) exon_const_cov
% (4) intron1_conf, (5) intron2_conf
info = [1, 0, 0, 0, 0];

verified = [0 0] ;

ho_offset = 0;
if isfield(CFG, 'is_half_open') && CFG.is_half_open,
    ho_offset = 1;
end;

%%% check validity of exon coordinates (>=0)
if any([event.exon_alt1 event.exon_alt2] <= 0),
    info(1) = 0 ;
    return ;
end ;

%%% check validity of intron coordinates (only one side is differing)
if (event.intron1(1) ~= event.intron2(1) && event.intron1(2) ~= event.intron2(2)),
    info(1) = 0 ;
    return ;
end ;

%%% find exons corresponding to event
sg = gene.splicegraph;
segs = gene.segmentgraph;
idx_exon_alt1  = find(sg{1}(1, :) == event.exon_alt1(1) & sg{1}(2, :) == event.exon_alt1(2));
idx_exon_alt2  = find(sg{1}(1, :) == event.exon_alt2(1) & sg{1}(2, :) == event.exon_alt2(2));
idx_exon_const = find(sg{1}(1, :) == event.exon_const(1) & sg{1}(2, :) == event.exon_const(2));
if isempty(idx_exon_alt1),
    seg_exon_alt1 = find(segs{1}(1, :) >= event.exon_alt1(1) & segs{1}(2, :) <= event.exon_alt1(2));
else
    seg_exon_alt1 = find(segs{2}(idx_exon_alt1, :));
end;
if isempty(idx_exon_alt2),
    seg_exon_alt2 = find(segs{1}(1, :) >= event.exon_alt2(1) & segs{1}(2, :) <= event.exon_alt2(2));
else
    seg_exon_alt2 = find(segs{2}(idx_exon_alt2, :));
end;
if isempty(idx_exon_const),
    seg_exon_const = find(segs{1}(1, :) >= event.exon_const(1) & segs{1}(2, :) <= event.exon_const(2));
else
    seg_exon_const = find(segs{2}(idx_exon_const, :));
end;
assert(~isempty(seg_exon_alt1));
assert(~isempty(seg_exon_alt2));
assert(~isempty(seg_exon_const));

seg_diff = setdiff(seg_exon_alt1, seg_exon_alt2);
seg_const = intersect(seg_exon_alt2, seg_exon_alt1);
if isempty(seg_diff),
    seg_diff = setdiff(seg_exon_alt2, seg_exon_alt1);
end;
seg_const = [seg_const, seg_exon_const];

seg_lens = segs{1}(2, :) - segs{1}(1, :) + 1 - ho_offset;

% exon_diff_cov
info(2) = sum(counts_segments(seg_diff) .* seg_lens(seg_diff)) / sum(seg_lens(seg_diff));
% exon_const_cov
info(3) = sum(counts_segments(seg_const) .* seg_lens(seg_const)) / sum(seg_lens(seg_const)); 

if info(2) >= CFG.alt_prime.min_diff_rel_cov * info(3),
    verified(1) = 1;
end;

%%% check intron confirmations as sum of valid intron scores
%%% intron score is the number of reads confirming this intron
if seg_exon_alt1(end) < seg_exon_const(1) && seg_exon_alt2(end) < seg_exon_const(1),
    % intron1_conf 
    idx = find(counts_edges(:, 1) == sub2ind(size(segs{3}), seg_exon_alt1(end), seg_exon_const(1)));
    info(4) = counts_edges(idx, 2);
    % intron2_conf 
    idx = find(counts_edges(:, 1) == sub2ind(size(segs{3}), seg_exon_alt2(end), seg_exon_const(1)));
    info(5) = counts_edges(idx, 2);
elseif seg_exon_alt1(1) > seg_exon_const(end) && seg_exon_alt2(1) > seg_exon_const(end),
    % intron1_conf 
    idx = find(counts_edges(:, 1) == sub2ind(size(segs{3}), seg_exon_const(end), seg_exon_alt1(1)));
    info(4) = counts_edges(idx, 2);
    % intron2_conf 
    idx = find(counts_edges(:, 1) == sub2ind(size(segs{3}), seg_exon_const(end), seg_exon_alt2(1)));
    info(5) = counts_edges(idx, 2);
else
    info(1) = 0;
    return;
end;

if min(info(4), info(5)) >= CFG.alt_prime.min_intron_count,
    verified(2) = 1 ;
end ;
