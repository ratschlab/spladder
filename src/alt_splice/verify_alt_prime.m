function [verified, info] = verify_alt_prime(event, genes, counts, CFG)
% [verified, info] = verify_exon_skip(event, genes, counts, CFG)
%

info.exon_diff_cov = 0;
info.exon_const_cov = 0;
info.intron1_conf = 0;
info.intron2_conf = 0;
info.valid = 1 ;

verified = [0 0] ;

ho_offset = 0;
if isfield(CFG, 'is_half_open') && CFG.is_half_open,
    ho_offset = 1;
end;

%%% check validity of exon coordinates (>=0)
if any([event.exon_alt1 event.exon_alt2] <= 0),
    info.valid = 0 ;
    return ;
end ;

%%% check validity of intron coordinates (only one side is differing)
if (event.intron1(1) ~= event.intron2(1) && event.intron1(2) ~= event.intron2(2)),
    info.valid = 0 ;
    return ;
end ;

%%% find exons corresponding to event
sg = genes(event.gene_idx).splicegraph;
segs = genes(event.gene_idx).segmentgraph;
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
%seg_const = [seg_const, seg_exon_const]; TODO

info.exon_diff_cov = sum(counts(event.gene_idx).segments(seg_diff) .* (segs{1}(2, seg_diff) - segs{1}(1, seg_diff) + 1 - ho_offset)) / (segs{1}(2, seg_diff(end)) - segs{1}(1, seg_diff(1)) + 1 - ho_offset);
info.exon_const_cov = sum(counts(event.gene_idx).segments(seg_const) .* (segs{1}(2, seg_const) - segs{1}(1, seg_const) + 1 - ho_offset)) / (segs{1}(2, seg_const(end)) - segs{1}(1, seg_const(1)) + 1 - ho_offset);

if info.exon_diff_cov >= CFG.alt_prime.min_diff_rel_cov * info.exon_const_cov,
    verified(1) = 1;
end;

%%% check intron confirmations as sum of valid intron scores
%%% intron score is the number of reads confirming this intron
if seg_exon_alt1(end) < seg_exon_const(1) && seg_exon_alt2(end) < seg_exon_const(1),
    idx = find(counts(event.gene_idx).edges(:, 1) == sub2ind(size(segs{3}), seg_exon_alt1(end), seg_exon_const(1)));
    info.intron1_conf = counts(event.gene_idx).edges(idx, 2);
    idx = find(counts(event.gene_idx).edges(:, 1) == sub2ind(size(segs{3}), seg_exon_alt2(end), seg_exon_const(1)));
    info.intron2_conf = counts(event.gene_idx).edges(idx, 2);
elseif seg_exon_alt1(1) > seg_exon_const(end) && seg_exon_alt2(1) > seg_exon_const(end),
    idx = find(counts(event.gene_idx).edges(:, 1) == sub2ind(size(segs{3}), seg_exon_const(end), seg_exon_alt1(1)));
    info.intron1_conf = counts(event.gene_idx).edges(idx, 2);
    idx = find(counts(event.gene_idx).edges(:, 1) == sub2ind(size(segs{3}), seg_exon_const(end), seg_exon_alt2(1)));
    info.intron2_conf = counts(event.gene_idx).edges(idx, 2);
else
    info.valid = 0;
    return;
end;

if min(info.intron1_conf, info.intron2_conf) >= CFG.alt_prime.min_intron_count,
    verified(2) = 1 ;
end ;
