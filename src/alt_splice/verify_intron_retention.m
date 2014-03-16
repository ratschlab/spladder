function [verified, info] = verify_intron_retention(event, genes, counts, CFG)
% function [verified, info] = verify_intron_retention(event, genes, counts, CFG)
%

verified = [0 0] ;

info.intron_cov = 0 ;
%info.intron_cov_region = 0 ;
info.exon1_cov = 0 ;
info.exon2_cov = 0 ;
info.intron_conf = 0 ;
info.valid = 1 ;

ho_offset = 0;
if isfield(CFG, 'is_half_open') && CFG.is_half_open,
    ho_offset = 1;
end;

%%% check validity of exon coordinates (>=0)
if any([event.exon1 event.exon2 event.intron] <= 0),
    info.valid = 0 ;
    return ;
%%% check validity of exon coordinates (start < stop && non-overlapping)
elseif event.exon1(1)>=event.exon1(2) || event.exon2(1)>=event.exon2(2) || ...
        event.intron(1)>=event.intron(2) || event.intron(1)<event.exon1(1) || event.intron(2)>event.exon2(2),
    info.valid = 0 ;
    return ;
end ;

%%% find exons corresponding to event
sg = genes(event.gene_idx).splicegraph;
segs = genes(event.gene_idx).segmentgraph;
idx_exon1  = find(sg{1}(1, :) == event.exon1(1) & sg{1}(2, :) == event.exon1(2));
idx_exon2 = find(sg{1}(1, :) == event.exon2(1) & sg{1}(2, :) == event.exon2(2));

%%% find segments corresponding to exons
seg_exon1 = sort(find(segs{2}(idx_exon1, :)));
seg_exon2 = sort(find(segs{2}(idx_exon2, :)));
seg_all = [seg_exon1(1):seg_exon2(end)];
seg_intron = setdiff(seg_all, seg_exon1);
seg_intron = setdiff(seg_intron, seg_exon2);
assert(~isempty(seg_intron));

%%% compute exon coverages as mean of position wise coverage
info.exon1_cov = sum(counts(event.gene_idx).segments(seg_exon1) .* (segs{1}(2, seg_exon1) - segs{1}(1, seg_exon1) + 1 - ho_offset)) / (segs{1}(2, seg_exon1(end)) - segs{1}(1, seg_exon1(1)) + 1 - ho_offset);
info.exon2_cov = sum(counts(event.gene_idx).segments(seg_exon2) .* (segs{1}(2, seg_exon2) - segs{1}(1, seg_exon2) + 1 - ho_offset)) / (segs{1}(2, seg_exon2(end)) - segs{1}(1, seg_exon2(1)) + 1 - ho_offset);
info.intron_cov = sum(counts(event.gene_idx).segments(seg_intron) .* (segs{1}(2, seg_intron) - segs{1}(1, seg_intron) + 1 - ho_offset)) / (segs{1}(2, seg_intron(end)) - segs{1}(1, seg_intron(1)) + 1 - ho_offset);
info.intron_cov_region = sum(counts(event.gene_idx).seg_pos(seg_intron)) / (segs{1}(2, seg_intron(end)) - segs{1}(1, seg_intron(1)) + 1 - ho_offset);

%%% check if counts match verification criteria
if info.intron_cov > CFG.intron_retention.min_retention_cov && ...
   info.intron_cov_region > CFG.intron_retention.min_retention_region && ...
   info.intron_cov >= CFG.intron_retention.min_retention_rel_cov*(info.exon1_cov + info.exon2_cov)/2 
    verified(1) = 1 ;
end ;

%%% check intron confirmation as sum of valid intron scores
%%% intron score is the number of reads confirming this intron
idx = find(counts(event.gene_idx).edges(:, 1) == sub2ind(size(segs{3}), seg_exon1(end), seg_exon2(1)));
info.intron_conf = counts(event.gene_idx).edges(idx, 2);

if info.intron_conf >= CFG.intron_retention.min_non_retention_count,
    verified(2) = 1 ;
end;
