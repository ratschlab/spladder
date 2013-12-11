function [verified, info] = verify_mult_exon_skip(event, fn_bam, CFG)
% [verified, info] = verify_mult_exon_skip(event, fn_bam, CFG) 
%

if ~isfield(CFG, 'is_half_open'),
    is_half_open = 0;
else
    is_half_open = CFG.is_half_open;
end;

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

if is_half_open == 1 && event.strand == '-'
    event.exon_pre = event.exon_pre + 1;
    event.exon_aft = event.exon_aft + 1;
    event.exons = event.exons + 1;
end;

gg.strand = event.strand ;
gg.chr = event.chr ;
gg.chr_num = event.chr_num ;
gg.start = event.exon_pre(1) ;
gg.stop = event.exon_aft(2) ;
gg.tracks = [] ;

assert(event.exon_pre(2)<event.exon_aft(1)) ;

%%% add RNA-seq evidence to the gene structure
gg = add_count_tracks(gg, fn_bam, CFG.read_filter);

%%% compute exon coverages as mean of position wise coverage
idx = [event.exon_pre(1):event.exon_pre(2)] - gg.start + 1 ;
exon_coverage_pre = mean(sum(gg.tracks(:, idx),1), 2) ;

idx = [event.exon_aft(1):event.exon_aft(2)] - gg.start + 1 ;
exon_coverage_aft = mean(sum(gg.tracks(:, idx),1), 2) ;

idx = [];
for i = 1:2:(size(event.exons, 2) - 1),
    idx = [idx, event.exons(i):event.exons(i+1)];
end;
idx = idx - gg.start + 1;
exons_coverage = mean(sum(gg.tracks(:, idx),1), 2) ;

info.exon_pre_cov = exon_coverage_pre ;
info.exon_aft_cov = exon_coverage_aft ;
info.exons_cov = exons_coverage ;

%%% check if coverage of skipped exon is >= than FACTOR times average of pre and after
if exons_coverage >= CFG.mult_exon_skip.min_skip_rel_cov * (exon_coverage_pre + exon_coverage_aft)/2 
    verified(1) = 1 ;
end ;

%%% set offset for half open intervals
if is_half_open,
    ho_offset = 1;
else
    ho_offset = 0;
end;

%%% check intron confirmation as sum of valid intron scores
%%% intron score is the number of reads confirming this intron
intron_tol = 0 ;
idx = find(abs(gg.introns(1,:) - (event.exon_pre(2) + 1 - ho_offset)) <= intron_tol & abs(gg.introns(2,:) - (event.exons(1) - 1)) <= intron_tol) ;
if ~isempty(idx),
    assert(length(idx)>=1) ;
    info.exon_pre_exon_conf = sum(gg.introns(3,idx)) ;
end ;
idx = find(abs(gg.introns(1,:) - (event.exons(end) + 1 - ho_offset)) <= intron_tol & abs(gg.introns(2,:) - (event.exon_aft(1) - 1)) <= intron_tol) ;
if ~isempty(idx),
    assert(length(idx)>=1) ;
    info.exon_exon_aft_conf = sum(gg.introns(3,idx)) ;
end ;
idx = find(abs(gg.introns(1,:) - (event.exon_pre(2) + 1 - ho_offset)) <= intron_tol & abs(gg.introns(2,:) - (event.exon_aft(1) - 1)) <= intron_tol) ;
if ~isempty(idx),
    assert(length(idx)>=1) ;
    info.exon_pre_exon_aft_conf = sum(gg.introns(3,idx)) ;
end ;
for i = 2:2:(size(event.exons, 2)-2),
    idx = find(abs(gg.introns(1,:) - (event.exons(i) + 1 - ho_offset)) <= intron_tol & abs(gg.introns(2,:) - (event.exons(i+1) - 1)) <= intron_tol) ;
    if ~isempty(idx),
        info.sum_inner_exon_conf = info.sum_inner_exon_conf + sum(gg.introns(3,idx));
    end;
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

if is_half_open == 1 && event.strand == '-'
    event.exon_pre = event.exon_pre - 1;
    event.exon_aft = event.exon_aft - 1;
    event.exons = event.exons - 1;
end;

