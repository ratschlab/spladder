function [verified, info] = verify_exon_skip(event, fn_bam, CFG)
% [verified, info] = verify_exon_skip(event, fn_bam, CFG)
%

if ~isfield(CFG, 'is_half_open'),
    is_half_open = 0;
else
    is_half_open = CFG.is_half_open;
end;

verified = [0 0 0 0] ;

info.exon_cov = 0;
info.exon_pre_cov = 0 ;
info.exon_aft_cov = 0 ;
info.exon_pre_exon_conf = 0 ;
info.exon_exon_aft_conf = 0 ;
info.exon_pre_exon_aft_conf = 0 ;
info.valid = 1 ;

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

if is_half_open == 1 && event.strand == '-'
    event.exon_pre = event.exon_pre + 1;
    event.exon_aft = event.exon_aft + 1;
    event.exon = event.exon + 1;
end;

gg.strand = event.strand ;
gg.chr = event.chr ;
gg.chr_num = event.chr_num ;
gg.start = event.exon_pre(1) ;
gg.stop = event.exon_aft(2) ;
gg.tracks = [] ;

assert(event.exon_pre(2)<event.exon_aft(1)) ;

%%% add RNA-seq evidence to the gene structure
gg = add_count_tracks(gg, fn_bam, CFG);

%%% compute exon coverages as mean of position wise coverage
idx = [event.exon_pre(1):event.exon_pre(2)] - gg.start + 1 ;
exon_coverage_pre = mean(gg.tracks(:, idx), 2) ;

idx = [event.exon_aft(1):event.exon_aft(2)] - gg.start + 1 ;
exon_coverage_aft = mean(gg.tracks(:, idx), 2) ;

idx = [event.exon(1):event.exon(2)] - gg.start + 1 ;
exon_coverage = mean(gg.tracks(:, idx), 2) ;

info.exon_pre_cov = exon_coverage_pre ;
info.exon_aft_cov = exon_coverage_aft ;
info.exon_cov = exon_coverage ;

%%% check if coverage of skipped exon is >= than FACTOR times average of pre and after
if exon_coverage >= CFG.exon_skip.min_skip_rel_cov * (exon_coverage_pre + exon_coverage_aft)/2 
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
idx = find(abs(gg.introns(1,:) - (event.exon_pre(2) + 1 - ho_offset)) <= CFG.exon_skip.intron_tolerance & abs(gg.introns(2,:) - (event.exon(1) - 1)) <= CFG.exon_skip.intron_tolerance) ;
if ~isempty(idx),
    assert(length(idx)>=1) ;
    info.exon_pre_exon_conf = sum(gg.introns(3,idx)) ;
end ;
idx = find(abs(gg.introns(1,:) - (event.exon(2) + 1 - ho_offset)) <= CFG.exon_skip.intron_tolerance & abs(gg.introns(2,:) - (event.exon_aft(1) - 1)) <= CFG.exon_skip.intron_tolerance) ;
if ~isempty(idx),
    assert(length(idx)>=1) ;
    info.exon_exon_aft_conf=sum(gg.introns(3,idx)) ;
end ;
idx = find(abs(gg.introns(1,:) - (event.exon_pre(2) + 1 - ho_offset)) <= CFG.exon_skip.intron_tolerance & abs(gg.introns(2,:) - (event.exon_aft(1) - 1)) <= CFG.exon_skip.intron_tolerance) ;
if ~isempty(idx),
    assert(length(idx)>=1) ;
    info.exon_pre_exon_aft_conf = sum(gg.introns(3,idx)) ;
end ;
if info.exon_pre_exon_conf >= CFG.exon_skip.min_non_skip_count,
    verified(2) = 1 ;
end ;
if info.exon_exon_aft_conf >= CFG.exon_skip.min_non_skip_count,
    verified(3) = 1 ;
end ;
if info.exon_pre_exon_aft_conf >= CFG.exon_skip.min_skip_count,
    verified(4) = 1 ;
end ;

if is_half_open == 1 && event.strand == '-'
    event.exon_pre = event.exon_pre - 1;
    event.exon_aft = event.exon_aft - 1;
    event.exon = event.exon - 1;
end;

