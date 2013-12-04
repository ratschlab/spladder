function [verified, info] = verify_intron_retention(event, fn_bam, CFG)
% function [verified, info] = verify_intron_retention(event, fn_bam, CFG)
%

if ~isfield(CFG, 'is_half_open'),
    is_half_open = 0;
else
    is_half_open = CFG.is_half_open;
end;

verified = [0 0] ;

info.intron_cov = 0 ;
info.intron_cov_region = 0 ;
info.exon1_cov = 0 ;
info.exon2_cov = 0 ;
info.intron_conf = 0 ;
info.valid = 1 ;

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

if is_half_open == 1 && event.strand == '-'
    event.exon1 = event.exon1 + 1;
    event.exon2 = event.exon2 + 1;
    event.intron = event.intron + 1;
end;

gg.strand = event.strand ;
gg.chr = event.chr ;
gg.chr_num = event.chr_num ;
gg.start = event.exon1(1) ;
gg.stop = event.exon2(2) ;
gg.tracks = [] ;
gg.segment_lists = {} ;
gg.segment_scores = {} ;

%%% add RNA-seq evidence to the gene structure
gg = add_count_tracks(gg, fn_bam, CFG.read_filter);

%%% compute exon coverages as mean of position wise coverage
idx = [event.exon1(1):event.exon1(2)] - gg.start + 1 ;
exon_coverage1 = mean(sum(gg.tracks(:, idx),1), 2) ;

idx = [event.exon2(1):event.exon2(2)] - gg.start + 1 ;
exon_coverage2 = mean(sum(gg.tracks(:, idx),1), 2) ;

idx = [event.intron(1):event.intron(2)] - gg.start + 1 ;
idx = idx(idx>0 & idx<=size(gg.tracks,2)) ; % hack
icov = gg.tracks(:,idx) ;

info.intron_cov = mean(icov) ; %mean(icov(1,:) + icov(2,:)) ;
info.intron_cov_region = mean(icov > 0) ; %mean((icov(1,:) + icov(2,:)) > 0) ;
info.exon1_cov = exon_coverage1 ;
info.exon2_cov = exon_coverage2 ;

%%% check if counts match verification criteria
if mean(icov) > CFG.intron_retention.min_retention_cov && ...
        mean(icov > 0) > CFG.intron_retention.min_retention_region && ...
        mean(icov) >= CFG.intron_retention.min_retention_rel_cov*(exon_coverage1+exon_coverage2)/2 
    verified(1) = 1 ;
end ;

%%% check intron confirmation as sum of valid intron scores
%%% intron score is the number of reads confirming this intron
intron_tol = 0;
idx = find(abs(gg.introns(1,:) - event.intron(1)) <= intron_tol & abs(gg.introns(2,:) - event.intron(2)) <= intron_tol) ;
if ~isempty(idx),
    assert(length(idx) >= 1) ;
    info.intron_conf = sum(gg.introns(3,idx)) ;
    if info.intron_conf >= CFG.intron_retention.min_non_retention_count,
        verified(2) = 1 ;
    end;
end ;

if is_half_open == 1 && event.strand == '-'
    event.exon1 = event.exon1 - 1;
    event.exon2 = event.exon2 - 1;
    event.intron = event.intron - 1;
end;
