function [verified,info] = verify_alt_prime(event, fn_bam, is_half_open, conf_filter)
% [verified,info]=verify_exon_skip(event, fn_bam, is_half_open, conf_filter)
%

if nargin < 4,
    conf_filter=[] ;
end ;

if nargin < 3,
    is_half_open = 0;
end;

%%% confidence filters for alignments 
if ~isfield(conf_filter, 'intron'),
    conf_filter.intron = 300000;
end 
if ~isfield(conf_filter, 'exon_len'),
    conf_filter.exon_len = 2; 
end
if ~isfield(conf_filter, 'mismatch'),
    conf_filter.mismatch = 1;
end ;
if ~isfield(conf_filter, 'mincount'),
    conf_filter.mincount = 1 ;
end ;

conf_alt_prime = [] ;
if ~isfield(conf_alt_prime, 'min_diff_rel_cov'),
    conf_alt_prime.min_diff_rel_cov = 0.05; 
end;
if ~isfield(conf_alt_prime, 'min_intron_count'),
    conf_alt_prime.min_intron_count = 3 ; 
end

info.exon_diff_cov = 0;
info.exon_const_cov = 0;
info.intron1_conf = 0;
info.intron2_conf = 0;
info.valid = 1 ;

verified = [0 0] ;

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

if is_half_open == 1 && event.strand == '-'
    event.exon_alt1 = event.exon_alt1 + 1;
    event.exon_alt2 = event.exon_alt2 + 1;
    event.exon_const = event.exon_const + 1;
end;

gg.strand = event.strand ;
gg.chr = event.chr ;
gg.chr_num = event.chr_num ;
gg.start = min([event.exon_alt1(1), event.exon_alt2(1), event.exon_const(1)]) ;
gg.stop = max([event.exon_alt1(2), event.exon_alt2(2), event.exon_const(2)]) ;
gg.tracks=[] ;

%%% add RNA-seq evidence to the gene structure
gg = add_count_tracks(gg, fn_bam, conf_filter);

%%% compute exon coverages as mean of position wise coverage
if (event.exon_alt1(1) == event.exon_alt2(1)),
    if event.exon_alt1(2) > event.exon_alt2(2),
        idx_diff = [event.exon_alt2(2) + 1:event.exon_alt1(2)] - gg.start + 1;
        idx_const = [event.exon_alt2(1):event.exon_alt2(2)] - gg.start + 1;
    else
        idx_diff = [event.exon_alt1(2) + 1:event.exon_alt2(2)] - gg.start + 1;
        idx_const = [event.exon_alt1(1):event.exon_alt1(2)] - gg.start + 1;
    end;
elseif (event.exon_alt1(2) == event.exon_alt2(2)),
    if event.exon_alt1(1) < event.exon_alt2(1),
        idx_diff = [event.exon_alt1(1):event.exon_alt2(1) - 1] - gg.start + 1;
        idx_const = [event.exon_alt2(1):event.exon_alt2(2)] - gg.start + 1;
    else
        idx_diff = [event.exon_alt2(1):event.exon_alt1(1) - 1] - gg.start + 1;
        idx_const = [event.exon_alt1(1):event.exon_alt1(2)] - gg.start + 1;
    end;
end;
%if event.exon_alt1(1) > event.exon_const(2),
%    idx_const = [[event.exon_const(1) : event.exon_const(2)] - gg.start + 1, idx_const];
%else
%    idx_const = [idx_const, [event.exon_const(1) : event.exon_const(2)] - gg.start + 1];
%end;

%keyboard;

exon_diff_coverage = mean(sum(gg.tracks(:, idx_diff),1), 2) ;
exon_const_coverage = mean(sum(gg.tracks(:, idx_const),1), 2) ;

info.exon_diff_cov = exon_diff_coverage ;
info.exon_const_cov = exon_const_coverage ;

if exon_diff_coverage >= conf_alt_prime.min_diff_rel_cov * exon_const_coverage
    verified(1) = 1;
end;

%%% check intron confirmations as sum of valid intron scores
%%% intron score is the number of reads confirming this intron
intron_tol = 0 ;
idx = find(abs(gg.introns(1,:) - (event.intron1(1))) <= intron_tol & abs(gg.introns(2,:) - (event.intron1(2))) <= intron_tol) ;
if ~isempty(idx),
    assert(length(idx)>=1) ;
    info.intron1_conf = sum(gg.introns(3,idx)) ;
end ;
idx = find(abs(gg.introns(1,:) - (event.intron2(1))) <= intron_tol & abs(gg.introns(2,:) - (event.intron2(2))) <= intron_tol) ;
if ~isempty(idx),
    assert(length(idx)>=1) ;
    info.intron2_conf = sum(gg.introns(3,idx)) ;
end ;

if min(info.intron1_conf, info.intron2_conf) >= conf_alt_prime.min_intron_count,
    verified(2) = 1 ;
end ;

if is_half_open == 1 && event.strand == '-'
    event.exon_alt1 = event.exon_alt1 - 1;
    event.exon_alt2 = event.exon_alt2 - 1;
    event.exon_const = event.exon_const - 1;
end;
