function gg = add_count_tracks(gg, fn_bam, conf_filter)

maxval = inf; 
if ~iscell(fn_bam)
    fn = fn_bam ;
    if ~isempty(fn),
        gg = add_reads_from_bam(gg, fn, 'exon_track,intron_list', '', maxval, conf_filter);
    else
        gg.tracks = zeros(1, gg.stop-gg.start) ;
        gg.segment_lists={} ;
        gg.segment_scores={} ;
    end ;
else
    % merge intron lists of several bam files
    segments = zeros(0,3) ;
    % iterate over input files
    for f = 1:length(fn_bam),
        fn = fn_bam{f} ;
        %%% count mincount over all files! (mincount only applied to intron_list)
        conf_filter_ = conf_filter;
        conf_filter_.mincount = 1;
        if ~isempty(fn),
            gg = add_reads_from_bam(gg, fn, 'exon_track,intron_list', '', maxval, conf_filter_);
            segments = [segments; gg.segment_lists{end} gg.segment_scores{end}] ;
        else
            if isempty(gg.tracks),
                gg.tracks = zeros(1, gg.stop-gg.start+1) ;
            else
                gg.tracks(end + 1,:) = 0 ;
            end ;
        end ;
    end ;
    % merge tracks of several bam files
    gg.tracks = sum(gg.tracks, 1); 

    %%% make segments unique
    segments = sortrows(segments) ;
    rm_idx = [] ;
    for i = 1:size(segments,1) - 1,
        if segments(i,1) == segments(i+1,1) && segments(i,2) == segments(i+1,2)
            rm_idx(end+1) = i ;
            segments(i+1,3) = segments(i,3) + segments(i+1,3) ;
        end ;
    end ;
    segments(rm_idx,:) = [] ;

    %%% filter segments by confidence level
    idx = find(segments(:,3) >= conf_filter.mincount) ;
    gg.segment_lists = {segments(idx,1:2)} ;
    gg.segment_scores = {segments(idx,3)} ;
end ;

%%% assemble intron coordinates according to strand
if gg.strand == '+'
    gg.introns = [double([gg.segment_lists{1}(:, 1)'; gg.segment_lists{1}(:, 2)' - 1] + gg.start - 1); gg.segment_scores{1}'];
else
    gg.introns = [double(gg.stop - [gg.segment_lists{1}(:, 2)' - 1; gg.segment_lists{1}(:, 1)'] + 1); gg.segment_scores{1}'];
end ;

if gg.strand=='-',
    gg.tracks = gg.tracks(:,end:-1:1) ;
end ;


