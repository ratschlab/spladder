function gg = add_count_tracks(gg, fn_bam, CFG)

if isempty(gg.start),
    %%% TODO!!!!
    gg = build_splice_graph(gg);
    gg.start = min(gg.splicegraph{1}(1, :));
    gg.stop = max(gg.splicegraph{1}(2, :));
end;

maxval = inf; 
if ~iscell(fn_bam)
    fn = fn_bam ;
    if ~isempty(fn),
        gg = add_reads_from_bam(gg, fn, 'exon_track,intron_list', '', maxval, CFG.read_filter, CFG.var_aware, CFG.only_primary);
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
        conf_filter = CFG.read_filter;
        conf_filter.mincount = 1;
        if ~isempty(fn),
            gg = add_reads_from_bam(gg, fn, 'exon_track,intron_list', '', maxval, conf_filter, CFG.var_aware, CFG.only_primary);
            if ~isempty(gg.segment_lists{end}),
                segments = [segments; gg.segment_lists{end} gg.segment_scores{end}] ;
            end;
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
    idx = find(segments(:,3) >= CFG.read_filter.mincount) ;
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


