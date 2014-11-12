function [counts] = count_intron_coverage(genes, fn_bam, CFG, fn_out)

    if nargin < 4,
        fn_out = '';
    end;
    if nargin < 2,
        PAR = genes;
        genes = PAR.genes;
        fn_bam = PAR.fn_bam;
        if isfield(PAR, 'fn_out'),
            fn_out = PAR.fn_out;
        end; 
        CFG = PAR.CFG;
    end;

    counts = struct();
    intron_tol = 0 ;

    if ~iscell(fn_bam),
        fn_bam = {fn_bam};
    end;

    for f = 1:length(fn_bam),
        %%% iterate over all genes and generate counts for
        %%% the segments in the segment graph
        %%% and the splice junctions in the splice graph
        for i = 1:length(genes),
            fprintf(1, '.');
            if mod(i, 50) == 0,
                fprintf('%i\n', i);
            end;
            gg = genes(i);
            if isempty(gg.splicegraph{1}),
                gg = build_splice_graph(gg);
            end;
            if isempty(gg.segmentgraph{1}),
                gg = create_segment_graph(gg, CFG);
            end;
            gg.tracks = [];
            gg.start = min(min(gg.segmentgraph{1}));
            gg.stop = max(max(gg.segmentgraph{1}));
            gg = add_count_tracks(gg, fn_bam{f}, CFG);
            %%% extract mean exon coverage for all segments
            counts(f, i).introns = [];
            seg_edges = triu(genes(i).segmentgraph{3});
            [k, l] = find(seg_edges);
            j = 1;
            for m = 1:length(k),
                %%% identify alternative introns (introns that overlap any segment or other intron)
                if (l(m) - k(m) > 1) || sum(sum(seg_edges(l(m):end, 1:k(m)))) > 1,
                    idx = [genes(i).segmentgraph{1}(2, k(m)) + 1 : genes(i).segmentgraph{1}(2, l(m)) - 1] - gg.start + 1;
                    counts(f, i).introns(j) = mean(sum(gg.tracks(:, idx), 1)); 
                    j = j + 1;
                end;
            end;
        end;
    end;

    if ~isempty(fn_out),
        save(fn_out, 'counts');
    end;
