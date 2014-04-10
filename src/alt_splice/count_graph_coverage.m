function [counts] = count_graph_coverage(genes, fn_bam, CFG, fn_out)

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

    ho_offset = 0;
    if isfield(CFG, 'is_half_open') && CFG.is_half_open,
        ho_offset = 1;
    end;

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
            gg.tracks = [];
            if ho_offset == 1,
                if gg.strand == '-',
                    gg.segmentgraph{1}(1, :) = gg.segmentgraph{1}(1, :) + 1;
                else,
                    gg.segmentgraph{1}(2, :) = gg.segmentgraph{1}(2, :) - 1;
                end;
            end;
            gg.start = min(min(gg.segmentgraph{1}));
            gg.stop = max(max(gg.segmentgraph{1}));
            gg = add_count_tracks(gg, fn_bam{f}, CFG);
            %%% extract mean exon coverage for all segments
            counts(f, i).segments = [];
            counts(f, i).seg_pos = [];
            for j = 1:size(gg.segmentgraph{1}, 2),
                idx = [gg.segmentgraph{1}(1, j) : gg.segmentgraph{1}(2, j)] - gg.start + 1;
                counts(f, i).segments(j) = mean(sum(gg.tracks(:, idx), 1)); 
                counts(f, i).seg_pos(j) = sum(sum(gg.tracks(:, idx), 1) > 0);
            end;
            %%% extract intron counts 
            counts(f, i).edges = [];
            [k, l] = find(gg.segmentgraph{3} == 1);
            for m = 1:length(k),
                idx = find(abs(gg.introns(1,:) - (gg.segmentgraph{1}(2, k(m)) + 1)) <= intron_tol & abs(gg.introns(2,:) - (gg.segmentgraph{1}(1, l(m)) - 1)) <= intron_tol) ;
                if isempty(counts(f, i).edges),
                    if ~isempty(idx),
                        counts(f, i).edges = [sub2ind(size(gg.segmentgraph{3}), k(m), l(m)), sum(gg.introns(3, idx))]; 
                    else,
                        counts(f, i).edges = [sub2ind(size(gg.segmentgraph{3}), k(m), l(m)), 0];
                    end;
                else,
                    if ~isempty(idx),
                        counts(f, i).edges = [counts(f, i).edges; [sub2ind(size(gg.segmentgraph{3}), k(m), l(m)), sum(gg.introns(3, idx))]]; 
                    else,
                        counts(f, i).edges = [counts(f, i).edges; [sub2ind(size(gg.segmentgraph{3}), k(m), l(m)), 0]];
                    end;
                end;
            end;
        end;
    end;

    if ~isempty(fn_out),
        save(fn_out, 'counts');
    end;
