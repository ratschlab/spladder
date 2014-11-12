function genes = create_segment_graph(genes, CFG),
    
    for i = 1:length(genes),
        sg = genes(i).splicegraph{1};
        sg(2, :) = sg(2, :) + 1;
        breakpoints = unique(sg(:));
        segments = [];
        for j = 2:length(breakpoints),
            s = sum(sg(1, :) < breakpoints(j));
            e = sum(sg(2, :) < breakpoints(j));
            if s > e,
                segments = [segments, [breakpoints(j-1); breakpoints(j) - 1]];
            end;
        end;
        sg(2, :) = sg(2, :) - 1;
        %%% match nodes to segments
        seg_match = [];
        for j = 1:size(sg, 2),
            tmp = [sg(1, j) <= segments(1, :) & sg(2, j) >= segments(2, :)];
            seg_match = [seg_match; tmp];
        end;
        %%% create edge graph between segments
        seg_edges = zeros(size(segments, 2));
        [k, l] = find(triu(genes(i).splicegraph{2}));
        for m = 1:length(k),
            %%% donor segment
            d = find(seg_match(k(m), :), 1, 'last');
            %%% acceptor segment
            a = find(seg_match(l(m), :), 1, 'first');
            seg_edges(d, a) = 1;
        end;
        genes(i).segmentgraph = {};
        genes(i).segmentgraph{1} = segments;
        genes(i).segmentgraph{2} = seg_match;
        genes(i).segmentgraph{3} = seg_edges;
    end;
