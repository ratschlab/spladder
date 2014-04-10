function count_graph_coverage_wrapper(fname_in, fname_out, CFG)

    load(fname_in);
    
    if ~isfield(genes, 'segmentgraph'),
        genes = create_segment_graph(genes, CFG);
        secsave(fname_in, genes, 'genes');
    end;

    counts.gene_ids_segs = {};
    counts.gene_ids_edges = {};
    counts.segments = {};
    counts.seg_pos = {};
    counts.edges = {};
    if ~CFG.rproc,
        for s_idx = 1:length(CFG.strains),
            fprintf('%i/%i\r', j, length(CFG.strains));
            counts_tmp(s_idx, :) = count_graph_coverage(genes, CFG.bam_fnames(s_idx), CFG);
        end ;
        for c = 1:size(counts_tmp, 2),
            counts.segments{c} = [vertcat(counts_tmp(:, c).segments)'];
            counts.seg_pos{c} = [vertcat(counts_tmp(:, c).seg_pos)'];
            counts.gene_ids_segs{c} = [repmat(c, size(counts_tmp(1, c).seg_pos, 2), 1)];
            tmp = [counts_tmp(:, c).edges];
            if ~isempty(tmp),
                counts.edges{c} = [tmp(:, 1), tmp(:, 2:2:end)];
                counts.gene_ids_edges{c} = [repmat(c, size(tmp, 1), 1)];
            end;
        end;
    else
        chunksize = 2;
        jobinfo = rproc_empty(0) ;
        PAR = struct();
        PAR.CFG = CFG;
        for c_idx = 1:chunksize:length(genes),
            cc_idx = min(length(genes), c_idx + chunksize - 1);
            fn = strrep(fname_out, '.mat', sprintf('.chunk_%i_%i.mat', c_idx, cc_idx));
            if exist(fn, 'file'),
                fprintf(1, '%s already exists\n', fn);
                continue;
            else
                fprintf(1, 'submitting chunk %i to %i\n', c_idx, cc_idx);
                PAR.genes = genes(c_idx:cc_idx);
                PAR.fn_bam = CFG.bam_fnames;
                PAR.fn_out = fn;
                PAR.CFG = CFG;
                jobinfo(end + 1) = rproc('count_graph_coverage', PAR, 8000, CFG.options_rproc, 48*60);
                %count_graph_coverage(PAR);
            end;
        end;

        jobinfo = rproc_wait(jobinfo, 30, 1, 1) ;

        %%% merge results
        for c_idx = 1:chunksize:length(genes),
            fprintf('%i/%i\r', c_idx, length(genes));
            cc_idx = min(length(genes), c_idx + chunksize - 1);
            fn = strrep(fname_out, '.mat', sprintf('.chunk_%i_%i.mat', c_idx, cc_idx));
            if ~exist(fn, 'file'),
                error('Not all chunks in counting graph coverage completed!');
            else
                L = load(fn);
                for c = 0:chunksize-1,
                    if c_idx + c <= length(genes),
                        counts.segments{c_idx + c} = [vertcat(L.counts(:, c+1).segments)'];
                        counts.seg_pos{c_idx + c} = [vertcat(L.counts(:, c+1).seg_pos)'];
                        counts.gene_ids_segs{c_idx + c} = [repmat(c_idx + c, size(L.counts(1, c+1).seg_pos, 2), 1)];
                        tmp = [L.counts(:, c+1).edges];
                        if ~isempty(tmp),
                            counts.edges{c_idx + c} = [tmp(:, 1), tmp(:, 2:2:end)];
                            counts.gene_ids_edges{c_idx + c} = [repmat(c_idx + c, size(tmp, 1), 1)];
                        end;
                    end;
                end;
            end;
        end;
    end;

    %%% collect data
    counts.gene_ids_segs = vertcat(counts.gene_ids_segs{:});
    counts.gene_ids_edges = vertcat(counts.gene_ids_edges{:});
    counts.segments = vertcat(counts.segments{:});
    counts.seg_pos = vertcat(counts.seg_pos{:});
    counts.edges = vertcat(counts.edges{:});
    counts.edge_idx = counts.edges(:, 1);
    counts.edges = counts.edges(:, 2:end);

    %%% write data to HDF5
    h5create(fname_out, '/gene_ids_segs', size(counts.gene_ids_segs), 'ChunkSize', min([size(counts.gene_ids_segs); [100 1]], [], 1), 'Deflate', 6);
    h5write(fname_out, '/gene_ids_segs', counts.gene_ids_segs);
    h5create(fname_out, '/gene_ids_edges', size(counts.gene_ids_edges), 'ChunkSize', min([size(counts.gene_ids_edges); [100 1]], [], 1), 'Deflate', 6);
    h5write(fname_out, '/gene_ids_edges', counts.gene_ids_edges);
    h5create(fname_out, '/segments', size(counts.segments), 'Chunksize', min([size(counts.segments); [100, 50]], [], 1), 'Deflate', 6);
    h5write(fname_out, '/segments', counts.segments);
    h5create(fname_out, '/seg_pos', size(counts.seg_pos), 'Chunksize', min([size(counts.seg_pos); [100, 50]], [], 1), 'Deflate', 6);
    h5write(fname_out, '/seg_pos', counts.seg_pos);
    h5create(fname_out, '/edges', size(counts.edges), 'Chunksize', min([size(counts.edges); [100, 50]], [], 1), 'Deflate', 6);
    h5write(fname_out, '/edges', counts.edges);
    h5create(fname_out, '/edge_idx', size(counts.edge_idx), 'Chunksize', min([size(counts.edge_idx); [100, 1]], [], 1), 'Deflate', 6);
    h5write(fname_out, '/edge_idx', counts.edge_idx);

    %secsave(fname_out, counts, 'counts');
