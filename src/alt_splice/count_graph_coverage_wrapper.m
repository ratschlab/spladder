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
        %%% have an adaptive chunk size, that takes into account the number of strains (take as many genes as it takes to have ~10K strains)
        chunksize = max(1, floor(10000 / length(CFG.strains)));
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
                jobinfo(end + 1) = rproc('count_graph_coverage', PAR, 8000, CFG.options_rproc, 60);
                %count_graph_coverage(PAR);
            end;
        end;

        jobinfo = rproc_wait(jobinfo, 30, 1, -1) ;

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
    seg_lens = [];
    if CFG.is_half_open,
        ho_offset = 1;
    else
        ho_offset = 1;
    end;
    for i = 1:length(genes),
        if ~isempty(genes(i).segmentgraph{1}),
            seg_lens = [seg_lens, genes(i).segmentgraph{1}(2, :) - genes(i).segmentgraph{1}(1, :) + 1 - ho_offset];
        end;
    end;

    %%% write gene_names to hdf5
    h5fid = H5F.create(fname_out, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
    type_id = H5T.copy('H5T_C_S1');
    H5T.set_size(type_id, 'H5T_VARIABLE');
    space_id = H5S.create_simple(1, [length(genes)], H5ML.get_constant_value('H5S_UNLIMITED'));
    plist = H5P.create('H5P_DATASET_CREATE');
    H5P.set_chunk(plist,2); % 2 strings per chunk
    dset_id = H5D.create(h5fid, 'gene_names', type_id, space_id, plist);
    H5D.write(dset_id, type_id, 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT', {genes.name});
    H5S.close(space_id);
    H5D.close(dset_id);

    %%% write strains to hdf5
    space_id = H5S.create_simple(1, [length(CFG.strains)], H5ML.get_constant_value('H5S_UNLIMITED'));
    plist = H5P.create('H5P_DATASET_CREATE');
    H5P.set_chunk(plist,2); % 2 strings per chunk
    dset_id = H5D.create(h5fid, 'strains', type_id, space_id, plist);
    H5D.write(dset_id, type_id, 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT', CFG.strains);
    H5S.close(space_id);
    H5D.close(dset_id);

    H5F.close(h5fid);

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
    h5create(fname_out, '/seg_len', size(seg_lens), 'Chunksize', min([size(seg_lens); [100, 50]], [], 1), 'Deflate', 6);
    h5write(fname_out, '/seg_len', seg_lens);


    %secsave(fname_out, counts, 'counts');
