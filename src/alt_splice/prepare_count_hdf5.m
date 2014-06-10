function prepare_count_hdf5(CFG, fname, events, event_features)

    %%% load gene info
    if isfield(CFG, 'spladder_infile') && exist(CFG.spladder_infile, 'file'),
        load(CFG.spladder_infile);
    else
        prune_tag = '';
        if CFG.do_prune,
            prune_tag = '_pruned';
        end;
        validate_tag = '';
        if CFG.validate_splicegraphs,
            validate_tag = '.validated';
        end;
        load(sprintf('%s/spladder/genes_graph_conf%i.%s%s%s.mat', CFG.out_dirname, CFG.confidence_level, CFG.merge_strategy, validate_tag, prune_tag));
    end;

    %%% write strain and gene indices to hdf5
    h5fid = H5F.create(fname, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
    type_id = H5T.copy('H5T_C_S1');
    H5T.set_size(type_id, 'H5T_VARIABLE');
    plist = H5P.create('H5P_DATASET_CREATE');
    H5P.set_chunk(plist,2); % 2 strings per chunk

    space_id = H5S.create_simple(1, [length(CFG.strains)], H5ML.get_constant_value('H5S_UNLIMITED'));
    dset_id = H5D.create(h5fid, 'strains', type_id, space_id, plist);
    H5D.write(dset_id, type_id, 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT', CFG.strains);
    H5S.close(space_id);
    H5D.close(dset_id);

    space_id = H5S.create_simple(1, [length(event_features)], H5ML.get_constant_value('H5S_UNLIMITED'));
    dset_id = H5D.create(h5fid, 'event_features', type_id, space_id, plist);
    H5D.write(dset_id, type_id, 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT', event_features);
    H5S.close(space_id);
    H5D.close(dset_id);

    space_id = H5S.create_simple(1, [length(genes)], H5ML.get_constant_value('H5S_UNLIMITED'));
    dset_id = H5D.create(h5fid, 'gene_names', type_id, space_id, plist);
    H5D.write(dset_id, type_id, 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT', {genes.name});
    H5S.close(space_id);
    H5D.close(dset_id);
    
    space_id = H5S.create_simple(1, [length(genes)], H5ML.get_constant_value('H5S_UNLIMITED'));
    dset_id = H5D.create(h5fid, 'gene_chr', type_id, space_id, plist);
    H5D.write(dset_id, type_id, 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT', {genes.chr});
    H5S.close(space_id);
    H5D.close(dset_id);

    space_id = H5S.create_simple(1, [length(genes)], H5ML.get_constant_value('H5S_UNLIMITED'));
    dset_id = H5D.create(h5fid, 'gene_strand', type_id, space_id, plist);
    H5D.write(dset_id, type_id, 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT', {genes.strand});
    H5S.close(space_id);
    H5D.close(dset_id);
    
    H5T.close(type_id);
    H5F.close(h5fid);

    h5create(fname, '/gene_pos', [length(genes) 2], 'ChunkSize', min([length(genes) 2; [50 2]], [], 1), 'Deflate', 6);
    h5write(fname, '/gene_pos', [vertcat(genes.start) vertcat(genes.stop)]);
