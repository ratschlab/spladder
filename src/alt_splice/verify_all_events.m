function ev = verify_all_events(ev, event_type, CFG, strain_idx) ;
% ev = verify_all_events(ev, event_type, CFG, strain_idx) ;

out_fn = '' ;

fieldnames = {'intron', 'exon1', 'exon2', 'exon_alt1', 'exon_alt2', 'exon_const', 'intron1', 'intron2', 'exon', 'exon_pre', 'exon_aft', 'exons'};

%%% set parameters if called by rproc
if nargin==1,
    PAR = ev ;
    ev = PAR.ev ;
    event_type = PAR.event_type;
    CFG = PAR.CFG;
    if isfield(PAR, 'out_fn'),
        out_fn = PAR.out_fn ;
    end ;
    if isfield(PAR, 'strain_idx'),
        strain_idx = PAR.strain_idx;
    else
        strain_idx = 1:length(CFG.strains);
    end;
elseif nargin < 4,
    strain_idx = 1:length(CFG.strains);
end ;

%%% load counts
prune_tag = '';
if CFG.do_prune,
    prune_tag = '_pruned';
end;
validate_tag = '';
if CFG.validate_splicegraphs,
    validate_tag = '.validated';
end;

%%% verify the events if demanded
if CFG.verify_alt_events,
    if isfield(CFG, 'spladder_infile') && exist(CFG.spladder_infile, 'file'),
        fn_count = strrep(CFG.spladder_infile, 'mat', 'count.mat');
        fn_genes = CFG.spladder_infile;
    else
        fn_count = sprintf('%s/spladder/genes_graph_conf%i.%s%s%s.count.mat', CFG.out_dirname, CFG.confidence_level, CFG.merge_strategy, validate_tag, prune_tag);
        fn_genes = sprintf('%s/spladder/genes_graph_conf%i.%s%s%s.mat', CFG.out_dirname, CFG.confidence_level, CFG.merge_strategy, validate_tag, prune_tag);
    end;
    %%% load count index data from hdf5
    gene_ids_segs = h5read(fn_count, '/gene_ids_segs');
    gene_ids_edges = h5read(fn_count, '/gene_ids_edges');

    load(fn_genes);
    
    %%% sort events by gene idx
    [tmp, s_idx] = sort([ev.gene_idx]);
    ev = ev(s_idx);
    [tmp, old_idx] = sort(s_idx);

    %%% find gene idx boundaries
    assert(isequal(gene_ids_segs, sort(gene_ids_segs)));
    assert(isequal(gene_ids_edges, sort(gene_ids_edges)));

    [tmp, genes_f_idx_segs, tmp] = unique(gene_ids_segs);
    genes_l_idx_segs = [genes_f_idx_segs(2:end) - 1; length(gene_ids_segs)];

    [tmp, genes_f_idx_edges, tmp] = unique(gene_ids_edges);
    genes_l_idx_edges = [genes_f_idx_edges(2:end) - 1; length(gene_ids_edges)];

    gr_idx_segs = 1;
    gr_idx_edges = 1;
    for i = 1:length(ev),
        fprintf('%i/%i\r', i, length(ev));
        g_idx = ev(i).gene_idx;
        while gene_ids_segs(genes_f_idx_segs(gr_idx_segs)) < g_idx,
            gr_idx_segs = gr_idx_segs + 1;
        end;
        assert(gene_ids_segs(genes_f_idx_segs(gr_idx_segs)) == g_idx);

        while gene_ids_edges(genes_f_idx_edges(gr_idx_edges)) < g_idx,
            gr_idx_edges = gr_idx_edges + 1;
        end;
        assert(gene_ids_edges(genes_f_idx_edges(gr_idx_edges)) == g_idx);

        span_segs = genes_l_idx_segs(gr_idx_segs) - genes_f_idx_segs(gr_idx_segs) + 1;
        span_edges = genes_l_idx_edges(gr_idx_edges) - genes_f_idx_edges(gr_idx_edges) + 1;
        
        %%% laod relevant count data from HDF5
        plist = 'H5P_DEFAULT';
        fid = H5F.open(fn_count);
        %%% segments
        dset_id = H5D.open(fid,'/segments');
        dims = fliplr([span_segs length(strain_idx)]);
        file_space_id = H5D.get_space(dset_id);
        mem_space_id = H5S.create_simple(2, dims, []);
        %%% H5S.select_hyperslab((file_space_id, operation, offset, [], [], block) 
        H5S.select_hyperslab(file_space_id, 'H5S_SELECT_SET', fliplr([genes_f_idx_segs(gr_idx_segs) - 1 strain_idx(1) - 1]),[], [], dims);
        segments = H5D.read(dset_id, 'H5ML_DEFAULT', mem_space_id, file_space_id, plist);
        H5D.close(dset_id);
        %%% seg_pos
        dset_id = H5D.open(fid,'/seg_pos');
        dims = fliplr([span_segs length(strain_idx)]);
        file_space_id = H5D.get_space(dset_id);
        mem_space_id = H5S.create_simple(2, dims, []);
        H5S.select_hyperslab(file_space_id, 'H5S_SELECT_SET', fliplr([genes_f_idx_segs(gr_idx_segs) - 1 strain_idx(1) - 1]),[], [], dims);
        seg_pos = H5D.read(dset_id, 'H5ML_DEFAULT', mem_space_id, file_space_id, plist);
        H5D.close(dset_id);
        %%% edges
        dset_id = H5D.open(fid,'/edges');
        dims = fliplr([span_edges length(strain_idx)]);
        file_space_id = H5D.get_space(dset_id);
        mem_space_id = H5S.create_simple(2, dims, []);
        H5S.select_hyperslab(file_space_id, 'H5S_SELECT_SET', fliplr([genes_f_idx_edges(gr_idx_edges) - 1 strain_idx(1) - 1]),[], [], dims);
        edges = H5D.read(dset_id, 'H5ML_DEFAULT', mem_space_id, file_space_id, plist);
        H5D.close(dset_id);
        %%% edge_idx
        dset_id = H5D.open(fid,'/edge_idx');
        dims = fliplr([span_edges 1]);
        file_space_id = H5D.get_space(dset_id);
        mem_space_id = H5S.create_simple(2, dims, []);
        H5S.select_hyperslab(file_space_id, 'H5S_SELECT_SET', fliplr([genes_f_idx_edges(gr_idx_edges) - 1 0]),[], [], dims);
        edge_idx = H5D.read(dset_id, 'H5ML_DEFAULT', mem_space_id, file_space_id, plist);
        H5D.close(dset_id);
        H5F.close(fid);

        for s_idx = 1:length(strain_idx),
            ev_tmp = ev(i) ;
            for f_idx = 1:length(fieldnames),
                field = fieldnames{f_idx};
                if isfield(ev, field) && size(ev_tmp.(field), 1) > 1,
                    ev_tmp.(field) = ev_tmp.(field)(s_idx, :);
                end;
            end;
            if strcmp(event_type, 'exon_skip'),
                [ev(i).verified(s_idx,:), ev(i).info(s_idx, :)] = verify_exon_skip(ev_tmp, genes(g_idx), segments(:, s_idx)', [edge_idx edges(:, s_idx)], CFG) ;
            elseif strcmp(event_type, 'alt_3prime') || strcmp(event_type, 'alt_5prime') || strcmp(event_type, 'alt_prime'),
                [ev(i).verified(s_idx,:), ev(i).info(s_idx, :)] = verify_alt_prime(ev_tmp, genes(g_idx), segments(:, s_idx)', [edge_idx edges(:, s_idx)], CFG) ;
            elseif strcmp(event_type, 'intron_retention'),
                [ev(i).verified(s_idx,:), ev(i).info(s_idx, :)] = verify_intron_retention(ev_tmp, genes(g_idx), segments(:, s_idx)', [edge_idx edges(:, s_idx)], seg_pos(:, s_idx)', CFG) ;
            elseif strcmp(event_type, 'mult_exon_skip'),
                [ev(i).verified(s_idx,:), ev(i).info(s_idx, :)] = verify_mult_exon_skip(ev_tmp, genes(g_idx), segments(:, s_idx)', [edge_idx edges(:, s_idx)], CFG) ;
            end;
        end ;
    end ;
    ev = ev(old_idx);
end ;

if ~isempty(out_fn),
    secsave(out_fn, ev, 'ev');
end;
