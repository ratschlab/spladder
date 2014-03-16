function count_graph_coverage_wrapper(fname_in, fname_out, CFG)

    load(fname_in);
    
    %if ~isfield(genes, 'segmentgraph'),
        genes = create_segment_graph(genes, CFG);
        secsave(fname_in, genes, 'genes');
    %end;

    if CFG.rproc,
        for s_idx = 1:length(CFG.strains),
            fprintf('%i/%i\r', j, length(CFG.strains));
            counts = count_graph_coverage(genes, CFG.list_bam(s_idx), CFG);
        end ;
    else
        chunksize = 10;
        jobinfo = rproc_empty(0) ;
        PAR = struct();
        PAR.CFG = CFG;
        for c_idx = 1:chunksize:length(genes),
            cc_idx = min(length(genes), c_idx + chunksize - 1);
            fn = strrep(fname_out, '.mat', sprintf('.chunk_%i_%i.mat', c_idx, cc_idx));
            if exist(fn, 'file'),
                continue;
            else
                fprintf(1, 'submitting chunk %i to %i\n', c_idx, cc_idx);
                PAR.genes = genes(c_idx:cc_idx);
                PAR.fn_bam = CFG.bam_fnames;
                PAR.fn_out = fn;
                PAR.CFG = CFG;
                %jobinfo(end + 1) = rproc('count_graph_coverage', PAR, 30000, CFG.options_rproc, 48*60);
                count_graph_coverage(PAR);
            end;
        end;

        jobinfo = rproc_wait(jobinfo, 30, 1, 1) ;

        %%% merge results
        counts = struct();
        for c_idx = 1:chunksize:length(genes),
            cc_idx = min(length(genes), c_idx + chunksize - 1);
            fn = strrep(fname_out, '.mat', sprintf('.chunk_%i_%i.mat', c_idx, cc_idx));
            if ~exist(fn, 'file'),
                error('Not all chunks in counting graph coverage completed!');
            else
                L = load(fn);
                if all(size(counts) == 1),
                    counts = L.counts;
                else
                    counts(:, c_idx:cc_idx) = L.counts;
                end;
            end;
        end;
    end;

    save(fname_out, 'counts');
