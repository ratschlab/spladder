function alt_genes_analyze(CFG, event_type)
% function alt_genes_analyze(CFG, event_type)

    if CFG.rproc,
        unix(sprintf('mkdir -p %s/event_count_chunks', CFG.out_dirname));
    end ;

    for replicate = CFG.replicate_idxs,
        
        fprintf(1, 'confidence %i / replicate %i\n', CFG.confidence_level, replicate);

        if length(CFG.replicate_idxs) > 1,
            rep_tag = sprintf('_R%i', r_idx);
        else
            rep_tag = '';
        end;

        fn_out = sprintf('%s/%s_%s%s_C%i.mat', CFG.out_dirname, CFG.merge_strategy, event_type, rep_tag, CFG.confidence_level) ;
        fn_out_conf = strrep(fn_out, '.mat', '.confirmed.mat');
        fn_out_count = strrep(fn_out, '.mat', '.counts.hdf5');

        %%% define result files
        fn_out_txt = strrep(fn_out, '.mat', '.txt') ;
        fn_out_conf_txt = strrep(fn_out_conf, '.mat', '.txt') ;
        fn_out_conf_tcga = strrep(fn_out_conf, '.mat', '.tcga.txt') ;
        fn_out_conf_gff3 = strrep(fn_out_conf, '.mat', '.gff3') ;
        %fn_out_conf_genes_alt = strrep(fn_out_conf, '.mat', '.genes.mat') ;

        %%% check if there is anything to do
        if exist(fn_out_txt, 'file') && exist(fn_out_conf_txt, 'file') && exist(fn_out_conf_tcga, 'file') && exist(fn_out_conf_gff3, 'file'),
            fprintf('All output files for %s exist.\n\n', event_type);
            continue;
        end;

        if strcmp(event_type, 'exon_skip'),
            event_features = {'valid', 'exon_cov', 'exon_pre_cov', 'exon_aft_cov', 'exon_pre_exon_conf', 'exon_exon_aft_conf', 'exon_pre_exon_aft_conf'};
        elseif strcmp(event_type, 'alt_3prime') || strcmp(event_type, 'alt_5prime'),
            event_features =  {'valid', 'exon_diff_cov', 'exon_const_cov', 'intron1_conf', 'intron2_conf'};
        elseif strcmp(event_type, 'intron_retention'),
            event_features = {'valid', 'intron_cov', 'exon1_cov', 'exon2_cov', 'intron_conf', 'intron_cov_region'}; 
        elseif strcmp(event_type, 'mult_exon_skip'),
            event_features = {'valid', 'exon_pre_cov', 'exons_cov', 'exon_aft_cov', 'exon_pre_exon_conf', 'exon_exon_aft_conf', 'exon_pre_exon_aft_conf', 'sum_inner_exon_conf', 'num_inner_exon'}; 
        end;

        %%% check, if confirmed version exists
        if ~exist(fn_out_count, 'file'),

            load(fn_out);

            %%% write strain and gene indices to hdf5
            prepare_count_hdf5(CFG, fn_out_count, events_all, event_features);

            %%% handle case where we did not find any event of this type
            if isempty([events_all.event_type]),
                h5create(fn_out_count, '/event_counts', [1]);
                h5write(fn_out_count, '/event_counts', 0);
            else
                %%% add strain information, so we can do two way chunking!
                if ~exist('events_all_strains', 'var'),
                    [events_all, events_all_strains] = add_strains(events_all, CFG);
                    secsave(fn_out, {events_all, events_all_strains}, {'events_all', 'events_all_strains'});
                end;
                
                if ~CFG.rproc,
                    events_all = verify_all_events(events_all, event_type, CFG) ;
                    dims = [length(CFG.strains) length(event_features) size(events_all, 2)];
                    h5create(fn_out_count, '/event_counts', dims, 'ChunkSize', min([dims; [50 50 10]], [], 1), 'Deflate', 6);
                    h5write(fn_out_count, '/event_counts', cat(3, events_all.info));
                    gene_idx_ = [events_all.gene_idx];
                    h5create(fn_out_count, '/gene_idx', [length(gene_idx_)]);
                    h5write(fn_out_count, '/gene_idx', gene_idx_);
                    events_all = rmfield(events_all, 'info');
                else
                    jobinfo = rproc_empty() ;
                    chunk_size_events = 1000 ;
                    chunk_size_strains = 500 ;
                    job_nr = 1;
                    for i = 1:chunk_size_events:length(events_all),
                        idx_events = i:min(i+chunk_size_events-1, length(events_all)) ;
                        for j = 1:chunk_size_strains:length(CFG.strains),
                            idx_strains = j:min(j+chunk_size_strains-1, length(CFG.strains));
                            PAR.ev = events_all(idx_events) ;
                            PAR.strain_idx = idx_strains ;
                            PAR.list_bam = CFG.bam_fnames(:, replicate) ;
                            PAR.out_fn = sprintf('%s/event_count_chunks/%s_%i_%i_R%i_C%i.mat', CFG.out_dirname, event_type, i, j, replicate, CFG.confidence_level);
                            PAR.event_type = event_type;
                            PAR.CFG = CFG;
                            if exist(PAR.out_fn, 'file')
                                fprintf('Chunk event %i, strain %i already completed\n', i, j);
                            else
                                fprintf('Submitting job %i, event chunk %i (%i), strain chunk %i (%i)\n', job_nr, i, length(events_all), j, length(CFG.strains));
                                jobinfo(job_nr) = rproc('verify_all_events', PAR, 10000, CFG.options_rproc, 60) ;
                                %verify_all_events(PAR);
                                job_nr = job_nr + 1;
                            end ;
                        end ;
                    end ;
                    
                    [jobinfo nr_crashed] = rproc_wait(jobinfo, 20, 1, -1) ;
                    
                    %%% open event count file to directly write counts into common hdf5
                    dims = [length(CFG.strains) length(event_features) size(events_all, 2)];
                    h5create(fn_out_count, '/event_counts', dims, 'ChunkSize', min([dims; [50 50 10]], [], 1), 'Deflate', 6);

                    events_all_ = [];
                    gene_idx_ = [];
                    fprintf('Collecting results from chunks ...\n');
                    for i = 1:chunk_size_events:length(events_all),
                        idx_events = i:min(i+chunk_size_events-1, length(events_all)) ;
                        ev_ = [];
                        for j = 1:chunk_size_strains:length(CFG.strains),
                            idx_strains = j:min(j+chunk_size_strains-1, length(CFG.strains));
                            fprintf('\r%i (%i), %i (%i)', i, length(events_all), j, length(CFG.strains));
                            out_fn = sprintf('%s/event_count_chunks/%s_%i_%i_R%i_C%i.mat', CFG.out_dirname, event_type, i, j, replicate, CFG.confidence_level);
                            if ~exist(out_fn, 'file'),
                                error(sprintf('not finished %s\n', out_fn)) ;
                            end ;
                            load(out_fn);
                            if j == 1,
                                ev_ = ev;
                            else 
                                for jj = 1:length(ev_),
                                    ev_(jj).verified(idx_strains, :) = ev(jj).verified;
                                    ev_(jj).info(idx_strains, :) = ev(jj).info;
                                end;
                            end;
                        end;
                        for jj = 1:length(ev_),
                            h5write(fn_out_count, '/event_counts', ev_(jj).info, [1 1 (i + jj - 1)], [size(ev_(jj).info) 1]);
                            gene_idx_(i + jj - 1) = ev_(jj).gene_idx;
                        end;
                        ev_ = rmfield(ev_, 'info');
                        events_all_ = [events_all_ ev_];
                    end ;
                    assert(length(events_all) == length(events_all_)) ;
                    if strcmp(event_type, 'alt_3prime') || strcmp(event_type, 'alt_5prime'),
                        assert(isequal([events_all.intron1_col], [events_all_.intron1_col])) ;
                        assert(isequal([events_all.intron2_col], [events_all_.intron2_col])) ;
                    end;
                    h5create(fn_out_count, '/gene_idx', [length(gene_idx_)]);
                    h5write(fn_out_count, '/gene_idx', gene_idx_);
                    events_all = events_all_ ;
                end ;
                
                %%% write more event infos to hdf5
                if strcmp(event_type, 'exon_skip'),
                    event_pos = [vertcat(events_all.exon_pre), vertcat(events_all.exon), vertcat(events_all.exon_aft)];
                elseif strcmp(event_type, 'intron_retention'),
                    event_pos = [vertcat(events_all.exon1), vertcat(events_all.exon2)];
                elseif strcmp(event_type, 'alt_3prime') || strcmp(event_type, 'alt_5prime'),
                    tmp = zeros(length(events_all), 6);
                    for tt = 1:length(events_all),
                        ttmp = sortrows([events_all(tt).exon_alt1; events_all(tt).exon_alt2; events_all(tt).exon_const])';
                        tmp(tt, :) = ttmp(:)';
                    end;
                    event_pos = tmp;
                elseif strcmp(event_type, 'mult_exon_skip'),
                    tmp = zeros(length(events_all), 2);
                    for tt = 1:length(events_all),
                        tmp(tt, 1) = events_all(tt).exons(1, 1);
                        tmp(tt, 2) = events_all(tt).exons(end, end);
                    end;
                    event_pos = [vertcat(events_all.exon_pre), tmp, vertcat(events_all.exon_aft)];
                end;
                h5create(fn_out_count, '/event_pos', size(event_pos));
                h5write(fn_out_count, '/event_pos', event_pos);

                for i = 1:length(events_all),
                    events_all(i).num_verified = sum(events_all(i).verified, 1) ;
                    events_all(i).confirmed = min(events_all(i).num_verified) ;
                end ;
                
                num_verified = vertcat(events_all.num_verified);

                %verified_count = [] ;
                %for min_verified = 1:length(CFG.strains),
                %    verified_count(min_verified) = sum([events_all.confirmed] >= min_verified) ;
                %end ;
                
                confirmed_idx = find([events_all.confirmed] >= 1) ;

                %%% save events
                secsave(fn_out, {events_all, events_all_strains}, {'events_all', 'events_all_strains'});
                secsave(fn_out_conf, confirmed_idx, 'confirmed_idx');
                if isempty(confirmed_idx),
                    h5create(fn_out_count, '/conf_idx', 1);
                    h5write(fn_out_count, '/conf_idx', -1);
                else
                    h5create(fn_out_count, '/conf_idx', [length(confirmed_idx)]);
                    h5write(fn_out_count, '/conf_idx', confirmed_idx);
                end;
                h5create(fn_out_count, '/verified', size(num_verified));
                h5write(fn_out_count, '/verified', num_verified);
            end;
        else
            load(fn_out);
            load(fn_out_conf);
        end;

        if isempty(events_all) || isempty([events_all.event_type]),
            fprintf('No %s event could be found. - Nothing to report.\n', event_type);
            continue;
        else
            fprintf('Reporting complete %s events.\n', event_type);
        end;
        if CFG.output_txt,
            if exist(fn_out_txt, 'file')
                fprintf('%s already exists\n', fn_out_txt);
            else
                write_events_txt(fn_out_txt, CFG.strains, events_all, fn_out_count, 1:length(events_all), CFG.global_gene_list) ;
            end;
        end;

        if isempty(confirmed_idx),
            fprintf('No %s event could be confirmed. - Nothing to report.\n', event_type);
            continue;
        else
            fprintf('Reporting confirmed %s events.\n', event_type);
        end;

        if CFG.output_confirmed_gff3,
            if exist(fn_out_conf_gff3, 'file')
                fprintf('%s already exists\n', fn_out_conf_gff3);
            else
                write_events_gff3(fn_out_conf_gff3, events_all, confirmed_idx);
            end;
        end;

        if CFG.output_confirmed_txt,
            if exist(fn_out_conf_txt, 'file')
                fprintf('%s already exists\n', fn_out_conf_txt);
            else
                write_events_txt(fn_out_conf_txt, CFG.strains, events_all, fn_out_count, confirmed_idx, CFG.global_gene_list) ;
            end;
        end;

        if CFG.output_confirmed_tcga,
            if exist(fn_out_conf_tcga, 'file')
                fprintf('%s already exists\n', fn_out_conf_tcga);
            else
                write_events_tcga(fn_out_conf_tcga, CFG.strains, events_all, fn_out_count, confirmed_idx) ;
            end;
        end;

        if CFG.output_filtered_txt,
            fn_out_conf_txt = strrep(fn_out_conf, '.mat', '.filt0.05.txt') ;
            if exist(fn_out_conf_txt, 'file')
                fprintf('%s already exists\n', fn_out_conf_txt);
            else
                fprintf('writing filtered events (sample freq 0.05)');
                cf_idx = ([events_all(confirmed_idx).confirmed] >= 0.05 * length(events_all(1).detected));
                write_events_txt(fn_out_conf_txt, CFG.strains, events_all, fn_out_count, confirmed_idx(cf_idx), CFG.global_gene_list);
            end;

            fn_out_conf_txt = strrep(fn_out_conf, '.mat', '.filt0.1.txt') ;
            if exist(fn_out_conf_txt, 'file')
                fprintf('%s already exists\n', fn_out_conf_txt);
            else
                fprintf('writing filtered events (sample freq 0.01)');
                cf_idx = ([events_all(confirmed_idx).confirmed] >= 0.1 * length(events_all(1).detected));
                write_events_txt(fn_out_conf_txt, CFG.strains, events_all, fn_out_count, confirmed_idx(cf_idx), CFG.global_gene_list);
            end;
        end;

        %genes_alt3 = get_ALT_END_genes(CFG.strains, events_confirmed) ;
        %fprintf('saving %s events to %s\n', event_type, fn_out_conf_genes_alt) ;
        %save(fn_out_conf_genes_alt, 'genes_alt3', '-v7.3') ;
    end;
