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
        fn_out_count = strrep(fn_out, '.mat', '.counts.mat');
        fn_out_info = strrep(fn_out, '.mat', '.info.mat');

        %%% check, if confirmed version exists
        if ~exist(fn_out_info, 'file'),

            load(fn_out);

            %%% add strain information, so we can do two way chunking!
            if ~exist('events_all_strains', 'var'),
                [events_all, events_all_strains] = add_strains(events_all, CFG);
                save(fn_out, 'events_all', 'events_all_strains');
            end;
            
            if ~CFG.rproc,
                events_all = verify_all_events(events_all, 1:length(CFG.strains), CFG.bam_fnames(replicate, :), 1, event_type, set_confidence_level(CFG)) ;
            else
                jobinfo = rproc_empty() ;
                chunk_size_events = 25 ;
                chunk_size_strains = 50 ;
                job_nr = 1;
                for i = 1:chunk_size_events:length(events_all),
                    idx_events = i:min(i+chunk_size_events-1, length(events_all)) ;
                    for j = 1:chunk_size_strains:length(strains),
                        idx_strains = j:min(j+chunk_size_strains-1, length(CFG.strains));
                        PAR.ev = events_all(idx_events) ;
                        PAR.strain_idx = idx_strains ;
                        PAR.list_bam = CFG.bam_fnames(replicate, :) ;
                        PAR.verify = 1 ;
                        PAR.out_fn = sprintf('%s/event_count_chunks/%s_%i_%i_R%i_C%i.mat', CFG.out_dirname, event_type, i, j, replicate, CFG.confidence_level);
                        PAR.event_type = event_type;
                        PAR.conf_filter = set_confidence_level(CFG);
                        if exist(PAR.out_fn, 'file')
                            fprintf('Chunk event %i, strain %i already completed\n', i, j);
                        else
                            fprintf('Submitting job %i, event chunk %i, strain chunk %i\n', job_nr, i, j);
                            jobinfo(job_nr) = rproc('verify_all_events', PAR, 8000, CFG.options_rproc, 60) ;
                            %verify_all_events(PAR);
                            job_nr = job_nr + 1;
                        end ;
                    end ;
                end ;
                
                [jobinfo nr_crashed] = rproc_wait(jobinfo, 20, 1, 1) ;
                
                events_all_ = [];
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
                                ev_(jj).info(idx_strains) = ev(jj).info;
                            end;
                        end;
                    end;
                    events_all_ = [events_all_ ev_];
                end ;
                assert(length(events_all) == length(events_all_)) ;
                if strcmp(event_type, 'alt_3prime') || strcmp(event_type, 'alt_5prime'),
                    assert(isequal([events_all.intron1_col], [events_all_.intron1_col])) ;
                    assert(isequal([events_all.intron2_col], [events_all_.intron2_col])) ;
                end;
                
                events_all = events_all_ ;
            end ;
            
            for i = 1:length(events_all),
                events_all(i).num_verified = sum(events_all(i).verified, 1) ;
                events_all(i).confirmed = min(events_all(i).num_verified) ;
            end ;
            
            verified_count = [] ;
            for min_verified = 1:length(CFG.strains),
                verified_count(min_verified) = sum([events_all.confirmed] >= min_verified) ;
            end ;
            
            min_verified = 1 ;
            events_confirmed = events_all([events_all.confirmed] >= min_verified) ;
            
            %%% save events per chromosome 
            chunksize = 200;
            events_confirmed_ = events_confirmed;
            events_all_ = events_all;
            max_len_all = length(events_all_);
            max_len_confirmed = length(events_confirmed_);

            for c1_idx = 1:chunksize:max_len_all,
                fn_out_count_chunk = strrep(fn_out_count, '.mat', sprintf('.chunk%i.mat', c1_idx));
                if ~exist(fn_out_count_chunk, 'file'),
                    tmp_idx = c1_idx:min(c1_idx + chunksize - 1, max_len_all);
                    events_all = events_all_(tmp_idx);
                    fprintf('saving counted %s events to %s\n', event_type, fn_out_count_chunk) ;
                    save(fn_out_count_chunk, 'events_all');
                else
                    fprintf('%s exists\n', fn_out_count_chunk);
                end;
            end;

            for c1_idx = 1:chunksize:max_len_confirmed,
                fn_out_conf_chunk = strrep(fn_out_conf, '.mat', sprintf('.chunk%i.mat', c1_idx));
                if ~exist(fn_out_conf_chunk, 'file'),
                    tmp_idx = c1_idx:min(c1_idx + chunksize - 1, max_len_confirmed);
                    events_confirmed = events_confirmed_(tmp_idx);
                    fprintf('saving confirmed %s events to %s\n', event_type, fn_out_conf_chunk) ;
                    save(fn_out_conf_chunk, 'events_confirmed');
                else
                    fprintf('%s exists\n', fn_out_conf_chunk);
                end;
            end;

            %%% save summary file
            save(fn_out_info, 'chunksize', 'max_len_all', 'max_len_confirmed');
            events_confirmed = events_confirmed_;
            event_all = events_all_;
        else
            fprintf('%s exists - loading chunk-wise!\n\n', fn_out_info);
            load(fn_out_info);
            events_confirmed_ = [];
            for c1_idx = 1:chunksize:max_len_confirmed,
                fprintf('...loading chunk %i\n', c1_idx);
                fn_out_conf_chunk = strrep(fn_out_conf, '.mat', sprintf('.chunk%i.mat', c1_idx));
                tic
                load(fn_out_conf_chunk);
                events_confirmed_ = [events_confirmed_ events_confirmed];
                toc
            end;
            events_confirmed = events_confirmed_;

            events_all_ = [];
            for c1_idx = 1:chunksize:max_len_all,
                fprintf('...loading chunk %i\n', c1_idx);
                fn_out_count_chunk = strrep(fn_out_count, '.mat', sprintf('.chunk%i.mat', c1_idx));
                tic
                load(fn_out_count_chunk);
                events_all_ = [events_all_ events_all];
                toc
            end;
            events_all = events_all_;
        end;

        fn_out_txt = strrep(fn_out, '.mat', '.txt') ;
        fn_out_conf_txt = strrep(fn_out_conf, '.mat', '.txt') ;
        fn_out_conf_tcga = strrep(fn_out_conf, '.mat', '.tcga.txt') ;
        fn_out_conf_gff3 = strrep(fn_out_conf, '.mat', '.gff3') ;
        fn_out_conf_genes_alt = strrep(fn_out_conf, '.mat', '.genes.mat') ;

        if isempty(events_all),
            fprintf('No %s event could be found. - Nothing  to report.\n', event_type);
            continue;
        else
            fprintf('Reporting complete %s events.\n', event_type);
        end;
        if exist(fn_out_txt, 'file')
            fprintf('%s already exists\n', fn_out_txt);
        else
            write_events_txt(fn_out_txt, CFG.strains, events_all) ;
        end;


        if isempty(events_confirmed),
            fprintf('No %s event could be confirmed. - Nothing  to report.\n', event_type);
            continue;
        else
            fprintf('Reporting confirmed %s events.\n', event_type);
        end;

        if exist(fn_out_conf_gff3, 'file')
            fprintf('%s already exists\n', fn_out_conf_gff3);
        else
            write_events_gff3(fn_out_conf_gff3, events_confirmed) ;
        end;

        if exist(fn_out_conf_txt, 'file')
            fprintf('%s already exists\n', fn_out_conf_txt);
        else
            write_events_txt(fn_out_conf_txt, CFG.strains, events_confirmed) ;
        end;

        if exist(fn_out_conf_tcga, 'file')
            fprintf('%s already exists\n', fn_out_conf_tcga);
        else
            write_events_tcga(fn_out_conf_tcga, CFG.strains, events_confirmed) ;
        end;

        fn_out_conf_txt = strrep(fn_out_conf, '.mat', '.filt0.05.txt') ;
        if exist(fn_out_conf_txt, 'file')
            fprintf('%s already exists\n', fn_out_conf_txt);
        else
            fprintf('writing filtered events (sample freq 0.05)');
            cf_idx = ([events_confirmed.confirmed] >= 0.05 * length(events_confirmed(1).detected));
            events_confirmed = events_confirmed(cf_idx);
            write_events_txt(fn_out_conf_txt, CFG.strains, events_confirmed);
        end;

        fn_out_conf_txt = strrep(fn_out_conf, '.mat', '.filt0.1.txt') ;
        if exist(fn_out_conf_txt, 'file')
            fprintf('%s already exists\n', fn_out_conf_txt);
        else
            fprintf('writing filtered events (sample freq 0.01)');
            cf_idx = ([events_confirmed.confirmed] >= 0.1 * length(events_confirmed(1).detected));
            events_confirmed = events_confirmed(cf_idx);
            write_events_txt(fn_out_conf_txt, CFG.strains, events_confirmed);
        end;

        %genes_alt3 = get_ALT_END_genes(CFG.strains, events_confirmed) ;
        %fprintf('saving %s events to %s\n', event_type, fn_out_conf_genes_alt) ;
        %save(fn_out_conf_genes_alt, 'genes_alt3', '-v7.3') ;
    end;
