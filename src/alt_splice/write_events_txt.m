function write_ES_txt(fn_out_txt, strains, events, anno_fn)
    
    if isempty(events),
        fprintf('No events present.\n');
        return
    end;

    if nargin > 3,
       load(anno_fn);
       anno = genes;
       [anno_names s_idx] = sort({genes.name});
       genes = genes(s_idx);
       clear genes
    end;

    fprintf('writing %s events in flat txt format to %s\n', events(1).event_type, fn_out_txt);

    fd = fopen(fn_out_txt, 'w+') ;
    if nargin > 3,
        gene_header = sprintf('\tgene_start\tgene_end');
    else
        gene_header = '';
    end;

    if strcmp(events(1).event_type, 'exon_skip'),
        fprintf(fd, 'contig\tstrand\tevent_id\tgene_name%s\texon_pre_start\texon_pre_end\texon_start\texon_end\texon_aft_start\texon_aft_end', gene_header) ;
        tmp = reshape(repmat(strains, 6, 1), 1, length(strains) * 6);
        fprintf(fd, '\t%s:exon_pre_cov\t%s:exon_cov\t%s:exon_aft_cov\t%s:intron_pre_conf\t%s:intron_aft_conf\t%s:intron_skip_conf', tmp{:}) ;
        fprintf(fd, '\n') ;
    elseif strcmp(events(1).event_type, 'alt_3prime') || strcmp(events(1).event_type, 'alt_5prime'),
        fprintf(fd, 'contig\tstrand\tevent_id\tgene_name%s\texon_const_start\texon_const_end\texon_alt1_start\texon_alt1_end\texon_alt2_start\texon_alt2_end', gene_header) ;
        tmp = reshape(repmat(strains, 4, 1), 1, length(strains) * 4);
        fprintf(fd, '\t%s:exon_diff_cov\t%s:exon_const_cov\t%s:intron1_conf\t%s:intron2_conf', tmp{:}) ;
        fprintf(fd, '\n') ;
    elseif strcmp(events(1).event_type, 'intron_retention'),
        fprintf(fd, 'contig\tstrand\tevent_id\tgene_name%s\texon1_start\texon1_end\tintron_start\tintron_end\texon2_start\texon2_end', gene_header) ;
        tmp = reshape(repmat(strains, 4, 1), 1, length(strains) * 4);
        fprintf(fd, '\t%s:exon1_cov\t%s:intron_cov\t%s:exon2_cov\t%s:intron_conf', tmp{:}) ;
        fprintf(fd, '\n') ;
    elseif strcmp(events(1).event_type, 'mult_exon_skip'),
        fprintf(fd, 'contig\tstrand\tevent_id\tgene_name%s\texon_pre_start\texon_pre_end\texon_starts\texon_ends\texon_aft_start\texon_aft_end', gene_header) ;
        tmp = reshape(repmat(strains, 8, 1), 1, length(strains) * 8);
        fprintf(fd, '\t%s:exon_pre_cov\t%s:exon_cov\t%s:exon_aft_cov\t%s:intron_pre_conf\t%s:intron_inner_conf\t%s:exon_inner_count\t%s:intron_aft_conf\t%s:intron_skip_conf', tmp{:}) ;
        fprintf(fd, '\n') ;
    else
        error(sprintf('Unknown event type: %s\n', events(1).event_type));
    end;

    for i = 1:length(events),
        fprintf(fd, '%s\t%c\t%s_%i\t%s', events(i).chr, events(i).strand, events(i).event_type, i, events(i).gene_name{1});
        if nargin > 3,
            a_idx = strmatch(events(i).gene_name{1}, anno_names, 'exact');
            assert(~isempty(a_idx));
            fprintf(fd, '\t%i\t%i', anno(a_idx).start, anno(a_idx).stop);
        end;
        ev = events(i);
        if strcmp(events(i).event_type, 'exon_skip'),
            fprintf(fd, '\t%i\t%i\t%i\t%i\t%i\t%i', ev.exon_pre_col(1), ev.exon_pre_col(2), ev.exon_col(1), ev.exon_col(2), ...
                    ev.exon_aft_col(1), ev.exon_aft_col(2)) ;
            for j = 1:length(strains),
                if ev.info(j).valid,
                    fprintf(fd, '\t%.1f\t%.1f\t%.1f\t%i\t%i\t%i', ev.info(j).exon_pre_cov, ev.info(j).exon_cov, ev.info(j).exon_aft_cov, ...
                            ev.info(j).exon_pre_exon_conf, ev.info(j).exon_exon_aft_conf, ev.info(j).exon_pre_exon_aft_conf) ;
                else
                    fprintf(fd, '\t-1\t-1\t-1\t-1') ;
                end ;
            end ;
        elseif strcmp(ev.event_type, 'intron_retention'),
            fprintf(fd, '\t%i\t%i\t%i\t%i\t%i\t%i', ev.exon1_col(1), ev.exon1_col(2), ev.intron_col(1), ev.intron_col(2), ...
                    ev.exon2_col(1), ev.exon2_col(2)) ;
            for j = 1:length(strains),
                if ev.info(j).valid,
                    fprintf(fd, '\t%.1f\t%.1f\t%.1f\t%i', ev.info(j).exon1_cov, ev.info(j).intron_cov, ev.info(j).exon2_cov, ev.info(j).intron_conf) ;
                else
                    fprintf(fd, '\t-1\t-1\t-1') ;
                end ;
            end ;
        elseif strcmp(ev.event_type, 'alt_3prime') || strcmp(ev.event_type, 'alt_5prime'),
            fprintf(fd, '\t%i\t%i\t%i\t%i\t%i\t%i', ev.exon_const_col(1), ev.exon_const_col(2), ev.exon_alt1_col(1), ...
                    ev.exon_alt1_col(2), ev.exon_alt2_col(1), ev.exon_alt2_col(2)) ;
            for j = 1:length(strains),
                if ev.info(j).valid,
                    fprintf(fd, '\t%1.1f\t%1.1f\t%i\t%i', ev.info(j).exon_diff_cov, ev.info(j).exon_const_cov, ...
                            ev.info(j).intron1_conf, ev.info(j).intron2_conf) ;
                else
                    fprintf(fd, '\t-1\t-1\t-1\t-1') ;
                end ;
            end ;
        elseif strcmp(ev.event_type, 'mult_exon_skip'),
            fprintf(fd, '\t%i\t%i', ev.exon_pre_col(1), ev.exon_pre_col(2));
            starts = sprintf('%i', ev.exons_col(1));
            ends = sprintf('%i', ev.exons_col(2));
            for k = 3:2:size(ev.exons_col, 2),
                starts = sprintf('%s:%i', starts, ev.exons_col(k));
                ends = sprintf('%s:%i', ends, ev.exons_col(k + 1));
            end;
            fprintf(fd, '%s\t%s\t%i\t%i', starts, ends, ev.exon_aft_col(1), ev.exon_aft_col(2)) ;
            for j = 1:length(strains),
                if ev.info(j).valid,
                    fprintf(fd, '\t%.1f\t%.1f\t%.1f\t%i\t%i\t%i\t%i\t%i', ev.info(j).exon_pre_cov, ev.info(j).exons_cov, ev.info(j).exon_aft_cov, ev.info(j).exon_pre_exon_conf, ...
                            ev.info(j).sum_inner_exon_conf, ev.info(j).num_inner_exon, ev.info(j).exon_exon_aft_conf, ev.info(j).exon_pre_exon_aft_conf) ;
                else
                    fprintf(fd, '\t-1\t-1\t-1\t-1\t-1\t-1\t-1') ;
                end ;
            end ;
        end;
        fprintf(fd, '\n') ;
    end ;
    fclose(fd) ;
