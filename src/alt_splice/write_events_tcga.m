function write_events_tcga(fn_out, strains, events)
% function write_events_tcga(fn_out, strains, events)

    if isempty(events),
        fprintf('No events present.\n');
        return
    end;

    event_type = events(1).event_type;

    fprintf('writing %s events in tcga format to %s\n', events(1).event_type, fn_out) ;

    fd = fopen(fn_out, 'w+') ;
    fprintf(fd, 'gene\teventtype\tcoordinates') ;
    fprintf(fd, '\t%s', strains{:});
    fprintf(fd, '\n') ;
    for i = 1:length(events),
        fprintf(fd, '%s\t%s\t%s:', events(i).gene_name{1}, event_type, events(i).chr);
        if strcmp(event_type, 'intron_reten'),
            fprintf(fd, ':%i-%i:%i-%i', ...
                events(i).exon1_col(1), events(i).exon1_col(2), ...
                events(i).exon2_col(1), events(i).exon2_col(2)) ;
        elseif strcmp(event_type, 'alt_3prime') || strcmp(event_type, 'alt_5prime'),
            fprintf(fd, ':%i-%i:%i-%i:%i-%i', ... 
                events(i).exon_const_col(1), events(i).exon_const_col(2), ...
                events(i).exon_alt1_col(1), events(i).exon_alt1_col(2), ...
                events(i).exon_alt2_col(1), events(i).exon_alt2_col(2)) ;
        elseif strcmp(event_type, 'exon_skip'),
            fprintf(fd, ':%i-%i:%i-%i:%i-%i', ... 
                events(i).exon_pre_col(1), events(i).exon_pre_col(2), ...
                events(i).exon_col(1), events(i).exon_col(2), ...
                events(i).exon_aft_col(1), events(i).exon_aft_col(2)) ;
        elseif strcmp(event_type, 'mult_exon_skip'),
            fprintf(fd, ':%i-%i', events(i).exon_pre_col(1), events(i).exon_pre_col(2));
            for j = 1:2:size(events(i).exons, 2)
                fprintf(fd, ':%i-%i:', events(i).exons_col(j), events(i).exons_col(j+1));
            end;
            fprintf(fd, ':%i-%i', events(i).exon_aft_col(1), events(i).exon_aft_col(2)) ;
        end;
        %%% TODO compute Splice Index !!!!
        for j = 1:length(strains),
            if events(i).info(j).valid,
                if (strcmp(events(i).event_type, 'alt_3prime') || strcmp(events(i).event_type, 'alt_5prime')),
                    if (events(i).intron1_col(2) - events(i).intron1_col(1)) < (events(i).intron2_col(2) - events(i).intron2_col(1)),
                        num = events(i).info(j).intron1_conf;
                    else
                        num = events(i).info(j).intron1_conf;
                    end;
                    denom = events(i).info(j).intron1_conf + events(i).info(j).intron2_conf;
                    confirmation = denom;
                elseif strcmp(events(i).event_type, 'exon_skip'),
                    num = events(i).info(j).exon_pre_exon_conf + events(i).info(j).exon_exon_aft_conf;
                    denom = events(i).info(j).exon_pre_exon_conf + events(i).info(j).exon_exon_aft_conf + (2 * events(i).info(j).exon_pre_exon_aft_conf);
                    confirmation = events(i).info(j).exon_pre_exon_conf + events(i).info(j).exon_exon_aft_conf + events(i).info(j).exon_pre_exon_aft_conf;
                elseif strcmp(events(i).event_type,  'mult_exon_skip'),
                    num = events(i).info(j).exon_pre_exon_conf + events(i).info(j).sum_inner_exon_conf + events(i).info(j).exon_exon_aft_conf;
                    denom = events(i).info(j).exon_pre_exon_conf + events(i).info(j).sum_inner_exon_conf + events(i).info(j).exon_exon_aft_conf + ((2 + events(i).info(j).num_inner_exon) * events(i).info(j).exon_pre_exon_aft_conf);
                    confirmation = events(i).info(j).exon_pre_exon_conf + events(i).info(j).sum_inner_exon_conf + events(i).info(j).exon_exon_aft_conf + events(i).info(j).exon_pre_exon_aft_conf;
                elseif strcmp(events(i).event_type,  'intron_retention'),
                    num = events(i).info(j).intron_conf;
                    denom = 1;
                    confirmation = num;
                else
                    error(sprintf('Unknown event type: %s\n', events(i).event_type));
                end;

                if confirmation < 10,
                    fprintf(fd, '\tNA') ;
                else
                    fprintf(fd, '\t%1.1f', num / denom);
                end;
            else
                fprintf(fd, '\tNA') ;
            end ;
        end ;
        fprintf(fd, '\n') ;
    end ;
    fclose(fd) ;
