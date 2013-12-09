function write_events_txt(fn_out, strains, events)
% function write_events_txt(fn_out, strains, events)

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
                fprintf(fd, ':%i-%:', events(i).exons_col(j), events(i).exons_col(j+1));
            end;
            fprintf(fd, ':%i-%i', events(i).exon_aft_col(1), events(i).exon_aft_col(2)) ;
        end;
        %%% TODO compute Splice Index !!!!
        for j = 1:length(strains),
            if events(i).info(j).valid,
                fprintf(fd, '\t%.1f', events(i).info(j).intron_conf) ;
            else
                fprintf(fd, '\tNA') ;
            end ;
        end ;
        fprintf(fd, '\n') ;
    end ;
    fclose(fd) ;
