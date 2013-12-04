function write_events_gff3(fn_out_gff3, events)
% function write_events_gff3(fn_out_gff3, events)
    
    if isempty(events),
        fprintf('No events present.\n');
        return
    end;

    fprintf('writing %s events in gff3 format to %s\n', events(1).event_type, fn_out_gff3) ;

    fd_out = fopen(fn_out_gff3, 'w+') ; 
    fprintf(fd_out, '##gff-version 3\n') ;

    %%% load gene structure
    for i = 1:size(events, 2),

        gene_name = events(i).gene_name{1};
        %start_pos = events(i).start;
        %stop_pos = events(i).stop;
        if strcmp(events(i).event_type, 'exon_skip'),
            start_pos = events(i).exon_pre(1);
            stop_pos = events(i).exon_aft(2);
        elseif strcmp(events(i).event_type, 'intron_retention'),
            start_pos = events(i).exon1(1);
            stop_pos = events(i).exon2(2);
        elseif strcmp(events(i).event_type, 'alt_3prime') || strcmp(events(i).event_type, 'alt_5prime'),
            start_pos = min([events(i).exon_const(1), events(i).exon_alt1(1), events(i).exons_alt2(1)]);
            stop_pos = min([events(i).exon_const(2), events(i).exon_alt1(2), events(i).exons_alt2(2)]);
        end;

        %%% get order of isoforms o_idx(1) -> iso1 and o_idx(2) -> iso2
        o_idx = [1 2];
        %%% for alt_prime the longer isoform (shorter intron) is iso1
        if (strcmp(events(i).event_type, 'alt_5prime') ||  strcmp(events(i).event_type, 'alt_3prime')) && ...
           (events(1, i).exons{1}(2, 1) - events(1, i).exons{1}(1, 2) > events(1, i).exons{2}(2, 1) - events(1, i).exons{2}(1, 2)),
            o_idx = [2 1];
        end;
        name = sprintf('%s.%i', events(i).event_type, i);

        fprintf(fd_out, '%s\t%s\tgene\t%i\t%i\t.\t%c\t.\tID=%s;GeneName="%s"\n', events(1, i).chr, events(1, i).event_type, start_pos,  ...
                stop_pos, events(1, i).strand, name, events(1, i).gene_name{1}) ;
        fprintf(fd_out, '%s\t%s\tmRNA\t%i\t%i\t.\t%c\t.\tID=%s_iso1;Parent=%s;GeneName="%s"\n', events(1, i).chr, events(1, i).event_type, ...
                start_pos, stop_pos, events(1, i).strand, name, name, events(1, i).gene_name{1}) ;
        fprintf(fd_out, '%s\t%s\texon\t%i\t%i\t.\t%c\t.\tParent=%s_iso1\n', events(1, i).chr, events(1, i).event_type, ...
                events(1, i).exons{o_idx(1)}(1, 1), events(1, i).exons{o_idx(1)}(1, 2), events(1, i).strand, name) ;
        fprintf(fd_out, '%s\t%s\texon\t%i\t%i\t.\t%c\t.\tParent=%s_iso1\n', events(1, i).chr, events(1, i).event_type, ...
                events(1, i).exons{o_idx(1)}(2, 1), events(1, i).exons{o_idx(1)}(2, 2), events(1, i).strand, name) ;
        if strcmp(events(1, i).event_type, 'exon_skip'),
            fprintf(fd_out, '%s\t%s\texon\t%i\t%i\t.\t%c\t.\tParent=%s_iso1\n', events(1, i).chr, events(1, i).event_type, ...
                    events(1, i).exons{o_idx(1)}(3, 1), events(1, i).exons{o_idx(1)}(3, 2), events(1, i).strand, name) ;
        end;
        fprintf(fd_out, '%s\t%s\tmRNA\t%i\t%i\t.\t%c\t.\tID=%s_iso2;Parent=%s;GeneName="%s"\n', events(1, i).chr, events(1, i).event_type, ...
                start_pos, stop_pos, events(1, i).strand, name, name, events(1, i).gene_name{1}) ;
        fprintf(fd_out, '%s\t%s\texon\t%i\t%i\t.\t%c\t.\tParent=%s_iso2\n', events(1, i).chr, events(1, i).event_type, ...
                events(1, i).exons{o_idx(2)}(1, 1), events(1, i).exons{o_idx(2)}(1, 2), events(1, i).strand, name) ;
        if ~strcmp(events(1, i).event_type, 'intron_ret') && ~strcmp(events(1, i).event_type, 'intron_retention'),
            fprintf(fd_out, '%s\t%s\texon\t%i\t%i\t.\t%c\t.\tParent=%s_iso2\n', events(1, i).chr, events(1, i).event_type, ...
                    events(1, i).exons{o_idx(2)}(2, 1), events(1, i).exons{o_idx(2)}(2, 2), events(1, i).strand, name) ;
        end;
    end ;
    fclose(fd_out) ;
