function write_events_gff3(fn_out_gff3, events, idx)
% function write_events_gff3(fn_out_gff3, events, idx)
    
    if isempty(events),
        fprintf('No events present.\n');
        return
    end;

    if nargin < 3,
        idx = 1:size(events, 2);
    end;

    fprintf('writing %s events in gff3 format to %s\n', events(1).event_type, fn_out_gff3) ;

    fd_out = fopen(fn_out_gff3, 'w+') ; 
    fprintf(fd_out, '##gff-version 3\n') ;

    %%% load gene structure
    for i = idx,

        ev = events(i);
        exons = {};

        gene_name = events(i).gene_name{1};
        if strcmp(ev.event_type, 'exon_skip'),
            exons{1} = [ev.exon_pre; ev.exon; ev.exon_aft];
            exons{2} = [ev.exon_pre; ev.exon_aft];
        elseif strcmp(ev.event_type, 'intron_retention'),
            exons{1} = [ev.exon1(1) ev.exon2(2)];
            exons{2} = [ev.exon1; ev.exon2];
        elseif strcmp(ev.event_type, 'alt_3prime') || strcmp(ev.event_type, 'alt_5prime'),
            exons{1} = sortrows([ev.exon_const; ev.exon_alt1]);
            exons{2} = sortrows([ev.exon_const; ev.exon_alt2]);
        elseif strcmp(ev.event_type, 'mult_exon_skip'),
            exons{1} = [ev.exon_pre; reshape(ev.exons, length(ev.exons) / 2, 2); ev.exon_aft];
            exons{2} = [ev.exon_pre; ev.exon_aft];
        end;

        start_pos = exons{1}(1, 1);
        stop_pos = exons{1}(end, end);

        %%% get order of isoforms o_idx(1) -> iso1 and o_idx(2) -> iso2
        o_idx = [1 2];
        %%% for alt_prime the longer isoform (shorter intron) is iso1
        if (strcmp(ev.event_type, 'alt_5prime') || strcmp(ev.event_type, 'alt_3prime')) && ...
           (exons{1}(2, 1) - exons{1}(1, 2) > exons{2}(2, 1) - exons{2}(1, 2)),
            o_idx = [2 1];
        end;
        name = sprintf('%s.%i', ev.event_type, ev.id);

        fprintf(fd_out, '%s\t%s\tgene\t%i\t%i\t.\t%c\t.\tID=%s;GeneName="%s"\n', ev.chr, ev.event_type, start_pos, stop_pos, ev.strand, name, ev.gene_name{1}) ;
        fprintf(fd_out, '%s\t%s\tmRNA\t%i\t%i\t.\t%c\t.\tID=%s_iso1;Parent=%s;GeneName="%s"\n', ev.chr, ev.event_type, ...
                start_pos, stop_pos, ev.strand, name, name, ev.gene_name{1}) ;
        fprintf(fd_out, '%s\t%s\texon\t%i\t%i\t.\t%c\t.\tParent=%s_iso1\n', ev.chr, ev.event_type, ...
                exons{o_idx(1)}(1, 1), exons{o_idx(1)}(1, 2), ev.strand, name) ;
        if ~strcmp(ev.event_type, 'intron_ret') && ~strcmp(ev.event_type, 'intron_retention'),
            fprintf(fd_out, '%s\t%s\texon\t%i\t%i\t.\t%c\t.\tParent=%s_iso1\n', ev.chr, ev.event_type, ...
                    exons{o_idx(1)}(2, 1), exons{o_idx(1)}(2, 2), ev.strand, name) ;
        end;
        if strcmp(ev.event_type, 'exon_skip'),
            fprintf(fd_out, '%s\t%s\texon\t%i\t%i\t.\t%c\t.\tParent=%s_iso1\n', ev.chr, ev.event_type, ...
                    exons{o_idx(1)}(3, 1), exons{o_idx(1)}(3, 2), events(1, i).strand, name) ;
        end;
        fprintf(fd_out, '%s\t%s\tmRNA\t%i\t%i\t.\t%c\t.\tID=%s_iso2;Parent=%s;GeneName="%s"\n', ev.chr, ev.event_type, ...
                start_pos, stop_pos, ev.strand, name, name, ev.gene_name{1}) ;
        fprintf(fd_out, '%s\t%s\texon\t%i\t%i\t.\t%c\t.\tParent=%s_iso2\n', ev.chr, ev.event_type, ...
                exons{o_idx(2)}(1, 1), exons{o_idx(2)}(1, 2), ev.strand, name) ;
        fprintf(fd_out, '%s\t%s\texon\t%i\t%i\t.\t%c\t.\tParent=%s_iso2\n', ev.chr, ev.event_type, ...
                exons{o_idx(2)}(2, 1), exons{o_idx(2)}(2, 2), ev.strand, name) ;
    end ;
    fclose(fd_out) ;
