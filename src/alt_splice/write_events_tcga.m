function write_events_tcga(fn_out, strains, events, fn_counts, event_idx)
% function write_events_tcga(fn_out, strains, events, fn_counts, event_idx)

    if isempty(events),
        fprintf('No events present.\n');
        return
    end;

    if nargin <5,
        event_idx = 1:length(events);
    end;

    event_type = events(1).event_type;

    %%% prepare chunked loading of count hdf5
    size_info = h5info(fn_counts, '/event_counts');
    events_size = size_info.Dataspace.Size;
    plist = 'H5P_DEFAULT';
    fid = H5F.open(fn_counts);
    dset_id = H5D.open(fid,'/event_counts');
    dims = fliplr([events_size(1:2) 1]);
    file_space_id = H5D.get_space(dset_id);
    mem_space_id = H5S.create_simple(3, dims, []);
 
    fprintf('writing %s events in tcga format to %s\n', events(1).event_type, fn_out) ;

    fd = fopen(fn_out, 'w+') ;
    fprintf(fd, 'gene\teventtype\tcoordinates') ;
    fprintf(fd, '\t%s', strains{:});
    fprintf(fd, '\n') ;
    for i = event_idx,
        %%% load current event-slice from count files
        H5S.select_hyperslab(file_space_id, 'H5S_SELECT_SET', fliplr([0 0 i - 1]),[], [], dims);
        counts = H5D.read(dset_id, 'H5ML_DEFAULT', mem_space_id, file_space_id, plist);

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
            if counts(j, 1),
                if (strcmp(events(i).event_type, 'alt_3prime') || strcmp(events(i).event_type, 'alt_5prime')),
                    if (events(i).intron1_col(2) - events(i).intron1_col(1)) < (events(i).intron2_col(2) - events(i).intron2_col(1)),
                        num = counts(j, 4);
                    else
                        num = counts(j, 5);
                    end;
                    denom = counts(j, 4) + counts(j, 5);
                    confirmation = denom;
                elseif strcmp(events(i).event_type, 'exon_skip'),
                    num = counts(j, 5) + counts(j, 6);
                    denom = counts(j, 5) + counts(j, 6) + (2 * counts(j, 7));
                    confirmation = counts(j, 5) + counts(j, 6) + counts(j, 7);
                elseif strcmp(events(i).event_type,  'mult_exon_skip'),
                    num = counts(j, 5) + counts(j, 8) + counts(j, 6);
                    denom = counts(j, 5) + counts(j, 8) + counts(j, 6) + ((2 + counts(j, 9)) * counts(j, 7));
                    confirmation = counts(j, 5) + counts(j, 8) + counts(j, 6) + counts(j, 7);
                elseif strcmp(events(i).event_type,  'intron_retention'),
                    num = counts(j, 5);
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
