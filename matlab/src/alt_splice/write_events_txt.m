function write_events_txt(fn_out_txt, strains, events, fn_counts, event_idx, anno_fn)
% function write_events_txt(fn_out_txt, strains, events, fn_counts, event_idx, anno_fn)
    
    if isempty(events),
        fprintf('No events present.\n');
        return
    end;

    if nargin < 5,
        event_idx = 1:length(events);
    end;

    anno_names = '';
    if nargin > 5 && ~isempty(anno_fn),
       load(anno_fn);
       anno = genes;
       [anno_names s_idx] = sort({genes.name});
       genes = genes(s_idx);
       clear genes
    end;

    fprintf('writing %s events in flat txt format to %s\n', events(1).event_type, fn_out_txt);

    fd = fopen(fn_out_txt, 'w+') ;
    if nargin > 5 && ~isempty(anno_names),
        gene_header = sprintf('\tgene_start\tgene_end');
    else
        gene_header = '';
    end;

    if strcmp(events(1).event_type, 'exon_skip'),
        fprintf(fd, 'contig\tstrand\tevent_id\tgene_name%s\texon_pre_start\texon_pre_end\texon_start\texon_end\texon_aft_start\texon_aft_end', gene_header) ;
        tmp = reshape(repmat(strains, 6, 1), 1, length(strains) * 6);
        fprintf(fd, '\t%s:exon_pre_cov\t%s:exon_cov\t%s:exon_aft_cov\t%s:intron_pre_conf\t%s:intron_aft_conf\t%s:intron_skip_conf', tmp{:}) ;
    elseif strcmp(events(1).event_type, 'alt_3prime') || strcmp(events(1).event_type, 'alt_5prime'),
        fprintf(fd, 'contig\tstrand\tevent_id\tgene_name%s\texon_const_start\texon_const_end\texon_alt1_start\texon_alt1_end\texon_alt2_start\texon_alt2_end', gene_header) ;
        tmp = reshape(repmat(strains, 4, 1), 1, length(strains) * 4);
        fprintf(fd, '\t%s:exon_diff_cov\t%s:exon_const_cov\t%s:intron1_conf\t%s:intron2_conf', tmp{:}) ;
    elseif strcmp(events(1).event_type, 'intron_retention'),
        fprintf(fd, 'contig\tstrand\tevent_id\tgene_name%s\texon1_start\texon1_end\tintron_start\tintron_end\texon2_start\texon2_end', gene_header) ;
        tmp = reshape(repmat(strains, 4, 1), 1, length(strains) * 4);
        fprintf(fd, '\t%s:exon1_cov\t%s:intron_cov\t%s:exon2_cov\t%s:intron_conf', tmp{:}) ;
    elseif strcmp(events(1).event_type, 'mult_exon_skip'),
        fprintf(fd, 'contig\tstrand\tevent_id\tgene_name%s\texon_pre_start\texon_pre_end\texon_starts\texon_ends\texon_aft_start\texon_aft_end', gene_header) ;
        tmp = reshape(repmat(strains, 8, 1), 1, length(strains) * 8);
        fprintf(fd, '\t%s:exon_pre_cov\t%s:exon_cov\t%s:exon_aft_cov\t%s:intron_pre_conf\t%s:intron_inner_conf\t%s:exon_inner_count\t%s:intron_aft_conf\t%s:intron_skip_conf', tmp{:}) ;
    elseif strcmp(events(1).event_type, 'mutex_exons'),
        fprintf(fd, 'contig\tstrand\tevent_id\tgene_name%s\texon_pre_start\texon_pre_end\texon1_start\texon1_end\texon2_start\texon2_end\texon_aft_start\texon_aft_end', gene_header) ;
        tmp = reshape(repmat(strains, 8, 1), 1, length(strains) * 8);
        fprintf(fd, '\t%s:exon_pre_cov\t%s:exon1_cov\t%s:exon2_cov\t%s:exon_aft_cov\t%s:pre_exon1_conf\t%s:pre_exon2_conf\t%s:exon1_aft_conf\t%s:exon2_aft_conf', tmp{:}) ;
    else
        error(sprintf('Unknown event type: %s\n', events(1).event_type));
    end;
    fprintf(fd, '\n') ;

    %%% prepare chunked loading of count hdf5
    size_info = h5info(fn_counts, '/event_counts');
    events_size = size_info.Dataspace.Size;
    plist = 'H5P_DEFAULT';
    fid = H5F.open(fn_counts);
    dset_id = H5D.open(fid,'/event_counts');
    dims = fliplr([events_size(1:2) 1]);
    file_space_id = H5D.get_space(dset_id);
    mem_space_id = H5S.create_simple(3, dims, []);
    for i = event_idx,
        %%% load current event-slice from count files
        H5S.select_hyperslab(file_space_id, 'H5S_SELECT_SET', fliplr([0 0 i - 1]),[], [], dims);
        counts = H5D.read(dset_id, 'H5ML_DEFAULT', mem_space_id, file_space_id, plist);

        fprintf(fd, '%s\t%c\t%s_%i\t%s', events(i).chr, events(i).strand, events(i).event_type, events(i).id, events(i).gene_name{1});
        if nargin > 5 && ~isempty(anno_names),
            a_idx = strmatch(events(i).gene_name{1}, anno_names, 'exact');
            assert(~isempty(a_idx));
            fprintf(fd, '\t%i\t%i', anno(a_idx).start, anno(a_idx).stop);
        end;
        ev = events(i);
        if strcmp(events(i).event_type, 'exon_skip'),
            fprintf(fd, '\t%i\t%i\t%i\t%i\t%i\t%i', ev.exon_pre_col(1), ev.exon_pre_col(2), ev.exon_col(1), ev.exon_col(2), ...
                    ev.exon_aft_col(1), ev.exon_aft_col(2)) ;
            for j = 1:length(strains),
                if counts(j, 1),
                    fprintf(fd, '\t%.1f\t%.1f\t%.1f\t%i\t%i\t%i', counts(j, 3), counts(j, 2), counts(j, 4), counts(j, 5), counts(j, 6), counts(j, 7)) ;
                else
                    fprintf(fd, '\t-1\t-1\t-1\t-1') ;
                end ;
            end ;
        elseif strcmp(ev.event_type, 'intron_retention'),
            fprintf(fd, '\t%i\t%i\t%i\t%i\t%i\t%i', ev.exon1_col(1), ev.exon1_col(2), ev.intron_col(1), ev.intron_col(2), ...
                    ev.exon2_col(1), ev.exon2_col(2)) ;
            for j = 1:length(strains),
                if counts(j, 1),
                    fprintf(fd, '\t%.1f\t%.1f\t%.1f\t%i', counts(j, 3), counts(j, 2), counts(j, 4), counts(j, 5)) ;
                else
                    fprintf(fd, '\t-1\t-1\t-1') ;
                end ;
            end ;
        elseif strcmp(ev.event_type, 'alt_3prime') || strcmp(ev.event_type, 'alt_5prime'),
            fprintf(fd, '\t%i\t%i\t%i\t%i\t%i\t%i', ev.exon_const_col(1), ev.exon_const_col(2), ev.exon_alt1_col(1), ...
                    ev.exon_alt1_col(2), ev.exon_alt2_col(1), ev.exon_alt2_col(2)) ;
            for j = 1:length(strains),
                if counts(j, 1),
                    fprintf(fd, '\t%1.1f\t%1.1f\t%i\t%i', counts(j, 2), counts(j, 3), counts(j, 4), counts(j, 5)) ;
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
            fprintf(fd, '\t%s\t%s\t%i\t%i', starts, ends, ev.exon_aft_col(1), ev.exon_aft_col(2)) ;
            for j = 1:length(strains),
                if counts(j, 1),
                    fprintf(fd, '\t%.1f\t%.1f\t%.1f\t%i\t%i\t%i\t%i\t%i', counts(j, 2), counts(j, 3), counts(j, 4), counts(j, 5), ...
                            counts(j, 8), counts(j, 9), counts(j, 6), counts(j, 7)) ;
                else
                    fprintf(fd, '\t-1\t-1\t-1\t-1\t-1\t-1\t-1') ;
                end ;
            end ;
        elseif strcmp(events(i).event_type, 'mutex_exons'),
            fprintf(fd, '\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i', ev.exon_pre_col(1), ev.exon_pre_col(2), ev.exon1_col(1), ev.exon1_col(2), ...
                    ev.exon2_col(1), ev.exon2_col(2), ev.exon_aft_col(1), ev.exon_aft_col(2)) ;
            for j = 1:length(strains),
                if counts(j, 1),
                    fprintf(fd, '\t%.1f\t%.1f\t%.1f\t%.1f\t%i\t%i\t%i\t%i', counts(j, 2), counts(j, 3), counts(j, 4), counts(j, 5), counts(j, 6), counts(j, 7), counts(j, 8), counts(j, 9)) ;
                else
                    fprintf(fd, '\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1') ;
                end ;
            end ;
        end;
        fprintf(fd, '\n') ;
    end ;
    fclose(fd) ;
    %%% close HDF5 files
    H5D.close(dset_id);
    H5F.close(fid);

