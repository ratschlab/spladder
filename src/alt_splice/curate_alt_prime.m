function event_list = curate_alt_prime(event_list, CFG)

    if isempty(event_list) || isempty([event_list.event_type]),
        return;
    end;

    rm_idx = [];
    corr_count = 0;

    for i = 1:length(event_list)
        
        %%% check if alt exons overlap
        if event_list(i).exon_alt1_col(1) > event_list(i).exon_alt2_col(2) || event_list(i).exon_alt1_col(2) < event_list(i).exon_alt2_col(1)
            rm_idx(end + 1) = i;
            continue;
        end;

        %%% check if we have introns of zero length
        if size(unique([event_list(i).intron1_col event_list(i).intron2_col]), 2) ~= 3,
            rm_idx(end + 1) = i;
            continue;
        end;
        
        if event_list(i).exon_alt1_col(2) < event_list(i).exon_const_col(1),
            if event_list(i).exon_alt1_col(1) < event_list(i).exon_alt2_col(1),
                event_list(i).exon_alt1_col(1) = event_list(i).exon_alt2_col(1);
                event_list(i).exon_alt1(:, 1) = event_list(i).exon_alt2(:, 1);
                corr_count = corr_count + 1;
            elseif event_list(i).exon_alt1_col(1) > event_list(i).exon_alt2_col(1),
                event_list(i).exon_alt2_col(1) = event_list(i).exon_alt1_col(1);
                event_list(i).exon_alt2(:, 1) = event_list(i).exon_alt1(:, 1);
                corr_count = corr_count + 1;
            end;
        else
            if event_list(i).exon_alt1_col(2) < event_list(i).exon_alt2_col(2),
                event_list(i).exon_alt2_col(2) = event_list(i).exon_alt1_col(2);
                event_list(i).exon_alt2(:, 2) = event_list(i).exon_alt1(:, 2);
                corr_count = corr_count + 1;
            elseif event_list(i).exon_alt2_col(2) < event_list(i).exon_alt1_col(2),
                event_list(i).exon_alt1_col(2) = event_list(i).exon_alt2_col(2);
                event_list(i).exon_alt1(:, 2) = event_list(i).exon_alt2(:, 2);
                corr_count = corr_count + 1;
            end;
        end;
    end;

    %%% remove events with non-overlapping alt_exons
    if ~isempty(rm_idx),
        event_list(rm_idx) = [];
    end;

    fprintf(1, 'Corrected %i events\n', corr_count);
    fprintf(1, 'Removed %i events\n', length(rm_idx));
