function event_list = make_unique_by_event(event_list)
% function event_list = make_unique_by_event(event_list)
%
% This script removes all events that share the sam alternative evnt coordinates
% but differ in the flanking size. The longest of several equal events is kept.

    rm_idx=[] ;
    field_names = {'exon', 'exons', 'exon1', 'exon2', 'exon_alt1', 'exon_alt2', 'exon_const', 'intron', 'intron1', 'intron2', ...
                   'exon_pre', 'exon_aft'} ;

    last_kept = 1;
    for i = 2:length(event_list),
        if mod(i, 1000) == 0
            fprintf(1, '.');
            if mod(i, 10000) == 0
                fprintf(1, '%i\n', i);
            end;
        end;
        if (strcmp(event_list(i).event_type, 'exon_skip') && ...
            isequal([event_list(last_kept).exon_pre_col(2) event_list(last_kept).exon_col event_list(last_kept).exon_aft_col(1)], ...
                    [event_list(i).exon_pre_col(2) event_list(i).exon_col event_list(i).exon_aft_col(1)])) || ...
           (strcmp(event_list(i).event_type, 'mult_exon_skip') && ...
            isequal([event_list(last_kept).exon_pre_col(2) event_list(last_kept).exons_col event_list(last_kept).exon_aft_col(1)], ...
                    [event_list(i).exon_pre_col(2) event_list(i).exons_col event_list(i).exon_aft_col(1)])) || ...
           (strcmp(event_list(i).event_type, 'intron_retention') && ...
            isequal(event_list(last_kept).intron_col, event_list(i).intron_col)) || ...
           ((strcmp(event_list(i).event_type, 'alt_5prime') || strcmp(event_list(i).event_type, 'alt_3prime')) && ...
            isequal(sort(unique([event_list(last_kept).intron1_col event_list(last_kept).intron2_col])), ...
                    sort(unique([event_list(i).intron1_col event_list(i).intron2_col])))),

            %%% assertion that we did everything right
            assert(event_list(last_kept).chr_num == event_list(i).chr_num) ;
            assert(event_list(last_kept).strand == event_list(i).strand) ;
            if isfield(event_list, 'exon_col_pos'),
                assert(isequal(event_list(last_kept).exon_col_pos, event_list(i).exon_col_pos)) ;
            end;
            if isfield(event_list, 'exons_col_pos'),
                assert(isequal(event_list(last_kept).exons_col_pos, event_list(i).exons_col_pos)) ;
            end;
            if isfield(event_list, 'intron_col_pos'),
                assert(isequal(event_list(last_kept).intron_col_pos, event_list(i).intron_col_pos)) ;
            end;

            %%% match and merge strains
            idx = setdiff(event_list(i).strain, event_list(last_kept).strain);
            
            %%% check, which event is longer -> keep longer event
            if strcmp(event_list(i).event_type, 'exon_skip') || strcmp(event_list(i).event_type, 'mult_exon_skip'),
                len1 = event_list(last_kept).exon_aft(2) - event_list(last_kept).exon_pre(1) + 1;   
                len2 = event_list(i).exon_aft(2) - event_list(i).exon_pre(1) + 1;   
            elseif strcmp(event_list(i).event_type, 'intron_retention'),
                len1 = event_list(last_kept).exon2(2) - event_list(last_kept).exon1(1) + 1;
                len2 = event_list(i).exon2(2) - event_list(i).exon1(1) + 1;
            elseif strcmp(event_list(i).event_type, 'alt_5prime') || strcmp(event_list(i).event_type, 'alt_3prime'),
                if event_list(i).exon_const(1) < event_list(i).exon_alt1(1),
                    len1 = max(event_list(last_kept).exon_alt1(2), event_list(last_kept).exon_alt2(2)) - event_list(last_kept).exon_const(1) + 1; 
                    len2 = max(event_list(i).exon_alt1(2), event_list(i).exon_alt2(2)) - event_list(i).exon_const(1) + 1; 
                else
                    len1 = event_list(last_kept).exon_const(2) - min(event_list(last_kept).exon_alt1(1), event_list(last_kept).exon_alt2(1)) + 1; 
                    len2 = event_list(i).exon_const(2) - min(event_list(i).exon_alt1(1), event_list(i).exon_alt2(1)) + 1; 
                end;
            end;

            if len1 > len2,
                keep_idx = last_kept;
                not_keep_idx = i;
            else
                keep_idx = i;
                not_keep_idx = last_kept;
            end;

            %%% check if we would loose strains 
            [dummy, idx] = setdiff(event_list(not_keep_idx).strain, event_list(keep_idx).strain);
            if ~isempty(idx),
                event_list(keep_idx).strain = [event_list(keep_idx).strain; event_list(not_keep_idx).strain(idx)];
                for f_idx = 1:length(field_names),
                    fn = [field_names{f_idx} '(idx, :)'];
                    %eval(sprintf('if isfield(event_list, ''%s''), event_list(keep_idx).%s = [event_list(keep_idx).%s; event_list(not_keep_idx).%s]; end;', fn, fn, fn, fn));
                    if isfield(event_list, fn), 
                        event_list(keep_idx).(fn) = [event_list(keep_idx).(fn); event_list(not_keep_idx).(fn)]; 
                    end;
                end;
                if isfield(event_list(i), 'gene_name'),
                    event_list(keep_idx).gene_name = union(event_list(keep_idx).gene_name, event_list(not_keep_idx).gene_name);
                end;
            end;
            rm_idx(end + 1) = not_keep_idx;
            last_kept = keep_idx;
        else
            last_kept = i;
        end ;
    end;
    fprintf('events dropped: %i\n', length(rm_idx));
    event_list(rm_idx) = [];
