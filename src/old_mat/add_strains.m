function [events, event_strains] = add_strains(events, CFG)

fieldnames = {'intron', 'exon1', 'exon2', 'exon_alt1', 'exon_alt2', 'exon_const', 'intron1', 'intron2', 'exon', 'exon_pre', 'exon_aft', 'exons'};

curr_fields = {};
for f_idx = 1:length(fieldnames),
    field = fieldnames{f_idx};
    if isfield(events, [field '_col']) 
        curr_fields{end + 1} = field;
    end;
end;

event_strains = {};
for i = 1:length(events),
    if mod(i,100) == 0,
        fprintf(1, '.');
        if mod(i, 1000) == 0,
            fprintf(1, '%i (%i)\n', i, length(events));
        end ;
    end ;

    new_strains = setdiff(CFG.strains, events(i).strain) ;

    if ~isfield(events, 'detected') || isempty(events(i).detected),
        events(i).detected = ones(1, length(events(i).strain)) ;
    end ;

    if isfield(CFG, 'reference_strain'),
        for s_idx = 1:length(new_strains),
            strain = new_strains{s_idx} ;

            for f_idx = 1:length(curr_fields),
                field = curr_fields{f_idx};

                %%% conversion to ref_strain coordinates
                tmp = convert_strain_intervals(events(i).chr_num, events(i).([field '_col']), CFG.reference_strain, strain) ;
                if isempty(tmp),
                    tmp = convert_strain_pos(events(i).chr_num, events(i).([field '_col']), CFG.reference_strain, strain) ;
                end;
                if s_idx == 1,
                    events(i).(field) = repmat(events(i).(field), length(new_strains) + 1, 1);
                end;
                events(i).(field)(s_idx + 1, :) = tmp;
            end;
        end ;
        event_strains(i, :) = [events(i).strain new_strains];
        events(i).detected = [events(i).detected zeros(1, length(new_strains))]; 
        [tmp, idx1, idx2] = intersect(CFG.strains, event_strains(i, :)) ;
        assert(length(idx2) == length(event_strains(i, :))) ;
        assert(length(idx2) == length(CFG.strains)) ;
        if ~isequal(idx1,idx2),
            event_strains(i, idx1) = event_strains(i, idx2) ;
            for f_idx = 1:length(curr_fields),
                field = curr_fields{f_idx};
                events(i).(field)(idx1, :) = events(i).(field)(idx2, :) ;
            end;
            events(i).detected(idx1) = events(i).detected(idx2) ;
        end ;
        assert(isequal(CFG.strains, event_strains(i, :).strain)) ;
    else
        if i == 1,
            event_strains(i, :) = [events(i).strain new_strains];
            [tmp, idx1, idx2] = intersect(CFG.strains, event_strains(i, :)) ;
            assert(length(idx2) == length(event_strains(i, :))) ;
            assert(length(idx2) == length(CFG.strains)) ;
            if ~isequal(idx1,idx2),
                event_strains(i, idx1) = event_strains(i, idx2) ;
            end ;
            assert(isequal(CFG.strains, event_strains(i, :))) ;
        end;
        events(i).detected = [events(i).detected zeros(1, length(new_strains))]; 
        events(i).detected(idx1) = events(i).detected(idx2) ;
    end ;
end ;
