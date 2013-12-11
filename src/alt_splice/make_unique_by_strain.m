function event_list = make_unique_by_strain(event_list)

    rm_idx=[] ;
    field_names = {'exon', 'exons', 'exon1', 'exon2', 'exon_alt1', 'exon_alt2', 'exon_const', 'intron', 'intron1', 'intron2', ...
                   'exon_pre', 'exon_aft'} ;

    for i = 2:length(event_list),
        if mod(i, 1000) == 0
            fprintf(1, '.');
            if mod(i, 10000) == 0
                fprintf(1, '%i\n', i);
            end;
        end;
        if (isfield(event_list, 'exon_col') && isfield(event_list, 'exon_pre_col') && isfield(event_list, 'exon_aft_col') && ...
            isequal([event_list(i-1).exon_pre_col event_list(i-1).exon_col event_list(i-1).exon_aft_col], ...
                    [event_list(i).exon_pre_col event_list(i).exon_col event_list(i).exon_aft_col])) || ...
           (isfield(event_list, 'exons_col') && isfield(event_list, 'exon_pre_col') && isfield(event_list, 'exon_aft_col') && ...
            isequal([event_list(i-1).exon_pre_col event_list(i-1).exons_col event_list(i-1).exon_aft_col], ...
                    [event_list(i).exon_pre_col event_list(i).exons_col event_list(i).exon_aft_col])) || ...
           (isfield(event_list, 'intron_col') && isfield(event_list, 'exon1_col') && isfield(event_list, 'exon2_col') && ...
            isequal([event_list(i-1).exon1_col event_list(i-1).intron_col event_list(i-1).exon2_col], ...
                    [event_list(i).exon1_col event_list(i).intron_col event_list(i).exon2_col])) || ...
           (isfield(event_list, 'intron1_col') && isfield(event_list, 'intron2_col') && isfield(event_list, 'exon_alt1_col') && isfield(event_list, 'exon_alt2_col') && ...
            isequal([event_list(i-1).intron1_col event_list(i-1).exon_alt1_col event_list(i-1).intron2_col event_list(i-1).exon_alt2_col], ...
                    [event_list(i).intron1_col event_list(i).exon_alt1_col event_list(i).intron2_col event_list(i).exon_alt2_col])),

            %%% assertion that we did everything right
            assert(event_list(i-1).chr_num == event_list(i).chr_num) ;
            assert(event_list(i-1).strand == event_list(i).strand) ;
            if isfield(event_list, 'exon_col_pos'),
                assert(isequal(event_list(i-1).exon_col_pos, event_list(i).exon_col_pos)) ;
            end;
            if isfield(event_list, 'exons_col_pos'),
                assert(isequal(event_list(i-1).exons_col_pos, event_list(i).exons_col_pos)) ;
            end;
            if isfield(event_list, 'intron_col_pos'),
                assert(isequal(event_list(i-1).intron_col_pos, event_list(i).intron_col_pos)) ;
            end;
            
            assert(length(event_list(i).strain) == 1);
            idx = strmatch(event_list(i).strain{1}, event_list(i-1).strain, 'exact') ;
            if ~isempty(idx),
                assert(length(idx) == 1) ;
                if isfield(event_list(i), 'exon'), assert(isequal(event_list(i).exon, event_list(i-1).exon(idx,:))); end;
                if isfield(event_list(i), 'exons'), assert(isequal(event_list(i).exons, event_list(i-1).exons(idx,:))); end;
                if isfield(event_list(i), 'intron'), assert(isequal(event_list(i).intron, event_list(i-1).intron(idx,:))); end;
                if isfield(event_list(i), 'intron1') && isfield(event_list(i), 'intron2'),
                    assert(isequal([event_list(i).intron1; event_list(i).intron2], [event_list(i-1).intron1(idx,:); event_list(i-1).intron2(idx, :)]));
                end;
                if isfield(event_list(i), 'gene_name') && isempty(strmatch(event_list(i).gene_name{1}, event_list(i-1).gene_name, 'exact')),
                    event_list(i-1).gene_name = {event_list(i-1).gene_name{:}, event_list(i).gene_name{:}};
                end;
                event_list(i) = event_list(i-1) ;
            else 
                event_list(i).strain = [event_list(i-1).strain event_list(i).strain] ;
                assert(isequal(sort(event_list(i).strain), sort(unique(event_list(i).strain)))) ; %%% TODO
                
                for f_idx = 1:length(field_names),
                    fn = field_names{f_idx};
                    eval(sprintf('if isfield(event_list, ''%s''), event_list(i).%s = [event_list(i-1).%s; event_list(i).%s]; end;', fn, fn, fn, fn));
                end;
                if isfield(event_list(i), 'gene_name') && isempty(strmatch(event_list(i).gene_name{1}, event_list(i-1).gene_name, 'exact')),
                    event_list(i).gene_name = {event_list(i-1).gene_name{:}, event_list(i).gene_name{:}};
                end;
            end ;
            rm_idx(end+1) = i-1 ;
        end ;
    end;
    fprintf('events dropped: %i\n', length(rm_idx));
    event_list(rm_idx) = [];
