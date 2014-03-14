function events = post_process_event_struct(events)

    if isempty([events.event_type]),
        return
    end;

    %%% filter out invalid coordinate projections
    idx_valid_col = zeros(1, length(events)) ;
    for i = 1:length(events),
        if strcmp(events(i).event_type, 'intron_retention'),
            idx_valid_col(i) = all(events(i).intron_col_pos > 0);
        elseif strcmp(events(i).event_type, 'alt_3prime') || strcmp(events(i).event_type, 'alt_5prime'),
            idx_valid_col(i) = all([events(i).exon_const_col events(i).exon_alt1_col events(i).exon_alt2_col] > 0);
        elseif strcmp(events(i).event_type, 'mult_exon_skip'),
            idx_valid_col(i) = all(events(i).exons_col_pos > 0) ;
        else
            idx_valid_col(i) = all(events(i).exon_col_pos > 0) ;
        end;
    end ;
    events = events(idx_valid_col~=0) ;
   
    %%% sort exon skip events by all coordinates
    events = sort_events_full(events) ;
    
    %%% make exon skip events unique by strain
    fprintf('\nMake %s events unique by strain\n', events(1).event_type);
    events = make_unique_by_strain(events);

    %%% sort exon skip events by event coordinates
    events = sort_events_by_event(events) ;
    
    %%% make exon skip events unique by strain
    fprintf('\nMake %s events unique by event\n', events(1).event_type);
    events = make_unique_by_event(events);

    %%% count detected strains
    for i = 1:length(events),
        events(i).num_detected = length(events(i).strain) ;
        events(i).id = i;
    end ;

 
    return
