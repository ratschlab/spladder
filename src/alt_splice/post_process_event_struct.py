def post_process_event_struct(events):
    # events = post_process_event_struct(events)

    ### filter out invalid coordinate projections
    idx_valid_col = sp.zeros((events.shape[0],))
    for i in range(events.shape[0]):
        idx_valid[i] = sp.all(events[i].get_coords(trafo=True) > 0)
    events = events[sp.where(idx_valid_col)[0]]
   
    ### sort exon skip events by all coordinates
    events = sort_events_full(events) 
    
    ### make exon skip events unique by strain
    print '\nMake %s events unique by strain', events[0].event_type
    events = make_unique_by_strain(events)

    ### sort exon skip events by event coordinates
    events = sort_events_by_event(events) 
    
    ### make exon skip events unique by strain
    print '\nMake %s events unique by event' % events[0].event_type
    events = make_unique_by_event(events)

    ### count detected strains
    for i in range(events.shape[0]):
        events[i].num_detected = len(events[i].strain)
        events[i].id = i

    return events
