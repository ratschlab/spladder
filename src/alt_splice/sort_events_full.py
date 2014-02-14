def sort_events_full(event_list):
# event_list = sort_events_full(event_list),

    if event_list.shape[0] == 0:
        return event_list

    coord_list = sp.array([x.get_coords() for x in event_list]) 
    chr_list = sp.array([x.chr_num for x in event_list])
    strand_list = sp.array([x.strand for x in event_list], dtype = 'double')
    sort_list = sp.c_[chr_list, strand_list, sort_list]
    tmp, idx = sort_rows(sort_list, index=True)
    return event_list[idx]
