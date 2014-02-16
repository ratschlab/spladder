def curate_alt_prime(event_list):
# event_list = curate_alt_prime(event_list)

    rm_idx = []
    corr_count = 0

    for i in range(event_list.shape[0]):
        
        ### check if we have introns of zero length
        if sp.any(event_list[i].exons1[:, 1] - event_list[i].exons1[:, 1] < 2) or sp.any(event_list[i].exons2[:, 1] - event_list[i].exons2[:, 1] < 2):
            rm_idx.append(i)
            continue

        ### check if alt exons overlap, otherwise we cannot curate (trim to shortest length)
        if (sp.all(event_list[i].exons1[0, :] == event_list[i].exons2[0, :]) and (event_list[i].exons1[1, 1] < event_list[i].exons2[1, 0] or event_list[i].exons1[1, 0] > event_list[i].exons2[1, 1])) or
           (sp.all(event_list[i].exons1[1, :] == event_list[i].exons2[1, :]) and (event_list[i].exons1[0, 1] < event_list[i].exons2[0, 0] or event_list[i].exons1[0, 0] > event_list[i].exons2[0, 1])):
            continue
         
        if sp.all(event_list[i].exons1[0, :] == event_list[i].exons2[0, :]):
            if event_list[i].exons1[1, 1] > event_list[i].exons2[1, 1]:
                event_list[i].exons1[1, 1] = event_list[i].exons2[1, 1]
                corr_count += 1
            elif event_list[i].exons1[1, 1] < event_list[i].exons2[1, 1]:
                event_list[i].exons2[1, 1] = event_list[i].exons1[1, 1]
                corr_count += 1
        elif sp.all(event_list[i].exons1[1, :] == event_list[i].exons2[1, :]):
            if event_list[i].exons1[0, 0] > event_list[i].exons2[0, 0]:
                event_list[i].exons2[0, 0] = event_list[i].exons1[0, 0]
                corr_count += 1
            elif event_list[i].exons1[0, 0] < event_list[i].exons2[0, 0]:
                event_list[i].exons1[0, 0] = event_list[i].exons2[0, 0]
                corr_count += 1

    ### remove events with non-overlapping alt_exons
    if len(rm_idx) > 0:
        keep_idx = sp.where(~sp.in1d(sp.arange(event_list.shape[0])), rm_idx)[0]
        event_list = event_list[keep_idx]

    print 'Corrected %i events' % corr_count
    print 'Removed %i events' % len(rm_idx)
