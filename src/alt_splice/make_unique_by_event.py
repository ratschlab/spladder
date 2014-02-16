def make_unique_by_event(event_list):
# function event_list = make_unique_by_event(event_list)
#
# This script removes all events that share the sam alternative evnt coordinates
# but differ in the flanking size. The longest of several equal events is kept.

    rm_idx = []
    last_kept = 0
    for i in range(1, event_list.shape[0]):
        if i % 1000 == 0:
            print '.',
            if i % 10000 == 0:
                print '%i' % i
        
        old_coords = event_list[last_kept].get_inner_coords(trafo=True)
        curr_coords = event_list[i].get_inner_coords(trafo=True) 

        if old_coords.shape[0] == curr_coords.shape[0] and sp.all(old_coords == curr_coords):

            ### assertion that we did everything right
            assert(event_list[last_kept].chr_num == event_list[i].chr_num)
            assert(event_list[last_kept].strand == event_list[i].strand)
            
            ### check, which event is longer -> keep longer event
            len1 = event_list[last_kept].get_len()
            len2 = event_list[i].get_len()

            if len1 > len2:
                keep_idx = last_kept
                not_keep_idx = i
            else:
                keep_idx = i
                not_keep_idx = last_kept

            ### check if we would loose strains 
            idx = sp.where(~sp.in1d(event_list[not_keep_idx].strain, event_list[keep_idx].strain))[0]
            if idx.shape[0] > 0:
                event_list[keep_idx].strain = sp.r_[event_list[keep_idx].strain, event_list[not_keep_idx].strain[idx]]
                ### TODO !!!!!!!!!!!!! make sure that we keep different coordinates if the strains differ ...
                event_list[keep_idx].gene_name = sp.union1d(event_list[keep_idx].gene_name, event_list[not_keep_idx].gene_name)

            rm_idx.append(not_keep_idx)
            last_kept = keep_idx
        else:
            last_kept = i

    print 'events dropped: %i' % len(rm_idx)
    keep_idx = sp.where(~sp.in1d(sp.arange(event_list.shape[0])), rm_idx)[0]
    event_list = event_list[keep_idx]
