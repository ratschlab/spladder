def make_unique_by_strain(event_list):
# event_list = make_unique_by_strain(event_list)

    rm_idx = []
    for i in range(1, event_list.shape[0]):
        if i % 1000 == 0:
            print '.',
            if i % 10000 == 0:
                print '%i' % i

        old_coords = event_list[i-1].get_coords(trafo=True)
        curr_coords = event_list[i].get_coords(trafo=True) 

        if old_coords.shape[0] == curr_coords.shape[0] and sp.all(old_coords == curr_coords):

            ### assertion that we did everything right
            assert(event_list[i - 1].chr_num == event_list[i].chr_num)
            assert(event_list[i - 1].strand == event_list[i].strand)
            assert(event_list[i].strain.shape[0] == 1)

            idx = sp.where(event_list[i-1].strain == event_list[i].strain[0])[0]
            if idx.shape[0] > 0:
                assert(idx.shape[0] == 1)
                assert(sp.all(event_list[i].get_coords(trafo=True) == event_list[i-1].get_coords(trafo=True)))
                if not event_list[i].gene_name[0] in event_list[i-1].gene_name:
                    event_list[i-1].gene_name = sp.c_[event_list[i-1].gene_name, event_list[i].gene_name[0]]
                event_list[i] = event_list[i-1]
            else: 
                event_list[i].strain = sp.c_[event_list[i-1].strain[0], event_list[i].strain] ;
                assert(sp.all(sp.sort(event_list[i].strain) == sp.sort(sp.unique1d(event_list[i].strain))))
                ### TODO !!!!!!!!!!!!! make sure that we keep different coordinates if the strains differ ...
                if not event_list[i].gene_name[0] in event_list[i-1].gene_name:
                    event_list[i].gene_name = sp.c_[event_list[i-1].gene_name, event_list[i].gene_name[0]]
            rm_idx.append(i - 1)

    print 'events dropped: %i' % len(rm_idx)
    keep_idx = sp.where(~sp.in1d(sp.arange(event_list.shape[0])), rm_idx)[0]
    event_list = event_list[keep_idx]
