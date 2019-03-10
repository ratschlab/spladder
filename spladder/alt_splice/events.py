import scipy as sp

if __package__ is None:
    __package__ = 'modules.alt_splice'

from ..utils import *

def sort_events_full(event_list, options):

    if event_list.shape[0] == 0:
        return event_list

    coord_list = sp.array([x.get_coords() for x in event_list]) 
    chr_list = sp.array([options.chrm_lookup[x.chr] for x in event_list])
    strand_list = sp.array([x.strand == '-' for x in event_list], dtype = 'double')
    sort_list = sp.c_[chr_list, strand_list, coord_list]
    tmp, idx = sort_rows(sort_list, index=True)

    return event_list[idx]



def sort_events_by_event(event_list, options):
   
    coord_list = sp.array([x.get_inner_coords() for x in event_list], dtype='double') 
    chr_list = sp.array([options.chrm_lookup[x.chr] for x in event_list], dtype='double')
    strand_list = sp.array([x.strand == '-' for x in event_list], dtype = 'double')
    sort_list = sp.c_[chr_list, strand_list, coord_list]
    tmp, idx = sort_rows(sort_list, index=True)

    return event_list[idx]



def post_process_event_struct(events, options):

    if events.shape[0] == 0:
        return events

    ### filter out invalid coordinate projections
    is_valid = sp.array([sp.all(_.get_coords(trafo=True) > 0) for _ in events], dtype='bool')
    events = events[is_valid]

    ### sort exons in events
    for e in events:
        e.exons1 = sort_rows(e.exons1) 
        e.exons2 = sort_rows(e.exons2) 

    ### remove all events that have 0-length introns
    print('\nRemove 0-length intron events')
    is_valid = sp.array([_.get_intron_lens().min() > 0 for _ in events], dtype='bool')
    events = events[is_valid]
   
    ### sort events by all coordinates
    events = sort_events_full(events, options) 
    
    ### make events unique by strain
    print('\nMake %s events unique by strain' % events[0].event_type)
    events = make_unique_by_strain(events)

    ### sort events by event coordinates
    events = sort_events_by_event(events, options) 
    
    ### make events unique by event
    print('\nMake %s events unique by event' % events[0].event_type)
    events = make_unique_by_event(events)

    ### count detected strains
    for i in range(events.shape[0]):
        events[i].num_detected = len(events[i].strain)
        events[i].id = (i + 1)

    return events


def make_unique_by_strain(event_list):
    # event_list = make_unique_by_strain(event_list)

    rm_idx = []
    for i in range(1, event_list.shape[0]):
        if i % 1000 == 0:
            print('.', end=' ')
            if i % 10000 == 0:
                print('%i' % i)

        old_coords = event_list[i-1].get_coords(trafo=True)
        curr_coords = event_list[i].get_coords(trafo=True) 

        if old_coords.shape[0] == curr_coords.shape[0] and sp.all(old_coords == curr_coords):

            ### assertion that we did everything right
            if event_list[i - 1].chr == event_list[i].chr:
                assert(event_list[i - 1].strand == event_list[i].strand)
                assert(event_list[i].strain.shape[0] == 1)
            else:
                assert(event_list[i - 1].gene_name != event_list[i].gene_name)

            idx = sp.where(event_list[i-1].strain == event_list[i].strain[0])[0]
            if idx.shape[0] > 0:
                assert(idx.shape[0] == 1)
                assert(sp.all(event_list[i].get_coords(trafo=True) == event_list[i-1].get_coords(trafo=True)))
                if not event_list[i].gene_name[0] in event_list[i-1].gene_name:
                    event_list[i-1].gene_name = sp.r_[event_list[i-1].gene_name, [event_list[i].gene_name[0]]]
                event_list[i] = event_list[i-1]
            else: 
                event_list[i].strain = sp.r_[[event_list[i-1].strain[0]], event_list[i].strain]
                assert(sp.all(sp.sort(event_list[i].strain) == sp.sort(sp.unique(event_list[i].strain))))
                ### TODO !!!!!!!!!!!!! make sure that we keep different coordinates if the strains differ ...
                if not event_list[i].gene_name[0] in event_list[i-1].gene_name:
                    event_list[i].gene_name = sp.r_[event_list[i-1].gene_name, [event_list[i].gene_name[0]]]
            rm_idx.append(i - 1)

    print('events dropped: %i' % len(rm_idx))
    keep_idx = sp.where(~sp.in1d(sp.arange(event_list.shape[0]), rm_idx))[0]
    event_list = event_list[keep_idx]

    return event_list


def make_unique_by_event(event_list):
    # function event_list = make_unique_by_event(event_list)
    #
    # This script removes all events that share the sam alternative evnt coordinates
    # but differ in the flanking size. The shortes of several equal events is kept.

    rm_idx = []
    last_kept = 0
    for i in range(1, event_list.shape[0]):
        if i % 1000 == 0:
            print('.', end=' ')
            if i % 10000 == 0:
                print('%i' % i)
        
        old_coords = event_list[last_kept].get_inner_coords(trafo=True)
        curr_coords = event_list[i].get_inner_coords(trafo=True) 

        if old_coords.shape[0] == curr_coords.shape[0] and sp.all(old_coords == curr_coords):

            ### assertion that we did everything right
            assert(event_list[last_kept].chr == event_list[i].chr)
            assert(event_list[last_kept].strand == event_list[i].strand)
            
            ### check, which event is longer -> keep shorter event
            len1 = event_list[last_kept].get_len()
            len2 = event_list[i].get_len()

            if len1 < len2:
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

    print('events dropped: %i' % len(rm_idx))
    keep_idx = sp.where(~sp.in1d(sp.arange(event_list.shape[0]), rm_idx))[0]
    event_list = event_list[keep_idx]

    return event_list


def curate_alt_prime(event_list, options):
    # event_list = curate_alt_prime(event_list)

    if event_list.shape[0] == 0:
        return event_list

    rm_idx = []
    corr_count = 0

    for i in range(event_list.shape[0]):
        
        ### check if we have introns of zero length
        #if sp.any(event_list[i].exons1[:, 1] - event_list[i].exons1[:, 1] < 2) or sp.any(event_list[i].exons2[:, 1] - event_list[i].exons2[:, 1] < 2):
        if (event_list[i].exons1[1, 0] - event_list[i].exons1[0, 1] < 1) or (event_list[i].exons2[1, 0] - event_list[i].exons2[0, 1] < 1):
            rm_idx.append(i)
            continue

        ### check if alt exons overlap, otherwise we cannot curate (trim to shortest length)
        if (sp.all(event_list[i].exons1[0, :] == event_list[i].exons2[0, :]) and (event_list[i].exons1[1, 1] <= event_list[i].exons2[1, 0] or event_list[i].exons1[1, 0] >= event_list[i].exons2[1, 1])) or \
           (sp.all(event_list[i].exons1[1, :] == event_list[i].exons2[1, :]) and (event_list[i].exons1[0, 1] <= event_list[i].exons2[0, 0] or event_list[i].exons1[0, 0] >= event_list[i].exons2[0, 1])):
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
        keep_idx = sp.where(~sp.in1d(sp.arange(event_list.shape[0]), rm_idx))[0]
        event_list = event_list[keep_idx]

    print('Corrected %i events' % corr_count)
    print('Removed %i events' % len(rm_idx))

    return event_list
