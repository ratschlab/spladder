import numpy as np

if __package__ is None:
    __package__ = 'modules.alt_splice'

from ..utils import *

def sort_events_full(event_list, options):

    if event_list.shape[0] == 0:
        return event_list

    coord_list = np.array([x.get_coords() for x in event_list]) 
    chr_list = np.array([options.chrm_lookup[x.chr] for x in event_list])
    strand_list = np.array([x.strand == '-' for x in event_list], dtype = 'double')
    sort_list = np.c_[chr_list, strand_list, coord_list]
    tmp, idx = sort_rows(sort_list, index=True)

    return event_list[idx]



def sort_events_by_event(event_list, options):
   
    coord_list = np.array([x.get_inner_coords() for x in event_list], dtype='double') 
    chr_list = np.array([options.chrm_lookup[x.chr] for x in event_list], dtype='double')
    strand_list = np.array([x.strand == '-' for x in event_list], dtype = 'double')
    sort_list = np.c_[chr_list, strand_list, coord_list]
    tmp, idx = sort_rows(sort_list, index=True)

    return event_list[idx]



def post_process_event_struct(events, options):

    if events.shape[0] == 0:
        return events

    ### filter out invalid coordinate projections
    is_valid = np.array([np.all(_.get_coords() > 0) for _ in events], dtype='bool')
    events = events[is_valid]

    ### sort exons in events
    for e in events:
        e.exons1 = sort_rows(e.exons1) 
        e.exons2 = sort_rows(e.exons2) 

    ### remove all events that have 0-length introns
    print('\nRemove 0-length intron events')
    is_valid = np.array([_.get_intron_lens().min() > 0 for _ in events], dtype='bool')
    events = events[is_valid]
   
    ### sort events by all coordinates
    events = sort_events_full(events, options) 
    
    ### sort events by event coordinates
    events = sort_events_by_event(events, options) 
    
    ### make events unique by event
    print('\nMake %s events unique by event' % events[0].event_type)
    events = make_unique_by_event(events)

    ### create event_ids
    for i in range(events.shape[0]):
        events[i].id = (i + 1)

    return events


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
        
        old_coords = event_list[last_kept].get_inner_coords()
        curr_coords = event_list[i].get_inner_coords()

        if old_coords.shape[0] == curr_coords.shape[0] and \
           np.all(old_coords == curr_coords) and \
           event_list[last_kept].chr == event_list[i].chr and \
           event_list[last_kept].strand == event_list[i].strand:

            ### check, which event is longer -> keep shorter event
            len1 = event_list[last_kept].get_len()
            len2 = event_list[i].get_len()

            if len1 < len2:
                keep_idx = last_kept
                not_keep_idx = i
            else:
                keep_idx = i
                not_keep_idx = last_kept

            rm_idx.append(not_keep_idx)
            last_kept = keep_idx
        else:
            last_kept = i

    print('events dropped: %i' % len(rm_idx))
    keep_idx = np.where(~np.in1d(np.arange(event_list.shape[0]), rm_idx))[0]
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
        #if np.any(event_list[i].exons1[:, 1] - event_list[i].exons1[:, 1] < 2) or np.any(event_list[i].exons2[:, 1] - event_list[i].exons2[:, 1] < 2):
        if (event_list[i].exons1[1, 0] - event_list[i].exons1[0, 1] < 1) or (event_list[i].exons2[1, 0] - event_list[i].exons2[0, 1] < 1):
            rm_idx.append(i)
            continue

        ### check if alt exons overlap, otherwise we cannot curate (trim to shortest length)
        if (np.all(event_list[i].exons1[0, :] == event_list[i].exons2[0, :]) and (event_list[i].exons1[1, 1] <= event_list[i].exons2[1, 0] or event_list[i].exons1[1, 0] >= event_list[i].exons2[1, 1])) or \
           (np.all(event_list[i].exons1[1, :] == event_list[i].exons2[1, :]) and (event_list[i].exons1[0, 1] <= event_list[i].exons2[0, 0] or event_list[i].exons1[0, 0] >= event_list[i].exons2[0, 1])):
            continue
         
        if np.all(event_list[i].exons1[0, :] == event_list[i].exons2[0, :]):
            if event_list[i].exons1[1, 1] > event_list[i].exons2[1, 1]:
                event_list[i].exons1[1, 1] = event_list[i].exons2[1, 1]
                corr_count += 1
            elif event_list[i].exons1[1, 1] < event_list[i].exons2[1, 1]:
                event_list[i].exons2[1, 1] = event_list[i].exons1[1, 1]
                corr_count += 1
            ### check whether our isoform convention is still met - if not, correct
            if np.sum(event_list[i].exons1[:, 1] - event_list[i].exons1[:, 0]) > np.sum(event_list[i].exons2[:, 1] - event_list[i].exons2[:, 0]):
                event_list[i].exons1, event_list[i].exons2 = event_list[i].exons2, event_list[i].exons1
        elif np.all(event_list[i].exons1[1, :] == event_list[i].exons2[1, :]):
            if event_list[i].exons1[0, 0] > event_list[i].exons2[0, 0]:
                event_list[i].exons2[0, 0] = event_list[i].exons1[0, 0]
                corr_count += 1
            elif event_list[i].exons1[0, 0] < event_list[i].exons2[0, 0]:
                event_list[i].exons1[0, 0] = event_list[i].exons2[0, 0]
                corr_count += 1
            ### check whether our isoform convention is still met - if not, correct
            if np.sum(event_list[i].exons1[:, 1] - event_list[i].exons1[:, 0]) > np.sum(event_list[i].exons2[:, 1] - event_list[i].exons2[:, 0]):
                event_list[i].exons1, event_list[i].exons2 = event_list[i].exons2, event_list[i].exons1

    ### remove events with non-overlapping alt_exons
    if len(rm_idx) > 0:
        keep_idx = np.where(~np.in1d(np.arange(event_list.shape[0]), rm_idx))[0]
        event_list = event_list[keep_idx]

    print('Corrected %i events' % corr_count)
    print('Removed %i events' % len(rm_idx))

    return event_list
