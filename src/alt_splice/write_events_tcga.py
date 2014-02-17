import sys

def write_events_tcga(fn_out, strains, events):
# write_events_tcga(fn_out, strains, events)

    if events.shape[0] == 0:
        print >> sys.stderr, 'WARNING: No events present.'
        return

    event_type = events[0].event_type

    print 'writing %s events in tcga format to %s' % (events[0].event_type, fn_out)

    fd = open(fn_out, 'w+')
    print >> fd, 'gene\teventtype\tcoordinates',
    for s in strains:
        print >> fd, '\t%s', s,
    print >> fd, '\n',

    for i in range(events.shape[0]):
        print >> fd, '%s\t%s\t%s:' % (events[i].gene_name[0], event_type, events[i].chr),
        if event_type == 'intron_retention':
            print >> fd, ':%i-%i:%i-%i' % (events[i].exons1[0, 0], events[i].exons1[0, 1], events[i].exons1[1, 0], exons1[1, 1]),
        elif event_type in ['alt_3prime', 'alt_5prime']:
            if sp.all(events[i].exons1[0, :] == events[i].exons2[0, :]):
                print >> fd, ':%i-%i:%i-%i:%i-%i' % (events[i].exons1[0, 0], events[i].exons1[0, 1], events[i].exons1[1, 0], events[i].exons1[1, 1], events[i].exons2[1, 0], events[i].exons2[1, 1]),
            else:
                print >> fd, ':%i-%i:%i-%i:%i-%i' % (events[i].exons1[1, 0], events[i].exons1[1, 1], events[i].exons1[0, 0], events[i].exons1[0, 1], events[i].exons2[0, 0], events[i].exons2[0, 1]),
        elif event_type == 'exon_skip':
            print >> fd, ':%i-%i:%i-%i:%i-%i' % (events[i].exons1[0, 0], events[i].exons1[0, 1], events[i].exons1[1, 0], events[i].exons1[1, 1], events[i].exons1[2, 0], events[i].exons1[2, 1]),
        elif event_type == 'mult_exon_skip':
            print >> fd, ':%i-%i' % (events[i].exons2[0, 0], events[i].exons2[0, 1]),
            for j in range(1, events.shape[0] - 1):
                print >> fd, ':%i-%i:' % (events[i].exons2[j, 0], events[i].exons2[j, 1]),
            print >> fd, ':%i-%i' % (events[i].exons2[-1, 0], events[i].exons2[-1, 1]),
        for j in range(len(strains)):
            if events[i].info[j]['valid']:
                if events[i].event_type in ['alt_3prime', 'alt_5prime']:
                    if (events[i].exons1[1, 0] - events[i].exons1[0, 1]) < (events[i].exons2[1, 0] - events[i].exons2[0, 1]):
                        num = events[i].info[j]['intron1_conf']
                    else:
                        num = events[i].info[j]['intron2_conf']
                    denom = events[i].info[j]['intron1_conf'] + events[i].info[j]['intron2_conf']
                    confirmation = denom
                elif events[i].event_type == 'exon_skip':
                    num = events[i].info[j]['exon_pre_exon_conf'] + events[i].info[j]['exon_exon_aft_conf']
                    denom = events[i].info[j]['exon_pre_exon_conf'] + events[i].info[j]['exon_exon_aft_conf'] + (2 * events[i].info[j]['exon_pre_exon_aft_conf'])
                    confirmation = events[i].info[j]['exon_pre_exon_conf'] + events[i].info[j]['exon_exon_aft_conf'] + events[i].info[j]['exon_pre_exon_aft_conf']
                elif events[i].event_type == 'mult_exon_skip':
                    num = events[i].info[j]['exon_pre_exon_conf'] + events[i].info[j]['sum_inner_exon_conf'] + events[i].info[j]['exon_exon_aft_conf']
                    denom = events[i].info[j]['exon_pre_exon_conf'] + events[i].info[j]['sum_inner_exon_conf'] + events[i].info[j]['exon_exon_aft_conf'] + ((2 + events[i].info[j]['num_inner_exon']) * events[i].info[j]['exon_pre_exon_aft_conf'])
                    confirmation = events[i].info[j]['exon_pre_exon_conf'] + events[i].info[j]['sum_inner_exon_conf'] + events[i].info[j]['exon_exon_aft_conf'] + events[i].info[j]['exon_pre_exon_aft_conf']
                elif events[i].event_type == 'intron_retention':
                    num = events[i].info[j]['intron_conf']
                    denom = 1
                    confirmation = num
                else:
                    print >> sys.stderr, 'Unknown event type: %s' % (events[i].event_type)
                    sys.exit(1)

                if confirmation < 10:
                    print >> fd, '\tNA',
                else:
                    print >> fd, '\t%1.1f', % (num / denom),
            else:
                print >> fd, '\tNA',
        print >> fd, '\n',
    fd.close()
