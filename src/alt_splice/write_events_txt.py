import sys

def write_events_txt(fn_out_txt, strains, events, anno_fn=None):
# write_events_txt(fn_out_txt, strains, events, anno_fn)
    
    if events.shape[0] == 0:
        print >> sys.stderr, 'WARNING: No events present.'
        return

    if anno_fn is not None:
       anno = load(anno_fn);
       anno_names, s_idx = sort([x.name for x in genes])
       anno = anno[s_idx]

    print 'writing %s events in flat txt format to %s' % (events[0].event_type, fn_out_txt)

    fd = open(fn_out_txt, 'w+')

    if anno_fn is not None:
        gene_header = '\tgene_start\tgene_end'
    else:
        gene_header = ''

    if events[0].event_type == 'exon_skip':
        print >> fd, 'contig\tstrand\tevent_id\tgene_name%s\texon_pre_start\texon_pre_end\texon_start\texon_end\texon_aft_start\texon_aft_end' % gene_header,
        for i in range(len(strains)):
            print >> fd, '\t%s:exon_pre_cov\t%s:exon_cov\t%s:exon_aft_cov\t%s:intron_pre_conf\t%s:intron_aft_conf\t%s:intron_skip_conf' % (strains[i], strains[i], strains[i], strains[i], strains[i], strains[i]),
            print >> fd, '\n',
    elif events[0].event_type in ['alt_3prime', 'alt_5prime']:
        print >> fd, 'contig\tstrand\tevent_id\tgene_name%s\texon_const_start\texon_const_end\texon_alt1_start\texon_alt1_end\texon_alt2_start\texon_alt2_end' % gene_header,
        for i in range(len(strains)):
            print >> fd, '\t%s:exon_diff_cov\t%s:exon_const_cov\t%s:intron1_conf\t%s:intron2_conf' % (strains[i], strains[i], strains[i], strains[i])
            print >> fd, '\n',
    elif events[0].event_type == 'intron_retention':
        print >> fd, 'contig\tstrand\tevent_id\tgene_name%s\texon1_start\texon1_end\tintron_start\tintron_end\texon2_start\texon2_end' % gene_header,
        for i in range(len(strains)):
            print >> fd, '\t%s:exon1_cov\t%s:intron_cov\t%s:exon2_cov\t%s:intron_conf' % (strains[i], strains[i], strains[i], strains[i]),
            print >> fd, '\n',
    elif strcmp(events(1).event_type, 'mult_exon_skip'),
        print >> fd, 'contig\tstrand\tevent_id\tgene_name%s\texon_pre_start\texon_pre_end\texon_starts\texon_ends\texon_aft_start\texon_aft_end' % gene_header,
            print >> fd, '\t%s:exon_pre_cov\t%s:exon_cov\t%s:exon_aft_cov\t%s:intron_pre_conf\t%s:intron_inner_conf\t%s:exon_inner_count\t%s:intron_aft_conf\t%s:intron_skip_conf' % (strains[i], strains[i], strains[i], strains[i], strains[i], strains[i], strains[i], strains[i]),
            print >> fd, '\n',
    else:
        print >> sys.stderr, 'Unknown event type: %s' % events[0].event_type

    for i in range(events.shape[0]):
        print >> fd, '%s\t%c\t%s_%i\t%s' % (events[i].chr, events[i].strand, events[i].event_type, events[i].id, events[i].gene_name[0]),
        if anno_fn is not None:
            a_idx = anno_names.index(events[i].gene_name[0])
            print >> fd, '\t%i\t%i' % (anno[a_idx].start, anno[a_idx].stop)

        ev = events[i]
        if ev.event_type == 'exon_skip':
            print >> fd, '\t%i\t%i\t%i\t%i\t%i\t%i' % (ev.exons2[0, 0], ev.exons2[0, 1], ev.exons2[1, 0], ev.exons2[1, 1], ev.exons2[2, 0], ev.exons2[2, 1]),
            for j in range(len(strains)):
                if ev.info[j]['valid']:
                    print >> fd, '\t%.1f\t%.1f\t%.1f\t%i\t%i\t%i' % (ev.info[j]['exon_pre_cov'], ev.info[j]['exon_cov'], ev.info[j]['exon_aft_cov'], ev.info[j]['exon_pre_exon_conf'], ev.info[j]['exon_exon_aft_conf'], ev.info[j]['exon_pre_exon_aft_conf']),
                else:
                    print >> fd, '\t-1\t-1\t-1\t-1',
        elif ev.event_type == 'intron_retention':
            print >> fd, '\t%i\t%i\t%i\t%i\t%i\t%i' % (ev.exons1[0, 0], ev.exons1[0, 1], ev.exons1[0, 1], ev.exons1[1, 0], ev.exons1[1, 0], ev.exons1[1, 1,]),
            for j in range(len(strains)):
                if ev.info[j]['valid']:
                    print >> fd, '\t%.1f\t%.1f\t%.1f\t%i' % (ev.info[j]['exon1_cov'], ev.info[j]['intron_cov'], ev.info[j]['exon2_cov'], ev.info[j]['intron_conf']),
                else:
                    print >> fd, '\t-1\t-1\t-1',
        elif ev.event_type in ['alt_3prime', 'alt_5prime']:
            if sp.all(ev.exons1[0, :] == ev.exons2[0, :]):
                print >> fd, '\t%i\t%i\t%i\t%i\t%i\t%i' % (ev.exons1[0, 0], ev.exons1[0, 1], ev.exons1[1, 0], ev.exons1[1, 1], ev.exons2[1, 0], ev.exons2[1, 1]),
            elif sp.all(ev.exons1[1, :] == ev.exons2[1, :]):
                print >> fd, '\t%i\t%i\t%i\t%i\t%i\t%i' % (ev.exons1[1, 0], ev.exons1[1, 1], ev.exons1[0, 0], ev.exons1[0, 1], ev.exons2[0, 0], ev.exons2[0, 1]),
            for j in range(len(strains)):
                if ev.info[j]['valid']:
                    print >> fd, '\t%1.1f\t%1.1f\t%i\t%i' % (ev.info[j]['exon_diff_cov'], ev.info[j]['exon_const_cov'], ev.info[j]['intron1_conf'], ev.info[j]['intron2_conf']),
                else:
                    print >> fd, '\t-1\t-1\t-1\t-1',
        elif ev.event_type == 'mult_exon_skip':
            print >> fd, '\t%i\t%i' % (ev.exons2[0, 0], ev.exons2[0, 1]),
            starts = '%i', ev.exons2[1, 0]
            ends = '%i', ev.exons2[1, 1]
            for k in range(2, ev.exons2.shape[0] - 1):
                starts = '%s:%i' % (starts, ev.exons2[k, 0])
                ends = '%s:%i' % (ends, ev.exons2[k, 1])
            print >> fd, '\t%s\t%s\t%i\t%i' % (starts, ends, ev.exons2[-1, 0], ev.exons2[-1, 1]),
            for j in range(len(strains)):
                if ev.info[j]['valid']:
                    print >> fd, '\t%.1f\t%.1f\t%.1f\t%i\t%i\t%i\t%i\t%i' % (ev.info[j]['exon_pre_cov'], ev.info[j]['exons_cov'], ev.info[j]['exon_aft_cov'], ev.info[j]['exon_pre_exon_conf'], ev.info[j]['sum_inner_exon_conf'], ev.info[j]['num_inner_exon'], ev.info[j]['exon_exon_aft_conf'], ev.info[j]['exon_pre_exon_aft_conf']),
                else:
                    print >> fd, '\t-1\t-1\t-1\t-1\t-1\t-1\t-1',
        print >> fd, '\n',
    fd.close()
