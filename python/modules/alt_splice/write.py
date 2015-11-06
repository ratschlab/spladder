import sys
import os
import scipy as sp
import h5py

def write_events_txt(fn_out_txt, strains, events, fn_counts, event_idx=None, anno_fn=None):
    # write_events_txt(fn_out_txt, strains, events, fn_counts, event_idx, anno_fn)
    
    if events.shape[0] == 0:
        print >> sys.stderr, 'WARNING: No events present.'
        return

    if event_idx is None:
        event_idx = sp.arange(events.shape[0])

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
        fd.write('contig\tstrand\tevent_id\tgene_name%s\texon_pre_start\texon_pre_end\texon_start\texon_end\texon_aft_start\texon_aft_end' % gene_header)
        for i in range(len(strains)):
            fd.write('\t%s:exon_pre_cov\t%s:exon_cov\t%s:exon_aft_cov\t%s:intron_pre_conf\t%s:intron_aft_conf\t%s:intron_skip_conf' % (strains[i], strains[i], strains[i], strains[i], strains[i], strains[i]))
    elif events[0].event_type in ['alt_3prime', 'alt_5prime']:
        fd.write('contig\tstrand\tevent_id\tgene_name%s\texon_const_start\texon_const_end\texon_alt1_start\texon_alt1_end\texon_alt2_start\texon_alt2_end' % gene_header)
        for i in range(len(strains)):
            fd.write('\t%s:exon_diff_cov\t%s:exon_const_cov\t%s:intron1_conf\t%s:intron2_conf' % (strains[i], strains[i], strains[i], strains[i]))
    elif events[0].event_type == 'intron_retention':
        fd.write('contig\tstrand\tevent_id\tgene_name%s\texon1_start\texon1_end\tintron_start\tintron_end\texon2_start\texon2_end' % gene_header)
        for i in range(len(strains)):
            fd.write('\t%s:exon1_cov\t%s:intron_cov\t%s:exon2_cov\t%s:intron_conf' % (strains[i], strains[i], strains[i], strains[i]))
    elif events[0].event_type == 'mult_exon_skip':
        fd.write('contig\tstrand\tevent_id\tgene_name%s\texon_pre_start\texon_pre_end\texon_starts\texon_ends\texon_aft_start\texon_aft_end' % gene_header)
        for i in range(len(strains)):
            fd.write('\t%s:exon_pre_cov\t%s:exon_cov\t%s:exon_aft_cov\t%s:intron_pre_conf\t%s:intron_inner_conf\t%s:exon_inner_count\t%s:intron_aft_conf\t%s:intron_skip_conf' % (strains[i], strains[i], strains[i], strains[i], strains[i], strains[i], strains[i], strains[i]))
    elif events[0].event_type == 'mutex_exons':
        fd.write('contig\tstrand\tevent_id\tgene_name%s\texon_pre_start\texon_pre_end\texon1_start\texon1_end\texon2_start\texon2_end\texon_aft_start\texon_aft_end' % gene_header)
        for i in range(len(strains)):
            fd.write('\t%s:exon_pre_cov\t%s:exon1_cov\t%s:exon2_cov\t%s:exon_aft_cov\t%s:pre_exon1_conf\t%s:pre_exon2_conf\t%s:exon1_aft_conf\t%s:exon2_aft_conf' % (strains[i], strains[i], strains[i], strains[i], strains[i], strains[i], strains[i], strains[i]))
    else:
        print >> sys.stderr, 'Unknown event type: %s' % events[0].event_type
    fd.write('\n')
    fd.flush()

    ### load data from count hdf5
    IN = h5py.File(fn_counts, 'r')

    for i in event_idx:
        fd.write('%s\t%c\t%s_%i\t%s' % (events[i].chr, events[i].strand, events[i].event_type, events[i].id, events[i].gene_name[0]))
        if anno_fn is not None:
            a_idx = anno_names.index(events[i].gene_name[0])
            fd.write('\t%i\t%i' % (anno[a_idx].start, anno[a_idx].stop))

        ev = events[i]
        counts = IN['event_counts'][:, :, i]
        if ev.event_type == 'exon_skip':
            fd.write('\t%i\t%i\t%i\t%i\t%i\t%i' % (ev.exons2[0, 0] + 1, ev.exons2[0, 1], ev.exons2[1, 0] + 1, ev.exons2[1, 1], ev.exons2[2, 0] + 1, ev.exons2[2, 1]))
            for j in range(len(strains)):
                if counts[j, 0] == 1:
                    fd.write('\t%.1f\t%.1f\t%.1f\t%i\t%i\t%i' % (counts[j, 2], counts[j, 1], counts[j, 3], counts[j, 4], counts[j, 5], counts[j, 6]))
                else:
                    fd.write('\t-1\t-1\t-1\t-1')
        elif ev.event_type == 'intron_retention':
            fd.write('\t%i\t%i\t%i\t%i\t%i\t%i' % (ev.exons1[0, 0] + 1, ev.exons1[0, 1], ev.exons1[0, 1] + 1, ev.exons1[1, 0], ev.exons1[1, 0] + 1, ev.exons1[1, 1,]))
            for j in range(len(strains)):
                if counts[j, 0] == 1:
                    fd.write('\t%.1f\t%.1f\t%.1f\t%i' % (counts[j, 2], counts[j, 1], counts[j, 3], counts[j, 4]))
                else:
                    fd.write('\t-1\t-1\t-1')
        elif ev.event_type in ['alt_3prime', 'alt_5prime']:
            if sp.all(ev.exons1[0, :] == ev.exons2[0, :]):
                fd.write('\t%i\t%i\t%i\t%i\t%i\t%i' % (ev.exons1[0, 0] + 1, ev.exons1[0, 1], ev.exons1[1, 0] + 1, ev.exons1[1, 1], ev.exons2[1, 0] + 1, ev.exons2[1, 1]))
            elif sp.all(ev.exons1[1, :] == ev.exons2[1, :]):
                fd.write('\t%i\t%i\t%i\t%i\t%i\t%i' % (ev.exons1[1, 0] + 1, ev.exons1[1, 1], ev.exons1[0, 0] + 1, ev.exons1[0, 1], ev.exons2[0, 0] + 1, ev.exons2[0, 1]))
            for j in range(len(strains)):
                if counts[j, 0]:
                    fd.write('\t%1.1f\t%1.1f\t%i\t%i' % (counts[j, 1], counts[j, 2], counts[j, 3], counts[j, 4]))
                else:
                    fd.write('\t-1\t-1\t-1\t-1')
        elif ev.event_type == 'mult_exon_skip':
            fd.write('\t%i\t%i' % (ev.exons2[0, 0] + 1, ev.exons2[0, 1]))
            starts = '%i' % (ev.exons2[1, 0] + 1)
            ends = '%i' % ev.exons2[1, 1]
            for k in range(2, ev.exons2.shape[0] - 1):
                starts = '%s:%i' % (starts, ev.exons2[k, 0] + 1)
                ends = '%s:%i' % (ends, ev.exons2[k, 1])
            fd.write('\t%s\t%s\t%i\t%i' % (starts, ends, ev.exons2[-1, 0] + 1, ev.exons2[-1, 1]))
            for j in range(len(strains)):
                if counts[j, 0]:
                    fd.write('\t%.1f\t%.1f\t%.1f\t%i\t%i\t%i\t%i\t%i' % (counts[j, 1], counts[j, 2], counts[j, 3], counts[j, 4], counts[j, 7], counts[j, 8], counts[j, 5], counts[j, 6]))
                else:
                    fd.write('\t-1\t-1\t-1\t-1\t-1\t-1\t-1')
        elif ev.event_type == 'mutex_exons':
            fd.write('\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i' % (ev.exons1[0, 0] + 1, ev.exons1[0, 1], ev.exons1[1, 0] + 1, ev.exons1[1, 1], ev.exons2[1, 0] + 1, ev.exons2[1, 1], ev.exons2[2, 0] + 1, ev.exons2[2, 1]))
            for j in range(len(strains)):
                if counts[j, 0] == 1:
                    fd.write('\t%.1f\t%.1f\t%.1f\t%.1f\t%i\t%i\t%i\t%i' % (counts[j, 1], counts[j, 2], counts[j, 3], counts[j, 4], counts[j, 5], counts[j, 6], counts[j, 7], counts[j, 8]))
                else:
                    fd.write('\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1')
        fd.write('\n')
        fd.flush()
    fd.close()
    IN.close()


def write_events_tcga(fn_out, strains, events, fn_counts, event_idx=None):
    # write_events_tcga(fn_out, strains, events, fn_counts, event_idx=None)

    if events.shape[0] == 0:
        print >> sys.stderr, 'WARNING: No events present.'
        return

    if event_idx is None:
        event_idx = sp.arange(events.shape[0])

    event_type = events[0].event_type

    print 'writing %s events in tcga format to %s' % (events[0].event_type, fn_out)

    ### load counts from hdf5
    IN = h5py.File(fn_counts, 'r')

    fd = open(fn_out, 'w+')
    print >> fd, 'gene\teventtype\tcoordinates',
    for s in strains:
        print >> fd, '\t%s', s,
    print >> fd, '\n',

    for i in range(events.shape[0]):
        counts = IN['event_counts'][:, :, i]

        print >> fd, '%s\t%s\t%s:' % (events[i].gene_name[0], event_type, events[i].chr),
        if event_type == 'intron_retention':
            print >> fd, ':%i-%i:%i-%i' % (events[i].exons1[0, 0], events[i].exons1[0, 1], events[i].exons1[1, 0], events[i].exons1[1, 1]),
        elif event_type in ['alt_3prime', 'alt_5prime']:
            if sp.all(events[i].exons1[0, :] == events[i].exons2[0, :]):
                print >> fd, ':%i-%i:%i-%i:%i-%i' % (events[i].exons1[0, 0], events[i].exons1[0, 1], events[i].exons1[1, 0], events[i].exons1[1, 1], events[i].exons2[1, 0], events[i].exons2[1, 1]),
            else:
                print >> fd, ':%i-%i:%i-%i:%i-%i' % (events[i].exons1[1, 0], events[i].exons1[1, 1], events[i].exons1[0, 0], events[i].exons1[0, 1], events[i].exons2[0, 0], events[i].exons2[0, 1]),
        elif event_type == 'exon_skip':
            print >> fd, ':%i-%i:%i-%i:%i-%i' % (events[i].exons2[0, 0], events[i].exons2[0, 1], events[i].exons2[1, 0], events[i].exons2[1, 1], events[i].exons2[2, 0], events[i].exons2[2, 1]),
        elif event_type == 'mult_exon_skip':
            print >> fd, ':%i-%i' % (events[i].exons2[0, 0], events[i].exons2[0, 1]),
            for j in range(1, events[i].exons2.shape[0] - 1):
                print >> fd, ':%i-%i:' % (events[i].exons2[j, 0], events[i].exons2[j, 1]),
            print >> fd, ':%i-%i' % (events[i].exons2[-1, 0], events[i].exons2[-1, 1]),
        elif event_type == 'mutex_exons':
            print >> fd, ':%i-%i:%i-%i:%i-%i:%i-%i' % (events[i].exons1[0, 0], events[i].exons1[0, 1], events[i].exons1[1, 0], events[i].exons1[1, 1], events[i].exons2[1, 0], events[i].exons2[1, 1], events[i].exons1[2, 0], events[i].exons1[2, 1]),
        for j in range(len(strains)):
            if counts[j, 0] == 1:
                if events[i].event_type in ['alt_3prime', 'alt_5prime']:
                    if (events[i].exons1[1, 0] - events[i].exons1[0, 1]) < (events[i].exons2[1, 0] - events[i].exons2[0, 1]):
                        num = counts[j, 3]
                    else:
                        num = counts[j, 4]
                    denom = counts[j, 3] + counts[j, 4]
                    confirmation = denom
                elif events[i].event_type == 'exon_skip':
                    num = counts[j, 4] + counts[j, 5]
                    denom = counts[j, 4] + counts[j, 5] + (2 * counts[j, 6])
                    confirmation = counts[j, 4] + counts[j, 5] + counts[j, 6]
                elif events[i].event_type == 'mult_exon_skip':
                    num = counts[j, 4] + counts[j, 7] + counts[j, 5]
                    denom = counts[j, 4] + counts[j, 7] + counts[j, 5] + ((2 + counts[j, 8]) * counts[j, 6]) 
                    confirmation = counts[j, 4] + counts[j, 7] + counts[j, 5] + counts[j, 6]
                elif events[i].event_type == 'intron_retention':
                    num = counts[j, 4]
                    denom = 1
                    confirmation = num
                elif events[i].event_type == 'mutex_exons':
                    num = counts[j, 5] + counts[j, 6]
                    denom = counts[j, 5] + counts[j, 6] + counts[j, 7] + counts[j, 8]
                    confirmation = denom
                else:
                    print >> sys.stderr, 'Unknown event type: %s' % (events[i].event_type)
                    sys.exit(1)

                if confirmation < 10:
                    print >> fd, '\tNA',
                else:
                    print >> fd, '\t%1.1f' % (num / float(denom)),
            else:
                print >> fd, '\tNA',
        print >> fd, '\n',
    fd.close()


def write_events_gff3(fn_out_gff3, events, idx=None, as_gtf=False):
    
    if events.shape[0] == 0:
        print >> sys.stderr, 'WARNING: No events present.'
        return

    if idx is None:
        idx = sp.arange(events.shape[0])

    print 'writing %s events in gff3 format to %s' % (events[0].event_type, fn_out_gff3)

    fd_out = open(fn_out_gff3, 'w+') 
    print >> fd_out, '##gff-version 3'

    ### load gene structure
    for i in idx:

        ev = events[i]
        gene_name = events[i].gene_name[0] ### TODO - why only first?
        start_pos = ev.exons1[0, 0]
        stop_pos = ev.exons1[-1, -1]

        ### get order of isoforms o_idx(1) -> iso1 and o_idx(2) -> iso2
        ### assert that first isoform is always the shorter one
        #if ev.event_type == 'intron_retention':
        #    assert(sp.sum(ev.exons1[:, 1] - ev.exons1[:, 0]) < sp.sum(ev.exons2[1] - ev.exons2[0]))
        #else:
        #    assert(sp.sum(ev.exons1[:, 1] - ev.exons1[:, 0]) < sp.sum(ev.exons2[:, 1] - ev.exons2[:, 0]))

        name = '%s.%i' % (ev.event_type, ev.id)

        if as_gtf:
            print >> fd_out, '%s\t%s\tgene\t%i\t%i\t.\t%c\t.\tgene_id "%s"; transcript_id "%s"; gene_name "%s";' % (ev.chr, ev.event_type, start_pos, stop_pos, ev.strand, name, name, gene_name)
            print >> fd_out, '%s\t%s\ttranscript\t%i\t%i\t.\t%c\t.\tgene_id "%s"; transcript_id "%s_iso1"; gene_name "%s";' % (ev.chr, ev.event_type, start_pos, stop_pos, ev.strand, name, name, gene_name) 
            for i in range(ev.exons1.shape[0]):
                print >> fd_out, '%s\t%s\texon\t%i\t%i\t.\t%c\t.\tgene_id "%s"; transcript_id "%s_iso1"; exon_id "%s_iso1_exon%i"' % (ev.chr, ev.event_type, ev.exons1[i, 0], ev.exons1[i, 1], ev.strand, name, name, name, i+1)
            print >> fd_out, '%s\t%s\tmRNA\t%i\t%i\t.\t%c\t.\tgene_id "%s"; transcript_id "%s_iso2"; gene_name "%s"' % (ev.chr, ev.event_type, start_pos, stop_pos, ev.strand, name, name, gene_name)
            ex_cnt = 1
            if ev.event_type == 'intron_retention':
                print >> fd_out, '%s\t%s\texon\t%i\t%i\t.\t%c\t.\tgene_id "%s"; transcript_id "%s_iso2"; exon_id "%s_iso2_exon1";' % (ev.chr, ev.event_type, ev.exons2[0], ev.exons2[1], ev.strand, name, name, name)
                ex_cnt += 1
            else:
                for i in range(ev.exons2.shape[0]):
                    print >> fd_out, '%s\t%s\texon\t%i\t%i\t.\t%c\t.\tgene_id "%s"; transcript_id "%s_iso2"; exon_id "%s_iso2_exon%i";' % (ev.chr, ev.event_type, ev.exons2[i, 0], ev.exons2[i, 1], ev.strand, name, name, name, ex_cnt)
                    ex_cnt += 1
        else:
            print >> fd_out, '%s\t%s\tgene\t%i\t%i\t.\t%c\t.\tID=%s;GeneName="%s"' % (ev.chr, ev.event_type, start_pos, stop_pos, ev.strand, name, ev.gene_name[0])
            print >> fd_out, '%s\t%s\tmRNA\t%i\t%i\t.\t%c\t.\tID=%s_iso1;Parent=%s;GeneName="%s"' % (ev.chr, ev.event_type, start_pos, stop_pos, ev.strand, name, name, ev.gene_name[0]) 
            for i in range(ev.exons1.shape[0]):
                print >> fd_out, '%s\t%s\texon\t%i\t%i\t.\t%c\t.\tParent=%s_iso1' % (ev.chr, ev.event_type, ev.exons1[i, 0], ev.exons1[i, 1], ev.strand, name)
            print >> fd_out, '%s\t%s\tmRNA\t%i\t%i\t.\t%c\t.\tID=%s_iso2;Parent=%s;GeneName="%s"' % (ev.chr, ev.event_type, start_pos, stop_pos, ev.strand, name, name, ev.gene_name[0])
            if ev.event_type == 'intron_retention':
                print >> fd_out, '%s\t%s\texon\t%i\t%i\t.\t%c\t.\tParent=%s_iso2' % (ev.chr, ev.event_type, ev.exons2[0], ev.exons2[1], ev.strand, name)
            else:
                for i in range(ev.exons2.shape[0]):
                    print >> fd_out, '%s\t%s\texon\t%i\t%i\t.\t%c\t.\tParent=%s_iso2' % (ev.chr, ev.event_type, ev.exons2[i, 0], ev.exons2[i, 1], ev.strand, name)
    fd_out.close()


def write_events_structured(fn_out_struc, events, idx=None):
    
    if events.shape[0] == 0:
        print >> sys.stderr, 'WARNING: No events present.'
        return

    if idx is None:
        idx = sp.arange(events.shape[0])

    print 'writing %s events in generic structured format to %s' % (events[0].event_type, fn_out_struc)
    mult_exon_skip_bool = True

    fd_out = open(fn_out_struc, 'w+') 
    
    ### load gene structure
    for i in idx:

        ev = events[i]
        gene_name = events[i].gene_name[0] ### TODO - why only first?
        start_pos = ev.exons1[0, 0]
        stop_pos = ev.exons1[-1, -1]

        name = '%s.%i' % (ev.event_type, ev.id)

        print >> fd_out, '%s\tundefined\t%s\t%i\t%i\t.\t%c\t.\tgene_id "%s"; transcript_id "%s"; gene_name "%s";' % (ev.chr, ev.event_type, start_pos, stop_pos, ev.strand, name, name, gene_name),
        ### handle + strand
        if ev.strand == '+':
            if ev.event_type == 'exon_skip':
                struc = '0,1-2^'
                flanks = '%i^,%i-' % (ev.exons1[0, 1], ev.exons1[-1, 0] + 1)
                schain = ',%i-%i^' % (ev.exons2[1, 0] + 1, ev.exons2[1, 1])
            elif ev.event_type in ['alt_3prime', 'alt_5prime']:
                if sp.all(ev.exons1[0, :] == ev.exons2[0, :]):
                    struc = '1-,2-'
                    flanks = '%i^,%i^' % (ev.exons1[0, 1], min(ev.exons1[-1, 1], ev.exons2[-1, 1]))
                    schain = '%i-,%i-' % (min(ev.exons1[-1, 0] + 1, ev.exons2[-1, 0]) + 1, max(ev.exons1[-1, 0], ev.exons2[-1, 0]) + 1)
                elif sp.all(ev.exons1[-1, :] == ev.exons2[-1, :]):
                    struc = '1^,2^'
                    flanks = '%i-,%i-' % (max(ev.exons1[0, 0], ev.exons2[0, 0]) + 1, ev.exons1[-1, 0] + 1)
                    schain = '%i^,%i^' % (min(ev.exons1[0, 1], ev.exons2[0, 1]), max(ev.exons1[0, 1], ev.exons2[0, 1]))
                else:
                    raise Exception("Misconfigured alt-prime event detected")
            elif ev.event_type == 'intron_retention':
                struc = '0,1^2-'
                flanks = '%i-,%i^' % (ev.exons2[0] + 1, ev.exons2[1])
                schain = ',%i^%i-' % (ev.exons1[0, 1], ev.exons1[1, 0] + 1)
            elif ev.event_type == 'mutex_exons':
                struc = '1-2^,3-4^'
                flanks = '%i^,%i-' % (ev.exons1[0, 1], ev.exons1[-1, 0] + 1) 
                schain = '%i-%i^,%i-%i^' % (ev.exons1[1, 0] + 1, ev.exons1[1, 1], ev.exons2[1, 0] + 1, ev.exons2[1, 1])
            elif ev.event_type == 'mult_exon_skip':
                if mult_exon_skip_bool:
                    mult_exon_skip_bool = False
                    print >> sys.stderr, 'WARNING: Event type mult_exon_skip not implemented yet for structured output'
                    break
            else:
                raise Exception("Unknown event type: %s" % ev.event_type)
        ### - strand - revert donor/acceptor
        else:
            if ev.event_type == 'exon_skip':
                struc = '0,1-2^'
                flanks = '%i^,%i-' % (ev.exons1[-1, 0] + 1, ev.exons1[0, 1])
                schain = ',%i-%i^' % (ev.exons2[1, 1], ev.exons2[1, 0] + 1)
            elif ev.event_type in ['alt_3prime', 'alt_5prime']:
                if sp.all(ev.exons1[0, :] == ev.exons2[0, :]):
                    struc = '1-,2-'
                    flanks = '%i-,%i-' % (min(ev.exons1[-1, 1], ev.exons2[-1, 1]), ev.exons1[0, 1])
                    schain = '%i^,%i^' % (max(ev.exons1[-1, 0], ev.exons2[-1, 0]) + 1, min(ev.exons1[-1, 0], ev.exons2[-1, 0]) + 1)
                elif sp.all(ev.exons1[-1, :] == ev.exons2[-1, :]):
                    struc = '1^,2^'
                    flanks = '%i^,%i^' % (ev.exons1[-1, 0] + 1, max(ev.exons1[0, 0], ev.exons2[0, 0]) + 1)
                    schain = '%i-,%i-' % (max(ev.exons1[0, 1], ev.exons2[0, 1]), min(ev.exons1[0, 1], ev.exons2[0, 1]))
                else:
                    raise Exception("Misconfigured alt-prime event detected")
            elif ev.event_type == 'intron_retention':
                struc = '0,1^2-'
                flanks = '%i-,%i^' % (ev.exons2[1], ev.exons2[0] + 1)
                schain = ',%i^%i-' % (ev.exons1[1, 0] + 1, ev.exons1[0, 1])
            elif ev.event_type == 'mutex_exons':
                struc = '1-2^,3-4^'
                flanks = '%i^,%i-' % (ev.exons1[-1, 0] + 1, ev.exons1[0, 1]) 
                schain = '%i-%i^,%i-%i^' % (ev.exons1[1, 1], ev.exons1[1, 0] + 1, ev.exons2[1, 1], ev.exons2[1, 0] + 1)
            elif ev.event_type == 'mult_exon_skip':
                if mult_exon_skip_bool:
                    mult_exon_skip_bool = False
                    print >> sys.stderr, 'WARNING: Event type mult_exon_skip not implemented yet for structured output'
                    break
            else:
                raise Exception("Unknown event type: %s" % ev.event_type)
        print >> fd_out, 'flanks "%s"; structure "%s"; splice_chain "%s";'  % (flanks, struc, schain)

    fd_out.close()


