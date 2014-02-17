import sys

def write_events_gff3(fn_out_gff3, events):
# function write_events_gff3(fn_out_gff3, events)
    
    if events.shape[0] == 0:
        print >> sys.stderr, 'WARNING: No events present.'
        return

    print 'writing %s events in gff3 format to %s' % (events[0].event_type, fn_out_gff3)

    fd_out = open(fn_out_gff3, 'w+') 
    print >> fd_out, '##gff-version 3'

    ### load gene structure
    for ev = in range(events.shape[0]):

        ev = events[i]
        gene_name = events[i].gene_name[0] ### TODO - why only first?
        start_pos = ev.exons1[0, 0]
        stop_pos = ev.exons1[-1, -1]

        ### get order of isoforms o_idx(1) -> iso1 and o_idx(2) -> iso2
        ### assert that first isoform is always the shorter one
        assert(sp.sum(ev.exons1[:, 1] - ev.exons1[:, 0]) < sp.sum(ev.exons2[:, 1] - ev.exons2[:, 0]))

        name = '%s.%i' % (ev.event_type, ev.id)

        print >> fd_out, '%s\t%s\tgene\t%i\t%i\t.\t%c\t.\tID=%s;GeneName="%s"' % (ev.chr, ev.event_type, start_pos, stop_pos, ev.strand, name, ev.gene_name[0])
        print >> fd_out, '%s\t%s\tmRNA\t%i\t%i\t.\t%c\t.\tID=%s_iso1;Parent=%s;GeneName="%s"' % (ev.chr, ev.event_type, start_pos, stop_pos, ev.strand, name, name, ev.gene_name[0]) 
        for i in range(ev.exons1.shape[0]):
            print >> fd_out, '%s\t%s\texon\t%i\t%i\t.\t%c\t.\tParent=%s_iso1' % (ev.chr, ev.event_type, ev.exons1[i, 0], exons1[i, 1], ev.strand, name)
        print >> fd_out, '%s\t%s\tmRNA\t%i\t%i\t.\t%c\t.\tID=%s_iso2;Parent=%s;GeneName="%s"' % (ev.chr, ev.event_type, start_pos, stop_pos, ev.strand, name, name, ev.gene_name[0])
        for i in range(ev.exons1.shape[0]):
            print >> fd_out, '%s\t%s\texon\t%i\t%i\t.\t%c\t.\tParent=%s_iso2' % (ev.chr, ev.event_type, ev.exons2[i, 0], exons2[i, 1], ev.strand, name)
    fd_out.close()
