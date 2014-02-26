def verify_all_events(ev, strain_idx=None, list_bam=None, event_type=None, CFG=None):
    # ev = verify_all_events(ev, strain_idx, list_bam, event_type, CFG) ;

    #fieldnames = {'intron', 'exon1', 'exon2', 'exon_alt1', 'exon_alt2', 'exon_const', 'intron1', 'intron2', 'exon', 'exon_pre', 'exon_aft', 'exons'};

    ### set parameters if called by rproc
    if strain_idx is None:
        PAR = ev
        ev = PAR['ev']
        strain_idx = PAR['strain_idx']
        list_bam = PAR['list_bam']
        if 'out_fn' in PAR:
            out_fn = PAR['out_fn']
        event_type = PAR['event_type']
        CFG = PAR['CFG']

    ### verify the events if demanded
    if CFG['verify_alt_events']:
        for j, s_idx in enumerate(strain_idx):
            print '%i/%i\r' % (j, len(strain_idx))
            for i in range(ev.shape[0]):
                ev_tmp = ev[i]
                ev_tmp.subset_strain(s_idx) ### TODO 
                if event_type == 'exon_skip':
                    ver, info = verify_exon_skip(ev_tmp, list_bam[s_idx], CFG)
                    ev[i].verified[j] = ver
                    ev[i].info[j] = info
                elif event_type in ['alt_3prime', 'alt_5prime']:
                    ver, info = verify_alt_prime(ev_tmp, list_bam[s_idx], CFG)
                    ev[i].verified[j] = ver
                    ev[i].info[j] = info
                elif event_type == 'intron_retention':
                    ver, info = verify_intron_retention(ev_tmp, list_bam[s_idx], CFG)
                    ev[i].verified[j] = ver
                    ev[i].info[j] = info
                elif event_type == 'mult_exon_skip':
                    ver, info = verify_mult_exon_skip(ev_tmp, list_bam[s_idx], CFG)
                    ev[i].verified[j] = ver
                    ev[i].info[j] = info

    if out_fn is not None:
        cPickle.dump(ev, open(out_fn, 'w'))

    return ev
