def verify_mult_exon_skip(event, fn_bam, CFG):
    # [verified, info] = verify_mult_exon_skip(event, fn_bam, CFG) 

    verified = [0, 0, 0, 0, 0]

    info['exons_cov'] = 0
    info['exon_pre_cov'] = 0
    info['exon_aft_cov'] = 0
    info['exon_pre_exon_conf'] = 0
    info['exon_exon_aft_conf'] = 0
    info['exon_pre_exon_aft_conf'] = 0
    info['sum_inner_exon_conf'] = 0
    info['num_inner_exon'] = 0
    info['valid'] = 1

    ### check validity of exon coordinates (>=0)
    if sp.any(event.exons1 < 0) or sp.any(event.exons2 < 0):
        info['valid'] = False
        return (verified, info)
    ### check validity of exon coordinates (start < stop && non-overlapping)
    elif sp.any(event.exons1[:, 1] - event.exons1[:, 0] < 1) or sp.any(event.exons2[:, 1] - event.exons2[:, 0] < 1):
        info['valid'] = False
        return (verified, info)

    gg = Gene()
    gg.strand = event.strand
    gg.chr = event.chr
    gg.chr_num = event.chr_num ;
    gg.start = event.exons1[0, 0]
    gg.stop = event.exon1[-1, -1]
    assert(event.exons1[0, 0] == event.exons2[0, 0])
    assert(event.exons1[-1, -1] == event.exons2[-1, -1])

    ### add RNA-seq evidence to the gene structure
    (tracks, intron_list) = add_reads_from_bam(gg, fn_bam, ['exon_track','intron_list'], CFG['read_filter'], CFG['var_aware']);

    ### compute exon coverages as mean of position wise coverage
    idx = sp.arange(event.exons2[0, 0], event.exons2[0, 1]) - gg.start
    exon_coverage_pre = sp.mean(sp.sum(tracks[:, idx], axis=0))

    idx = sp.arange(event.exons2[-1, 0], event.exons2[-1, 1]) - gg.start
    exon_coverage_aft = sp.mean(sp.sum(tracks[:, idx], axis=0))

    idx = sp.array([i for sublist in [range(event.exons2[x, 0], event.exons2[x, 1]) for x in range(1, event.exons2.shape[0] - 1)] for i in sublist]) - gg.start
    exons_coverage = sp.mean(sp.sum(tracks[:, idx], axis=0))

    info['exon_pre_cov'] = exon_coverage_pre
    info['exon_aft_cov'] = exon_coverage_aft
    info['exons_cov'] = exons_coverage

    ### check if coverage of skipped exon is >= than FACTOR times average of pre and after
    if exons_coverage >= CFG['mult_exon_skip']['min_skip_rel_cov'] * (exon_coverage_pre + exon_coverage_aft) / 2 
        verified[0] = 1 

    ### check intron confirmation as sum of valid intron scores
    ### intron score is the number of reads confirming this intron
    intron_tol = 0
    idx = sp.where((sp.absolute(intron_list[:, 0] - event.exons2[0, 1]) <= intron_tol) & (sp.absolute(intron_list[:, 1] - event.exons2[1, 0]) <= intron_tol))[0]
    if idx.shape[0] > 0:
        info['exon_pre_exon_conf'] = sp.sum(intron_list[idx, 2])
    idx = sp.where((sp.absolute(intron_list[:, 0] - event.exons2[-2, 1]) <= intron_tol) & (sp.absolute(intron_list[:, 1] - event.exons2[-1, 0]) <= intron_tol))[0]
    if idx.shape[0] > 0:
        info['exon_exon_aft_conf'] = sp.sum(intron_list[idx, 2])
    idx = sp.where((sp.absolute(intron_list[:, 0] - event.exons2[0, 1]) <= intron_tol) & (sp.absolute(intron_list[:, 1] - event.exons2[-1, 0]) <= intron_tol))[0]
    if idx.shape[0] > 0:
        info['exon_pre_exon_aft_conf'] = sp.sum(intron_list[idx, 2])
    for i in range(1, event.exons2.shape[0] - 2):
        idx = sp.where((sp.absolute(intron_list[:, 0] - event.exons2[i, 1]) <= intron_tol) & (sp.absolute(intron_list[:, 1] - event.exons2[i + 1, 0]) <= intron_tol))[0]
        if idx.shape[0] > 0:
            info['sum_inner_exon_conf'] += sp.sum(intron_list[idx, 2])
    info['num_inner_exon'] = event.exons2.shape[0] - 2
    if info['exon_pre_exon_conf'] >= CFG['mult_exon_skip']['min_non_skip_count']:
        verified[1] = 1 
    if info['exon_exon_aft_conf'] >= CFG['mult_exon_skip']['min_non_skip_count']:
        verified[2] = 1 
    if (info['sum_inner_exon_conf'] / float(info['num_inner_exon'])) >= CFG['mult_exon_skip']['min_non_skip_count']:
        verified[3] = 1 
    if info['exon_pre_exon_aft_conf'] >= CFG['mult_exon_skip']['min_skip_count']:
        verified[4] = 1 

    return (verified, info)


def verify_intron_retention(event, fn_bam, CFG):
    # [verified, info] = verify_intron_retention(event, fn_bam, CFG)

    verified = [0, 0]

    info['intron_cov'] = 0
    info['intron_cov_region'] = 0
    info['exon1_cov'] = 0
    info['exon2_cov'] = 0
    info['intron_conf'] = 0
    info['valid'] = True

    ### check validity of exon coordinates (>=0)
    if sp.any(event.exons1 < 0) or sp.any(event.exons2 < 0):
        info['valid'] = False
        return (verified, info)
    ### check validity of exon coordinates (start < stop && non-overlapping)
    elif sp.any(event.exons1[:, 1] - event.exons1[:, 0] < 1) or sp.any(event.exons2[:, 1] - event.exons2[:, 0] < 1):
        info['valid'] = False
        return (verified, info)

    gg = Gene()
    gg.strand = event.strand
    gg.chr = event.chr
    gg.chr_num = event.chr_num
    gg.start = event.exons1[0, 0]
    gg.stop = event.exons1[-1, -1]
    assert(event.exons1[0, 0] == event.exons2[0, 0])
    assert(event.exons1[-1, -1] == event.exons2[-1, -1])

    ### add RNA-seq evidence to the gene structure
    (tracks, intron_list) = add_reads_from_bam(gg, fn_bam, ['exon_track','intron_list'], CFG['read_filter'], CFG['var_aware']);

    ### compute exon coverages as mean of position wise coverage
    idx = sp.arange(event.exons1[0, 0], event.exons1[0, 1]) - gg.start
    exon_coverage1 = sp.mean(sp.sum(tracks[:, idx], axis=0))

    idx = sp.arange(event.exons1[1, 0], event.exons1[1, 1]) - gg.start
    exon_coverage2 = sp.mean(sp.sum(tracks[:, idx], axis=0))

    idx = sp.arange(event.exons1[0, 1], event.exons1[1, 0]) - gg.start
    icov = sp.sum(tracks[:, idx], axis=0)

    info['intron_cov'] = sp.mean(icov)
    info['intron_cov_region'] = sp.mean(icov > 0)
    info['exon1_cov'] = exon_coverage1
    info['exon2_cov'] = exon_coverage2

    ### check if counts match verification criteria
    if sp.mean(icov) > CFG['intron_retention']['min_retention_cov'] and
       sp.mean(icov > 0) > CFG['intron_retention']['min_retention_region'] and
       sp.mean(icov) >= (CFG['intron_retention']['min_retention_rel_cov'] * (exon_coverage1 + exon_coverage2) / 2): 
        verified[0] = 1

    ### check intron confirmation as sum of valid intron scores
    ### intron score is the number of reads confirming this intron
    intron_tol = 0
    idx = sp.where((sp.absolute(intron_list[:, 0] - event.exons1[0, 1]) <= intron_tol) & (sp.absolute(intron_list[:, 1] - event.exons1[1, 0]) <= intron_tol))[0]
    if idx.shape[0] > 0:
        info['intron_conf'] = sp.sum(intron_list[idx, 2])
        if info['intron_conf'] >= CFG['intron_retention']['min_non_retention_count']:
            verified[1] = 1

    return (verified, info)



def verify_exon_skip(event, fn_bam, CFG):
    # [verified, info] = verify_exon_skip(event, fn_bam, CFG)

    verified = [0, 0, 0, 0]

    info['exon_cov'] = 0
    info['exon_pre_cov'] = 0
    info['exon_aft_cov'] = 0 
    info['exon_pre_exon_conf'] = 0 
    info['exon_exon_aft_conf'] = 0 
    info['exon_pre_exon_aft_conf'] = 0 
    info['valid'] = True

    ### check validity of exon coordinates (>=0)
    if sp.any(event.exons1 < 0) or sp.any(event.exons2 < 0):
        info['valid'] = False
        return (verified, info)
    ### check validity of exon coordinates (start < stop && non-overlapping)
    elif sp.any(event.exons1[:, 1] - event.exons1[:, 0] < 1) or sp.any(event.exons2[:, 1] - event.exons2[:, 0] < 1):
        info['valid'] = False
        return (verified, info)

    gg = Gene()
    gg.strand = event.strand
    gg.chr = event.chr
    gg.chr_num = event.chr_num
    gg.start = event.exons1[0, 0]
    gg.stop = event.exon1[-1, -1]
    assert(event.exons1[0, 0] == event.exons2[0, 0])
    assert(event.exons1[-1, -1] == event.exons2[-1, -1])

    ### add RNA-seq evidence to the gene structure
    (tracks, intron_list) = add_reads_from_bam(gg, fn_bam, ['exon_track','intron_list'], CFG['read_filter'], CFG['var_aware']);

    ### compute exon coverages as mean of position wise coverage
    idx = sp.arange(exons2[0, 0], exons2[0, 1]) - gg.start
    exon_coverage_pre = sp.mean(sp.sum(tracks[:, idx], axis=0))

    idx = sp.arange(exons2[2, 0], exons2[2, 1]) - gg.start
    exon_coverage_aft = sp.mean(sp.sum(tracks[:, idx], axis=0))

    idx = sp.arange(exons2[1, 0], exons2[1, 1]) - gg.start
    exon_coverage = sp.mean(sp.sum(tracks[:, idx], axis=0))

    info['exon_pre_cov'] = exon_coverage_pre
    info['exon_aft_cov'] = exon_coverage_aft
    info['exon_cov'] = exon_coverage

    ### check if coverage of skipped exon is >= than FACTOR times average of pre and after
    if exon_coverage >= CFG['exon_skip']['min_skip_rel_cov'] * (exon_coverage_pre + exon_coverage_aft)/2 
        verified[0] = 1

    ### check intron confirmation as sum of valid intron scores
    ### intron score is the number of reads confirming this intron
    idx = sp.where((sp.absolute(intron_list[:, 0] - event.exons2[0, 1]) <= CFG['exon_skip']['intron_tolerance']) & (sp.absolute(intron_list[:, 1] - event.exons2[1, 0]) <= CFG['exon_skip']['intron_tolerance']))[0]
    if idx.shape[0] > 0:
        info['exon_pre_exon_conf'] = sp.sum(intron_list[idx, 2])
    idx = sp.where((sp.absolute(intron_list[:, 0] - event.exons2[1, 1]) <= CFG['exon_skip']['intron_tolerance']) & (sp.absolute(intron_list[:, 1] - event.exons2[2, 0]) <= CFG['exon_skip']['intron_tolerance']))[0]
    if idx.shape[0] > 0:
        info['exon_exon_aft_conf'] = sp.sum(intron_list[idx, 2])
    idx = sp.where((sp.absolute(intron_list[:, 0] - event.exons2[0, 1]) <= CFG['exon_skip']['intron_tolerance']) & (sp.absolute(intron_list[:, 1] - event.exons2[2, 0]) <= CFG['exon_skip']['intron_tolerance'])))[0]
    if idx.shape[0] > 0:
        info['exon_pre_exon_aft_conf'] = sp.sum(intron_list[idx, 2])
    if info['exon_pre_exon_conf'] >= CFG['exon_skip']['min_non_skip_count']:
        verified[1] = 1
    if info['exon_exon_aft_conf'] >= CFG['exon_skip']['min_non_skip_count']:
        verified[2] = 1
    if info['exon_pre_exon_aft_conf'] >= CFG['exon_skip']['min_skip_count']:
        verified[3] = 1

    return (verified, info)


def verify_alt_prime(event, fn_bam, cfg):
    # [verified, info] = verify_exon_skip(event, fn_bam, cfg)

    info['exon_diff_cov'] = 0
    info['exon_const_cov'] = 0
    info['intron1_conf'] = 0
    info['intron2_conf'] = 0
    info['valid'] = 1 

    verified = [0, 0]

    ### check validity of exon coordinates (>=0)
    if sp.any(event.exons1 < 0) or sp.any(event.exons2 < 0):
        info['valid'] = 0 
        return (verified, info)

    ### check validity of intron coordinates (only one side is differing)
    if (event.exons1[0, 1] != event.exons2[0, 1]) and (event.exons1[1, 0] != event.exons2[1, 0]):
        info['valid'] = 0 
        return (verified, info)

    gg = gene()
    gg.strand = event.strand 
    gg.chr = event.chr 
    gg.chr_num = event.chr_num 
    gg.start = sp.array([event.exon_alt1[0], event.exon_alt2[0], event.exon_const[0]]).min()
    gg.stop = max([event.exon_alt1[1], event.exon_alt2[1], event.exon_const[1]]).max()

    ### add rna-seq evidence to the gene structure
    (tracks, intron_list) = add_reads_from_bam(gg, fn_bam, ['exon_track','intron_list'], cfg['read_filter'], cfg['var_aware']);

    ### compute exon coverages as mean of position wise coverage
    if (event.exon_alt1[0] == event.exon_alt2[0]):
        if event.exon_alt1[1] > event.exon_alt2[1]:
            idx_diff = sp.arange(event.exon_alt[1], event.exon_alt1[1]) - gg.start
            idx_const = sp.arange(event.exon_alt2[0], event.exon_alt2[1]) - gg.start
        else:
            idx_diff = sp.arange(event.exon_alt1[1], event.exon_alt2[1]) - gg.start
            idx_const = sp.arange(event.exon_alt1[0], event.exon_alt1[1]) - gg.start
    elif (event.exon_alt1[1] == event.exon_alt2[1]):
        if event.exon_alt1[0] < event.exon_alt2[0]:
            idx_diff = sp.arange(event.exon_alt1[0], event.exon_alt2[0]) - gg.start
            idx_const = sp.arange(event.exon_alt2[0], event.exon_alt2[1]) - gg.start
        else:
            idx_diff = sp.arange(event.exon_alt2[0], event.exon_alt1[0]) - gg.start
            idx_const = sp.arange(event.exon_alt1[0], event.exon_alt1[1]) - gg.start

    exon_diff_coverage = sp.mean(sp.sum(tracks[:, idx_diff], axis=0)) 
    exon_const_coverage = sp.mean(sp.sum(tracks[:, idx_const], axis=0)) 

    info['exon_diff_cov'] = exon_diff_coverage 
    info['exon_const_cov'] = exon_const_coverage 

    if exon_diff_coverage >= cfg['alt_prime']['min_diff_rel_cov'] * exon_const_coverage:
        verified[0] = 1

    ### check intron confirmations as sum of valid intron scores
    ### intron score is the number of reads confirming this intron
    intron_tol = 0
    idx = sp.where((abs(intron_list[0,:] - event.exons1[0, 1]) <= intron_tol) & (abs(intron_list[1,:] - event.exons1[1, 0]) <= intron_tol))[0]
    if idx.shape[0] > 0:
        info['intron1_conf'] = sp.sum(intron_list[2, idx])

    idx = sp.where((abs(intron_list[0,:] - event.exons2[0, 1]) <= intron_tol) & (abs(intron_list[1,:] - event.exons2[1, 0]) <= intron_tol))[0]
    if idx.shape[0] > 0:
        info['intron2_conf'] = sp.sum(intron_list[2, idx])

    if min(info['intron1_conf'], info['intron2_conf']) >= cfg['alt_prime']['min_intron_count']:
        verified[1] = 1

    return (verified, info)


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
