def verify_alt_prime(event, fn_bam, CFG):
# [verified, info] = verify_exon_skip(event, fn_bam, CFG)

if not 'is_half_open' in CFG:
    is_half_open = False
else:
    is_half_open = CFG['is_half_open']

info['exon_diff_cov'] = 0
info['exon_const_cov'] = 0
info['intron1_conf'] = 0
info['intron2_conf'] = 0
info['valid'] = 1 

verified = [0, 0]

### check validity of exon coordinates (>=0)
if sp.any(event.exon_alt1) <= 0 or sp.any(event.exon_alt2) <= 0:
    info.valid = 0 
    return (verified, info)

### check validity of intron coordinates (only one side is differing)
if (event.intron1[0] != event.intron2[0]) and (event.intron1[1] != event.intron2[1]):
    info.valid = 0 
    return (verified, info)

if is_half_open and event.strand == '-':
    event.exon_alt1 += 1
    event.exon_alt2 += 1 
    event.exon_const += 1

gg.strand = event.strand 
gg.chr = event.chr 
gg.chr_num = event.chr_num 
gg.start = sp.array([event.exon_alt1[0], event.exon_alt2[0], event.exon_const[0]).min()
gg.stop = max([event.exon_alt1[1], event.exon_alt2[1], event.exon_const[1]]).max()

### add RNA-seq evidence to the gene structure
(tracks, intron_list) = add_reads_from_bam(gg, fn_bam, ['exon_track','intron_list'], CFG['read_filter'], CFG['var_aware']);

### set offset for half open intervals
if is_half_open:
    ho_offset = 1
else
    ho_offset = 0

### compute exon coverages as mean of position wise coverage
if (event.exon_alt1[0] == event.exon_alt2[0]):
    if event.exon_alt1[1] > event.exon_alt2[1]:
        idx_diff = sp.arange(event.exon_alt[1], event.exon_alt1[1]) - gg.start - ho_offset
        idx_const = sp.arange(event.exon_alt2[0], event.exon_alt2[1] - ho_offset) - gg.start
    else:
        idx_diff = sp.arange(event.exon_alt1[1], event.exon_alt2[1]) - gg.start - ho_offset
        idx_const = sp.arange(event.exon_alt1[0], event.exon_alt1[1] - ho_offset) - gg.start
elif (event.exon_alt1[1] == event.exon_alt2[1]):
    if event.exon_alt1[0] < event.exon_alt2[0]:
        idx_diff = sp.arange(event.exon_alt1[0], event.exon_alt2[0]) - gg.start
        idx_const = sp.arange(event.exon_alt2[0], event.exon_alt2[1] - ho_offset) - gg.start
    else:
        idx_diff = sp.arange(event.exon_alt2[0], event.exon_alt1[0]] - gg.start
        idx_const = sp.arange(event.exon_alt1[0], event.exon_alt1[1] - ho_offset) - gg.start

exon_diff_coverage = sp.mean(sp.sum(tracks[:, idx_diff], axis=0)) 
exon_const_coverage = sp.mean(sp.sum(tracks[:, idx_const], axis=0)) 

info['exon_diff_cov'] = exon_diff_coverage 
info['exon_const_cov'] = exon_const_coverage 

if exon_diff_coverage >= CFG['alt_prime']['min_diff_rel_cov'] * exon_const_coverage:
    verified[0] = 1

### check intron confirmations as sum of valid intron scores
### intron score is the number of reads confirming this intron
intron_tol = 0
idx = sp.where((abs(intron_list[0,:] - event.intron1[0]) <= intron_tol) & (abs(intron_list[1,:] - (event.intron1[1] - ho_offset)) <= intron_tol))[0]
if idx.shape[0] > 0:
    info['intron1_conf'] = sp.sum(intron_list[2, idx])

idx = sp.where((abs(intron_list[0,:] - event.intron2[0]) <= intron_tol) & (abs(intron_list[1,:] - (event.intron2[1] - ho_offset)) <= intron_tol))[0]
if idx.shape[0] > 0:
    info['intron2_conf'] = sp.sum(intron_list[2, idx])

if min(info['intron1_conf'], info['intron2_conf']) >= CFG['alt_prime']['min_intron_count']:
    verified[1] = 1

if is_half_open and event.strand == '-':
    event.exon_alt1 -= 1
    event.exon_alt2 -= 1
    event.exon_const -= 1
