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

