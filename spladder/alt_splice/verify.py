import scipy as sp
import pickle
import h5py
import sys
import os

if __package__ is None:
    __package__ = 'modules.alt_splice'

from ..utils import *

def verify_mult_exon_skip(event, gene, counts_segments, counts_edges, options):

    verified = [0, 0, 0, 0, 0]
    # (0) exon coordinates are valid (>= 0 && start < stop && non-overlapping) & skipped exon coverage >= FACTOR * mean(pre, after)
    # (1) inclusion count first intron >= threshold
    # (2) inclusion count last intron >= threshold
    # (3) avg inclusion count inner exons >= threshold
    # (4) skip count >= threshold

    info = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    # (0) valid, (1) exon_pre_cov, (2) exons_cov, (3) exon_aft_cov
    # (4) exon_pre_exon_conf, (5) exon_exon_aft_conf, (6) exon_pre_exon_aft_conf
    # (7) sum_inner_exon_conf, (8) num_inner_exon, (9) len_inner_exon

    ### check validity of exon coordinates (>=0)
    if sp.any(event.exons1 < 0) or sp.any(event.exons2 < 0):
        info[0] = 0
        return (verified, info)
    ### check validity of exon coordinates (start < stop && non-overlapping)
    elif sp.any(event.exons1[:, 1] - event.exons1[:, 0] < 1) or sp.any(event.exons2[:, 1] - event.exons2[:, 0] < 1):
        info[0] = 0
        return (verified, info)

    sg = gene.splicegraph
    segs = gene.segmentgraph

    ### find exons corresponding to event
    idx_exon_pre  = sp.where((sg.vertices[0, :] == event.exons2[0, 0]) & (sg.vertices[1, :] == event.exons2[0, 1]))[0]
    idx_exon_aft  = sp.where((sg.vertices[0, :] == event.exons2[-1, 0]) & (sg.vertices[1, :] == event.exons2[-1, 1]))[0]
    seg_exons = []
    for i in range(1, event.exons2.shape[0] - 1):
        tmp = sp.where((sg.vertices[0, :] == event.exons2[i, 0]) & (sg.vertices[1, :] == event.exons2[i, 1]))[0]
        seg_exons.append(sp.where(segs.seg_match[tmp, :])[1])
    
    ### find segments corresponding to exons
    seg_exon_pre = sp.sort(sp.where(segs.seg_match[idx_exon_pre, :])[1])
    seg_exon_aft = sp.sort(sp.where(segs.seg_match[idx_exon_aft, :])[1])
    seg_exons_u = sp.sort(sp.unique([x for sublist in seg_exons for x in sublist]))

    seg_lens = segs.segments[1, :] - segs.segments[0, :]

    # exon_pre_cov
    info[1] = sp.sum(counts_segments[seg_exon_pre] * seg_lens[seg_exon_pre]) / sp.sum(seg_lens[seg_exon_pre])
    # exon_aft_cov
    info[3] = sp.sum(counts_segments[seg_exon_aft] * seg_lens[seg_exon_aft]) / sp.sum(seg_lens[seg_exon_aft])
    # exons_cov
    info[2] = sp.sum(counts_segments[seg_exons_u] * seg_lens[seg_exons_u]) / sp.sum(seg_lens[seg_exons_u])

    ### check if coverage of skipped exon is >= than FACTOR times average of pre and after
    if info[2] >= options.mult_exon_skip['min_skip_rel_cov'] * (info[1] + info[3]) / 2:
        verified[0] = 1

    ### check intron confirmation as sum of valid intron scores
    ### intron score is the number of reads confirming this intron
    # exon_pre_exon_conf
    idx = sp.where(counts_edges[:, 0] == sp.ravel_multi_index([seg_exon_pre[-1], seg_exons[0][0]], segs.seg_edges.shape))[0]
    if len(idx.shape) > 0 and idx.shape[0] > 0:
        info[4] = counts_edges[idx[0], 1]
    # exon_exon_aft_conf
    idx = sp.where(counts_edges[:, 0] == sp.ravel_multi_index([seg_exons[-1][-1], seg_exon_aft[0]], segs.seg_edges.shape))[0]
    if len(idx.shape) > 0 and idx.shape[0] > 0:
        info[5] = counts_edges[idx[0], 1]
    # exon_pre_exon_aft_conf
    idx = sp.where(counts_edges[:, 0] == sp.ravel_multi_index([seg_exon_pre[-1], seg_exon_aft[0]], segs.seg_edges.shape))[0]
    if len(idx.shape) > 0 and idx.shape[0] > 0:
        info[6] = counts_edges[idx[0], 1]
    for i in range(len(seg_exons) - 1):
        # sum_inner_exon_conf
        idx = sp.where(counts_edges[:, 0] == sp.ravel_multi_index([seg_exons[i][-1], seg_exons[i+1][0]], segs.seg_edges.shape))[0]
        if len(idx.shape) > 0 and idx.shape[0] > 0:
            info[7] += counts_edges[idx[0], 1]

    # num_inner_exon
    info[8] = event.exons2.shape[0] - 2
    info[9] = sp.sum(event.exons2[1:-1, 1] - event.exons2[1:-1, 0])
    if info[4] >= options.mult_exon_skip['min_non_skip_count']:
        verified[1] = 1
    if info[5] >= options.mult_exon_skip['min_non_skip_count']:
        verified[2] = 1
    if (info[7] / info[8]) >= options.mult_exon_skip['min_non_skip_count']:
        verified[3] = 1 
    if info[6] >= options.mult_exon_skip['min_skip_count']:
        verified[4] = 1 

    return (verified, info)


def verify_intron_retention(event, gene, counts_segments, counts_edges, counts_seg_pos, options):

    verified = [0, 0]
    # (0) counts meet criteria for min_retention_cov, min_retention_region and min_retetion_rel_cov 
    # (1) min_non_retention_count >= threshold

    info = [1, 0, 0, 0, 0, 0]
    # (0) valid, (1) intron_cov, (2) exon1_cov, (3), exon2_cov
    # (4) intron_conf, (5) intron_cov_region

    ### check validity of exon coordinates (>=0)
    if sp.any(event.exons1 < 0) or sp.any(event.exons2 < 0):
        info[0] = 0
        return (verified, info)
    ### check validity of exon coordinates (start < stop && non-overlapping)
    elif sp.any(event.exons1[:, 1] - event.exons1[:, 0] < 1) or sp.any((event.exons2[1] - event.exons2[0]) < 1):
        info[0] = 0
        return (verified, info)

    sg = gene.splicegraph
    segs = gene.segmentgraph

    ### find exons corresponding to event
    idx_exon1  = sp.where((sg.vertices[0, :] == event.exons1[0, 0]) & (sg.vertices[1, :] == event.exons1[0, 1]))[0]
    idx_exon2  = sp.where((sg.vertices[0, :] == event.exons1[1, 0]) & (sg.vertices[1, :] == event.exons1[1, 1]))[0]

    ### find segments corresponding to exons
    seg_exon1 = sp.sort(sp.where(segs.seg_match[idx_exon1, :])[1])
    seg_exon2 = sp.sort(sp.where(segs.seg_match[idx_exon2, :])[1])
    seg_all = sp.arange(seg_exon1[0], seg_exon2[-1])

    seg_intron = sp.setdiff1d(seg_all, seg_exon1)
    seg_intron = sp.setdiff1d(seg_intron, seg_exon2)
    assert(seg_intron.shape[0] > 0)

    seg_lens = segs.segments[1, :] - segs.segments[0, :]

    ### compute exon coverages as mean of position wise coverage
    # exon1_cov
    info[2] = sp.sum(counts_segments[seg_exon1] * seg_lens[seg_exon1]) / sp.sum(seg_lens[seg_exon1])
    # exon2_cov
    info[3] = sp.sum(counts_segments[seg_exon2] * seg_lens[seg_exon2]) / sp.sum(seg_lens[seg_exon2])
    # intron_cov
    info[1] = sp.sum(counts_segments[seg_intron] * seg_lens[seg_intron]) / sp.sum(seg_lens[seg_intron])
    # intron_cov_region
    info[5] = sp.sum(counts_seg_pos[seg_intron]) / sp.sum(seg_lens[seg_intron])

    ### check if counts match verification criteria
    if info[1] > options.intron_retention['min_retention_cov'] and \
       info[5] > options.intron_retention['min_retention_region'] and \
       info[1] >= options.intron_retention['min_retention_rel_cov'] * (info[2] + info[3]) / 2:
        verified[0] = 1

    ### check intron confirmation as sum of valid intron scores
    ### intron score is the number of reads confirming this intron
    # intron conf
    idx = sp.where(counts_edges[:, 0] == sp.ravel_multi_index([seg_exon1[-1], seg_exon2[0]], segs.seg_edges.shape))[0]
    info[4] = counts_edges[idx, 1]

    if info[4] >= options.intron_retention['min_non_retention_count']:
        verified[1] = 1

    return (verified, info)


def verify_exon_skip(event, gene, counts_segments, counts_edges, options):

    verified = [0, 0, 0, 0]
    # (0) coverage of skipped exon is >= than FACTOR * mean(pre, after)
    # (1) inclusion count of first intron >= threshold 
    # (2) inclusion count of second intron >= threshold
    # (3) skip count of exon >= threshold

    info = [1, 0, 0, 0, 0, 0, 0]
    # (0) valid, (1) exon_cov, (2) exon_pre_cov, (3) exon_aft_cov, 
    # (4) exon_pre_exon_conf, (5) exon_exon_aft_conf, (6) exon_pre_exon_aft_conf

    ### check validity of exon coordinates (>=0)
    if sp.any(event.exons1 < 0) or sp.any(event.exons2 < 0):
        info[0] = False
        return (verified, info)
    ### check validity of exon coordinates (start < stop && non-overlapping)
    elif sp.any(event.exons1[:, 1] - event.exons1[:, 0] < 1) or sp.any(event.exons2[:, 1] - event.exons2[:, 0] < 1):
        info[0] = False
        return (verified, info)

    sg = gene.splicegraph
    segs = gene.segmentgraph

    ### find exons corresponding to event
    idx_exon_pre = sp.where((sg.vertices[0, :] == event.exons2[0, 0]) & (sg.vertices[1, :] == event.exons2[0, 1]))[0]
    idx_exon = sp.where((sg.vertices[0, :] == event.exons2[1, 0]) & (sg.vertices[1, :] == event.exons2[1, 1]))[0]
    idx_exon_aft = sp.where((sg.vertices[0, :] == event.exons2[2, 0]) & (sg.vertices[1, :] == event.exons2[2, 1]))[0]

    ### find segments corresponding to exons
    seg_exon_pre = sp.sort(sp.where(segs.seg_match[idx_exon_pre, :])[1])
    seg_exon_aft = sp.sort(sp.where(segs.seg_match[idx_exon_aft, :])[1])
    seg_exon = sp.sort(sp.where(segs.seg_match[idx_exon, :])[1])

    seg_lens = segs.segments[1, :] - segs.segments[0, :]

    # exon pre cov
    info[2] = sp.sum(counts_segments[seg_exon_pre] * seg_lens[seg_exon_pre]) /sp.sum(seg_lens[seg_exon_pre])
    # exon aft cov
    info[3] = sp.sum(counts_segments[seg_exon_aft] * seg_lens[seg_exon_aft]) /sp.sum(seg_lens[seg_exon_aft])
    # exon cov
    info[1] = sp.sum(counts_segments[seg_exon] * seg_lens[seg_exon]) /sp.sum(seg_lens[seg_exon])

    ### check if coverage of skipped exon is >= than FACTOR times average of pre and after
    if info[1] >= options.exon_skip['min_skip_rel_cov'] * (info[2] + info[3]) / 2: 
        verified[0] = 1

    ### check intron confirmation as sum of valid intron scores
    ### intron score is the number of reads confirming this intron
    # exon_pre_exon_conf
    idx = sp.where(counts_edges[:, 0] == sp.ravel_multi_index([seg_exon_pre[-1], seg_exon[0]], segs.seg_edges.shape))[0]
    info[4] = counts_edges[idx, 1]
    if info[4] >= options.exon_skip['min_non_skip_count']:
        verified[1] = 1
    # exon_exon_aft_conf
    idx = sp.where(counts_edges[:, 0] == sp.ravel_multi_index([seg_exon[-1], seg_exon_aft[0]], segs.seg_edges.shape))[0]
    info[5] = counts_edges[idx, 1]
    if info[5] >= options.exon_skip['min_non_skip_count']:
        verified[2] = 1
    # exon_pre_exon_aft_conf
    idx = sp.where(counts_edges[:, 0] == sp.ravel_multi_index([seg_exon_pre[-1], seg_exon_aft[0]], segs.seg_edges.shape))[0]
    info[6] = counts_edges[idx, 1]
    if info[6] >= options.exon_skip['min_skip_count']:
        verified[3] = 1

    return (verified, info)


def verify_alt_prime(event, gene, counts_segments, counts_edges, options):
    # [verified, info] = verify_exon_skip(event, fn_bam, cfg)

    verified = [0, 0]
    # (0) coverage of diff region is at least FACTOR * coverage constant region
    # (1) both alternative introns are >= threshold 

    info = [1, 0, 0, 0, 0]
    # (0) valid, (1) exon_diff_cov, (2) exon_const_cov
    # (3) intron1_conf, (4) intron2_conf

    ### check validity of exon coordinates (>=0)
    if sp.any(event.exons1 < 0) or sp.any(event.exons2 < 0):
        info[0] = 0 
        return (verified, info)

    ### check validity of intron coordinates (only one side is differing)
    if (event.exons1[0, 1] != event.exons2[0, 1]) and (event.exons1[1, 0] != event.exons2[1, 0]):
        info[0] = 0 
        return (verified, info)

    sg = gene.splicegraph
    segs = gene.segmentgraph

    ### find exons corresponding to event
    idx_exon11 = sp.where((sg.vertices[0, :] == event.exons1[0, 0]) & (sg.vertices[1, :] == event.exons1[0, 1]))[0]
    if idx_exon11.shape[0] == 0:
        segs_exon11 = sp.where((segs.segments[0, :] >= event.exons1[0, 0]) & (segs.segments[1, :] <= event.exons1[0, 1]))[0]
    else:
        segs_exon11 = sp.where(segs.seg_match[idx_exon11, :])[1]
    idx_exon12 = sp.where((sg.vertices[0, :] == event.exons1[1, 0]) & (sg.vertices[1, :] == event.exons1[1, 1]))[0]
    if idx_exon12.shape[0] == 0:
        segs_exon12 = sp.where((segs.segments[0, :] >= event.exons1[1, 0]) & (segs.segments[1, :] <= event.exons1[1, 1]))[0]
    else:
        segs_exon12 = sp.where(segs.seg_match[idx_exon12, :])[1]
    idx_exon21 = sp.where((sg.vertices[0, :] == event.exons2[0, 0]) & (sg.vertices[1, :] == event.exons2[0, 1]))[0]
    if idx_exon21.shape[0] == 0:
        segs_exon21 = sp.where((segs.segments[0, :] >= event.exons2[0, 0]) & (segs.segments[1, :] <= event.exons2[0, 1]))[0]
    else:
        segs_exon21 = sp.where(segs.seg_match[idx_exon21, :])[1]
    idx_exon22 = sp.where((sg.vertices[0, :] == event.exons2[1, 0]) & (sg.vertices[1, :] == event.exons2[1, 1]))[0]
    if idx_exon22.shape[0] == 0:
        segs_exon22 = sp.where((segs.segments[0, :] >= event.exons2[1, 0]) & (segs.segments[1, :] <= event.exons2[1, 1]))[0]
    else:
        segs_exon22 = sp.where(segs.seg_match[idx_exon22, :] > 0)[1]

    assert(segs_exon11.shape[0] > 0)
    assert(segs_exon12.shape[0] > 0)
    assert(segs_exon21.shape[0] > 0)
    assert(segs_exon22.shape[0] > 0)

    if segs_exon11.shape == segs_exon21.shape and sp.all(segs_exon11 == segs_exon21):
        seg_exon_const = segs_exon11
        seg_diff = sp.setdiff1d(segs_exon12, segs_exon22)
        if seg_diff.shape[0] == 0:
            seg_diff = sp.setdiff1d(segs_exon22, segs_exon12)
        seg_const = sp.intersect1d(segs_exon12, segs_exon22)
    elif segs_exon12.shape == segs_exon22.shape and sp.all(segs_exon12 == segs_exon22):
        seg_exon_const = segs_exon12
        seg_diff = sp.setdiff1d(segs_exon11, segs_exon21)
        if seg_diff.shape[0] == 0:
            seg_diff = sp.setdiff1d(segs_exon21, segs_exon11)
        seg_const = sp.intersect1d(segs_exon21, segs_exon11)
    else:
        print("ERROR: both exons differ in alt prime event in verify_alt_prime", file=sys.stderr)
        sys.exit(1)
    seg_const = sp.r_[seg_exon_const, seg_const]

    seg_lens = segs.segments[1, :] - segs.segments[0, :]

    # exon_diff_cov
    info[1] = sp.sum(counts_segments[seg_diff] * seg_lens[seg_diff]) / sp.sum(seg_lens[seg_diff])
    # exon_const_cov
    info[2] = sp.sum(counts_segments[seg_const] * seg_lens[seg_const]) / sp.sum(seg_lens[seg_const])

    if info[1] >= options.alt_prime['min_diff_rel_cov'] * info[2]:
        verified[0] = 1

    ### check intron confirmations as sum of valid intron scores
    ### intron score is the number of reads confirming this intron
    # intron1_conf 
    idx = sp.where(counts_edges[:, 0] == sp.ravel_multi_index([segs_exon11[-1], segs_exon12[0]], segs.seg_edges.shape))[0]
    assert(idx.shape[0] > 0)
    info[3] = counts_edges[idx, 1]
    # intron2_conf 
    idx = sp.where(counts_edges[:, 0] == sp.ravel_multi_index([segs_exon21[-1], segs_exon22[0]], segs.seg_edges.shape))[0]
    assert(idx.shape[0] > 0)
    info[4] = counts_edges[idx, 1]

    if min(info[3], info[4]) >= options.alt_prime['min_intron_count']:
        verified[1] = 1

    return (verified, info)

def verify_mutex_exons(event, gene, counts_segments, counts_edges, options):
    #

    verified = [0, 0, 0, 0]
    # (0) coverage of first alt exon is >= than FACTOR times average of pre and after 
    # (1) coverage of second alt exon is >= than FACTOR times average of pre and after 
    # (2) both introns neighboring first alt exon are confirmed >= threshold
    # (3) both introns neighboring second alt exon are confirmed >= threshold

    info = [1, 0, 0, 0, 0, 0, 0, 0, 0]
    # (0) valid, (1) exon_pre_cov, (2) exon1_cov, (3) exon2_cov, (4) exon_aft_cov, 
    # (5) exon_pre_exon1_conf, (6) exon_pre_exon2_conf, (7) exon1_exon_aft_conf, (8) exon2_exon_aft_conf

    ### check validity of exon coordinates (>=0)
    if sp.any(event.exons1 < 0) or sp.any(event.exons2 < 0):
        info[0] = 0
        return (verified, info)
    ### check validity of exon coordinates (start < stop && non-overlapping)
    elif sp.any(event.exons1[:, 1] - event.exons1[:, 0] < 1) or sp.any(event.exons2[:, 1] - event.exons2[:, 0] < 1) or \
         (event.exons1[1, 1] > event.exons2[1, 0] and event.exons1[1, 0] < event.exons2[1, 0]) or \
         (event.exons2[1, 1] > event.exons1[1, 0] and event.exons2[1, 0] < event.exons1[1, 0]):
        info[0] = 0
        return (verified, info)

    sg = gene.splicegraph
    segs = gene.segmentgraph

    ### find exons corresponding to event
    idx_exon_pre  = sp.where((sg.vertices[0, :] == event.exons1[0, 0]) & (sg.vertices[1, :] == event.exons1[0, 1]))[0]
    idx_exon_aft  = sp.where((sg.vertices[0, :] == event.exons1[-1, 0]) & (sg.vertices[1, :] == event.exons1[-1, 1]))[0]
    idx_exon1  = sp.where((sg.vertices[0, :] == event.exons1[1, 0]) & (sg.vertices[1, :] == event.exons1[1, 1]))[0]
    idx_exon2  = sp.where((sg.vertices[0, :] == event.exons2[1, 0]) & (sg.vertices[1, :] == event.exons2[1, 1]))[0]
    
    ### find segments corresponding to exons
    seg_exon_pre = sp.sort(sp.where(segs.seg_match[idx_exon_pre, :])[1])
    seg_exon_aft = sp.sort(sp.where(segs.seg_match[idx_exon_aft, :])[1])
    seg_exon1 = sp.sort(sp.where(segs.seg_match[idx_exon1, :])[1])
    seg_exon2 = sp.sort(sp.where(segs.seg_match[idx_exon2, :])[1])

    seg_lens = segs.segments[1, :] - segs.segments[0, :]

    # exon pre cov
    info[1] = sp.sum(counts_segments[seg_exon_pre] * seg_lens[seg_exon_pre]) / sp.sum(seg_lens[seg_exon_pre])
    # exon1 cov
    info[2] = sp.sum(counts_segments[seg_exon1] * seg_lens[seg_exon1]) / sp.sum(seg_lens[seg_exon1])
    # exon2 cov
    info[3] = sp.sum(counts_segments[seg_exon2] * seg_lens[seg_exon2]) / sp.sum(seg_lens[seg_exon2])
    # exon aft cov
    info[4] = sp.sum(counts_segments[seg_exon_aft] * seg_lens[seg_exon_aft]) / sp.sum(seg_lens[seg_exon_aft])

    ### check if coverage of first exon is >= than FACTOR times average of pre and after
    if info[2] >= options.mutex_exons['min_skip_rel_cov'] * (info[1] + info[4])/2:
        verified[0] = 1
    if info[3] >= options.mutex_exons['min_skip_rel_cov'] * (info[1] + info[4])/2:
        verified[1] = 1

    ### check intron confirmation as sum of valid intron scores
    ### intron score is the number of reads confirming this intron
    # exon_pre_exon1_conf
    idx = sp.where(counts_edges[:, 0] == sp.ravel_multi_index([seg_exon_pre[-1], seg_exon1[0]], segs.seg_edges.shape))[0]
    if len(idx.shape) > 0 and idx.shape[0] > 0:
        info[5] = counts_edges[idx[0], 1]
    # exon_pre_exon2_conf
    idx = sp.where(counts_edges[:, 0] == sp.ravel_multi_index([seg_exon_pre[-1], seg_exon2[0]], segs.seg_edges.shape))[0]
    if len(idx.shape) > 0 and idx.shape[0] > 0:
        info[6] = counts_edges[idx[0], 1]
    # exon1_exon_aft_conf
    idx = sp.where(counts_edges[:, 0] == sp.ravel_multi_index([seg_exon1[-1], seg_exon_aft[0]], segs.seg_edges.shape))[0]
    if len(idx.shape) > 0 and idx.shape[0] > 0:
        info[7] = counts_edges[idx[0], 1]
    # exon2_exon_aft_conf
    idx = sp.where(counts_edges[:, 0] == sp.ravel_multi_index([seg_exon2[-1], seg_exon_aft[0]], segs.seg_edges.shape))[0]
    if len(idx.shape) > 0 and idx.shape[0] > 0:
        info[8] = counts_edges[idx[0], 1]

    # set verification flags for intron confirmation
    if min(info[5], info[6]) >= options.mutex_exons['min_conf_count']:
        verified[2] = 1
    if min(info[7], info[8]) >= options.mutex_exons['min_conf_count']:
        verified[3] = 1

    return (verified, info)


def verify_all_events(ev, strain_idx=None, list_bam=None, event_type=None, options=None, out_fn=None):

    ### set parameters if called by rproc
    if strain_idx is None:
        PAR = ev
        ev = PAR['ev']
        strain_idx = PAR['strain_idx']
        list_bam = PAR['list_bam']
        if 'out_fn' in PAR:
            out_fn = PAR['out_fn']
        event_type = PAR['event_type']
        options = PAR['options']

    ### verify the events if demanded
    if options.verify_alt_events:

        prune_tag = ''
        if options.do_prune:
            prune_tag = '_pruned'
        validate_tag = ''
        if options.validate_sg:
            validate_tag = '.validated'

        if options.merge == 'single':
            (genes, inserted) = pickle.load(open('%s/spladder/genes_graph_conf%i.%s%s%s.pickle' % (options.outdir, options.confidence, options.samples[strain_idx], validate_tag, prune_tag), 'rb'))
            fn_count = '%s/spladder/genes_graph_conf%i.%s%s%s.count.hdf5' % (options.outdir, options.confidence, options.samples[strain_idx], validate_tag, prune_tag)
        else:
            (genes, inserted) = pickle.load(open('%s/spladder/genes_graph_conf%i.%s%s%s.pickle' % (options.outdir, options.confidence, options.merge, validate_tag, prune_tag), 'rb'))
            fn_count = '%s/spladder/genes_graph_conf%i.%s%s%s.count.hdf5' % (options.outdir, options.confidence, options.merge, validate_tag, prune_tag)

        ### load count index data from hdf5
        IN = h5py.File(fn_count, 'r')
        if os.path.exists(fn_count + '.quick_ids_segs'):
            gene_ids_segs = pickle.load(open(fn_count + '.quick_ids_segs', 'rb'))
        else:
            gene_ids_segs = IN['gene_ids_segs'][:]
            pickle.dump(gene_ids_segs, open(fn_count + '.quick_ids_segs', 'wb'), -1)
        if os.path.exists(fn_count + '.quick_ids_edges'):
            gene_ids_edges = pickle.load(open(fn_count + '.quick_ids_edges', 'rb'))
        else:
            gene_ids_edges = IN['gene_ids_edges'][:]
            pickle.dump(gene_ids_edges, open(fn_count + '.quick_ids_edges', 'wb'), -1)
        if os.path.exists(fn_count + '.quick_edge_idx'):
            edge_idx = pickle.load(open(fn_count + '.quick_edge_idx', 'rb'))
        else:
            edge_idx = IN['edge_idx'][:]
            pickle.dump(edge_idx, open(fn_count + '.quick_edge_idx', 'wb'), -1)

        ### sort events by gene idx
        s_idx = sp.argsort([x.gene_idx for x in ev])
        ev = ev[s_idx]
        old_idx = sp.argsort(s_idx)

        counts = []
        for i in range(ev.shape[0]):
            
            sys.stdout.write('.')
            if i > 0 and i % 50 == 0:
                sys.stdout.write('%i (%i)\n' % (i, ev.shape[0]))
            sys.stdout.flush()
        
            g_idx = ev[i].gene_idx
            ev[i].verified = [] ### TODO: maybe solve that differently

            ### there are no edges present in the event
            if gene_ids_edges.shape[0] == 0:
                ver, info = verify_empty(event_type)
                counts.append(sp.array([info]))
                ev[i].verified = sp.array(ev[i].verified, dtype='bool')
                continue

            gr_idx_segs = sp.where(gene_ids_segs == g_idx)[0]
            gr_idx_edges = sp.where(gene_ids_edges == g_idx)[0]
            if gr_idx_edges.shape[0] == 0:
                ver, info = verify_empty(event_type)
                counts.append(sp.array([info]))
                ev[i].verified = sp.array(ev[i].verified, dtype='bool')
                continue

            if isinstance(strain_idx, int):
                strain_idx = [strain_idx]

            ### laod relevant count data from HDF5
            segments = sp.atleast_2d(IN['segments'][gr_idx_segs, :])[:, strain_idx]
            seg_pos = sp.atleast_2d(IN['seg_pos'][gr_idx_segs, :])[:, strain_idx]
            edges = sp.atleast_2d(IN['edges'][gr_idx_edges, :])[:, strain_idx]
            curr_edge_idx = edge_idx[gr_idx_edges]

            for s_idx in range(len(strain_idx)):
                #if s_idx > 0 and s_idx % 50 == 0:
                #    sys.stdout.write('%i (%i)\n' % (s_idx, len(strain_idx)))
               # ev_tmp.subset_strain(s_idx) ### TODO 
                #sys.stdout.flush()
                if event_type == 'exon_skip':
                    ver, info = verify_exon_skip(ev[i], genes[g_idx], segments[:, s_idx].T,  sp.c_[curr_edge_idx, edges[:, s_idx]], options)
                elif event_type in ['alt_3prime', 'alt_5prime']:
                    ver, info = verify_alt_prime(ev[i], genes[g_idx], segments[:, s_idx].T,  sp.c_[curr_edge_idx, edges[:, s_idx]], options)
                elif event_type == 'intron_retention':
                    ver, info = verify_intron_retention(ev[i], genes[g_idx], segments[:, s_idx].T,  sp.c_[curr_edge_idx, edges[:, s_idx]], seg_pos[:, s_idx].T, options)
                elif event_type == 'mult_exon_skip':
                    ver, info = verify_mult_exon_skip(ev[i], genes[g_idx], segments[:, s_idx].T,  sp.c_[curr_edge_idx, edges[:, s_idx]], options)
                elif event_type == 'mutex_exons':
                    ver, info = verify_mutex_exons(ev[i], genes[g_idx], segments[:, s_idx].T,  sp.c_[curr_edge_idx, edges[:, s_idx]], options)

                ev[i].verified.append(ver)
                if s_idx == 0:
                    counts.append(sp.array([info], dtype='float'))
                else:
                    counts[-1] = sp.r_[counts[-1], sp.array([info], dtype='float')]
            ev[i].verified = sp.array(ev[i].verified, dtype='bool')

        IN.close()
        counts = sp.dstack(counts)
    ev = ev[old_idx]
    counts = counts[:, :, old_idx]

    if out_fn is not None:
        pickle.dump((ev, counts), open(out_fn, 'wb'), -1)

    return (ev, counts)

def verify_empty(event_type):
    
    if event_type == 'exon_skip':
        verified = [0, 0, 0, 0]
        info = [1, 0, 0, 0, 0, 0, 0]
    elif event_type in ['alt_3prime', 'alt_5prime']:
        verified = [0, 0]
        info = [1, 0, 0, 0, 0]
    elif event_type == 'intron_retention':
        verified = [0, 0]
        info = [1, 0, 0, 0, 0, 0]
    elif event_type == 'mult_exon_skip':
        verified = [0, 0, 0, 0, 0]
        info = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    elif event_type == 'mutex_exons':
        verified = [0, 0, 0, 0]
        info = [1, 0, 0, 0, 0, 0, 0, 0, 0]

    return (verified, info)
