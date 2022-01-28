import numpy as np
import pickle
import h5py
import sys
import os

import multiprocessing as mp 
import signal as sig

if __package__ is None:
    __package__ = 'modules.alt_splice'

from ..utils import *
from ..helpers import log_progress
from ..classes.segmentgraph import Segmentgraph

def verify_mult_exon_skip(event, gene, counts_segments, counts_edges, options):

    verified = [0, 0, 0, 0, 0]
    # (0) exon coordinates are valid (>= 0 && start < stop && non-overlapping) & skipped exon coverage >= FACTOR * mean(pre, after)
    # (1) inclusion count first intron >= threshold
    # (2) inclusion count last intron >= threshold
    # (3) avg inclusion count inner exons >= threshold
    # (4) skip count >= threshold

    info = np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float')
    # (0) valid, (1) exon_pre_cov (e1_cov), (2) exons_cov (e2_cov), (3) exon_aft_cov (e3_cov),
    # (4) exon_pre_exon_conf (e1e2_conf), (5) exon_exon_aft_conf (e2e3_conf), (6) exon_pre_exon_aft_conf (e1e3_conf),
    # (7) sum_inner_exon_conf (sum_e2_conf), (8) num_inner_exon (num_e2), (9) len_inner_exon (len_e2)

    ### check validity of exon coordinates (>=0)
    if np.any(event.exons1 < 0) or np.any(event.exons2 < 0):
        info[0] = 0
        return (verified, info)
    ### check validity of exon coordinates (start < stop && non-overlapping)
    elif np.any(event.exons1[:, 1] - event.exons1[:, 0] < 1) or np.any(event.exons2[:, 1] - event.exons2[:, 0] < 1):
        info[0] = 0
        return (verified, info)

    sg = gene.splicegraph
    segs = gene.segmentgraph

    ### find exons corresponding to event
    idx_exon_pre  = np.where((sg.vertices[0, :] == event.exons2[0, 0]) & (sg.vertices[1, :] == event.exons2[0, 1]))[0]
    idx_exon_aft  = np.where((sg.vertices[0, :] == event.exons2[-1, 0]) & (sg.vertices[1, :] == event.exons2[-1, 1]))[0]
    seg_exons = []
    for i in range(1, event.exons2.shape[0] - 1):
        tmp = np.where((sg.vertices[0, :] == event.exons2[i, 0]) & (sg.vertices[1, :] == event.exons2[i, 1]))[0]
        seg_exons.append(np.where(segs.seg_match[tmp, :])[1])
    
    ### find introns corresponding to event
    intron_pre = (event.exons2[0, 1], event.exons2[1, 0])
    intron_aft = (event.exons2[-2, 1], event.exons2[-1, 0])
    intron_skip = (event.exons2[0, 1], event.exons2[-1, 0])
    introns_inner = set()
    for i in range(1, event.exons2.shape[0] - 2):
        introns_inner.add((event.exons2[i, 1], event.exons2[i + 1, 0]))

    ### find segments corresponding to exons
    seg_exon_pre = np.sort(np.where(segs.seg_match[idx_exon_pre, :])[1])
    seg_exon_aft = np.sort(np.where(segs.seg_match[idx_exon_aft, :])[1])
    seg_exons_u = np.sort(np.unique([x for sublist in seg_exons for x in sublist]))

    seg_lens = segs.segments[1, :] - segs.segments[0, :]

    # exon_pre_cov
    info[1] = np.sum(counts_segments[seg_exon_pre] * seg_lens[seg_exon_pre]) / np.sum(seg_lens[seg_exon_pre])
    # exon_aft_cov
    info[3] = np.sum(counts_segments[seg_exon_aft] * seg_lens[seg_exon_aft]) / np.sum(seg_lens[seg_exon_aft])
    # exons_cov
    info[2] = np.sum(counts_segments[seg_exons_u] * seg_lens[seg_exons_u]) / np.sum(seg_lens[seg_exons_u])

    ### check if coverage of skipped exon is >= than FACTOR times average of pre and after
    if (info[2] >= options.mult_exon_skip['min_skip_rel_cov'] * (info[1] + info[3]) / 2) or \
       (options.use_anno_support and intron_pre in gene.introns_anno and intron_aft in gene.introns_anno): 
        verified[0] = 1

    ### check intron confirmation as sum of valid intron scores
    ### intron score is the number of reads confirming this intron
    # exon_pre_exon_conf
    idx = np.where(counts_edges[:, 0] == np.ravel_multi_index([seg_exon_pre[-1], seg_exons[0][0]], segs.seg_edges.shape))[0]
    if len(idx.shape) > 0 and idx.shape[0] > 0:
        info[4] = counts_edges[idx[0], 1]
    # exon_exon_aft_conf
    idx = np.where(counts_edges[:, 0] == np.ravel_multi_index([seg_exons[-1][-1], seg_exon_aft[0]], segs.seg_edges.shape))[0]
    if len(idx.shape) > 0 and idx.shape[0] > 0:
        info[5] = counts_edges[idx[0], 1]
    # exon_pre_exon_aft_conf
    idx = np.where(counts_edges[:, 0] == np.ravel_multi_index([seg_exon_pre[-1], seg_exon_aft[0]], segs.seg_edges.shape))[0]
    if len(idx.shape) > 0 and idx.shape[0] > 0:
        info[6] = counts_edges[idx[0], 1]
    for i in range(len(seg_exons) - 1):
        # sum_inner_exon_conf
        idx = np.where(counts_edges[:, 0] == np.ravel_multi_index([seg_exons[i][-1], seg_exons[i+1][0]], segs.seg_edges.shape))[0]
        if len(idx.shape) > 0 and idx.shape[0] > 0:
            info[7] += counts_edges[idx[0], 1]

    # num_inner_exon
    info[8] = event.exons2.shape[0] - 2
    info[9] = np.sum(event.exons2[1:-1, 1] - event.exons2[1:-1, 0])
    # exon_pre_exon_conf
    if (info[4] >= options.mult_exon_skip['min_non_skip_count']) or \
       (options.use_anno_support and intron_pre in gene.introns_anno):
        verified[1] = 1
    # exon_exon_aft_conf
    if (info[5] >= options.mult_exon_skip['min_non_skip_count']) or \
       (options.use_anno_support and intron_aft in gene.introns_anno):
        verified[2] = 1
    # sum_inner_exon_conf
    if ((info[7] / info[8]) >= options.mult_exon_skip['min_non_skip_count']) or \
       (options.use_anno_support and len(introns_inner.intersection(gene.introns_anno)) == len(introns_inner)):
        verified[3] = 1 
    # exon_pre_exon_aft_conf
    if (info[6] >= options.mult_exon_skip['min_skip_count']) or \
       (options.use_anno_support and intron_skip in gene.introns_anno):
        verified[4] = 1 

    return (verified, info)


def verify_intron_retention(event, gene, counts_segments, counts_edges, counts_seg_pos, options):

    verified = [0, 0]
    # (0) counts meet criteria for min_retention_cov, min_retention_region and min_retetion_rel_cov 
    # (1) min_non_retention_count >= threshold

    info = np.array([1, 0, 0, 0, 0, 0], dtype='float')
    # (0) valid, (1) exon1_cov (e1_cov), (2) intron_cov (e2_cov), (3) exon2_cov (e3_cov),
    # (4) intron_conf (e1e3_conf), (5) intron_cov_region (e2_cov_region)

    ### check validity of exon coordinates (>=0)
    if np.any(event.exons1 < 0) or np.any(event.exons2 < 0):
        info[0] = 0
        return (verified, info)
    ### check validity of exon coordinates (start < stop && non-overlapping)
    elif np.any(event.exons1[:, 1] - event.exons1[:, 0] < 1) or np.any(event.exons2[:, 1] - event.exons2[:, 0] < 1):
        info[0] = 0
        return (verified, info)

    sg = gene.splicegraph
    segs = gene.segmentgraph

    ### find exons corresponding to event
    idx_exon1  = np.where((sg.vertices[0, :] == event.exons1[0, 0]) & (sg.vertices[1, :] == event.exons1[0, 1]))[0]
    idx_exon2  = np.where((sg.vertices[0, :] == event.exons1[1, 0]) & (sg.vertices[1, :] == event.exons1[1, 1]))[0]

    ### find intron and exon corresponding to event
    intron = (event.exons1[0, 1], event.exons1[1, 0])
    exon_long =  (event.exons1[0, 0], event.exons1[1, 1])

    ### find segments corresponding to exons
    seg_exon1 = np.sort(np.where(segs.seg_match[idx_exon1, :])[1])
    seg_exon2 = np.sort(np.where(segs.seg_match[idx_exon2, :])[1])
    seg_all = np.arange(seg_exon1[0], seg_exon2[-1])

    seg_intron = np.setdiff1d(seg_all, seg_exon1)
    seg_intron = np.setdiff1d(seg_intron, seg_exon2)
    assert(seg_intron.shape[0] > 0)

    seg_lens = segs.segments[1, :] - segs.segments[0, :]

    ### compute exon coverages as mean of position wise coverage
    # exon1_cov
    info[1] = np.sum(counts_segments[seg_exon1] * seg_lens[seg_exon1]) / np.sum(seg_lens[seg_exon1])
    # exon2_cov
    info[3] = np.sum(counts_segments[seg_exon2] * seg_lens[seg_exon2]) / np.sum(seg_lens[seg_exon2])
    # intron_cov
    info[2] = np.sum(counts_segments[seg_intron] * seg_lens[seg_intron]) / np.sum(seg_lens[seg_intron])
    # intron_cov_region
    info[5] = np.sum(counts_seg_pos[seg_intron]) / np.sum(seg_lens[seg_intron])

    ### check if counts match verification criteria
    if (info[2] > options.intron_retention['min_retention_cov'] and \
        info[5] > options.intron_retention['min_retention_region'] and \
        info[2] >= options.intron_retention['min_retention_rel_cov'] * (info[1] + info[3]) / 2) or \
       (options.use_anno_support and exon_long in set([(_[0], _[1]) for _ in np.vstack(gene.exons)])):
        verified[0] = 1

    ### check intron confirmation as sum of valid intron scores
    ### intron score is the number of reads confirming this intron
    # intron conf
    idx = np.where(counts_edges[:, 0] == np.ravel_multi_index([seg_exon1[-1], seg_exon2[0]], segs.seg_edges.shape))[0]
    info[4] = counts_edges[idx, 1]

    if (info[4] >= options.intron_retention['min_non_retention_count']) or \
       (options.use_anno_support and intron in gene.introns_anno):
        verified[1] = 1

    return (verified, info)


def verify_exon_skip(event, gene, counts_segments, counts_edges, options):

    verified = [0, 0, 0, 0]
    # (0) coverage of skipped exon is >= than FACTOR * mean(pre, after)
    # (1) inclusion count of first intron >= threshold 
    # (2) inclusion count of second intron >= threshold
    # (3) skip count of exon >= threshold

    info = np.array([1, 0, 0, 0, 0, 0, 0], dtype='float')
    # (0) valid, (1) exon_pre_cov (e1_cov), (2) exon_cov (e2_cov), (3) exon_aft_cov (e3_cov), 
    # (4) exon_pre_exon_conf (e1e2_conf), (5) exon_exon_aft_conf (e2e3_conf), (6) exon_pre_exon_aft_conf (e1e3_conf)

    ### check validity of exon coordinates (>=0)
    if np.any(event.exons1 < 0) or np.any(event.exons2 < 0):
        info[0] = False
        return (verified, info)
    ### check validity of exon coordinates (start < stop && non-overlapping)
    elif np.any(event.exons1[:, 1] - event.exons1[:, 0] < 1) or np.any(event.exons2[:, 1] - event.exons2[:, 0] < 1):
        info[0] = False
        return (verified, info)

    sg = gene.splicegraph
    segs = gene.segmentgraph

    ### find exons corresponding to event
    idx_exon_pre = np.where((sg.vertices[0, :] == event.exons2[0, 0]) & (sg.vertices[1, :] == event.exons2[0, 1]))[0]
    idx_exon = np.where((sg.vertices[0, :] == event.exons2[1, 0]) & (sg.vertices[1, :] == event.exons2[1, 1]))[0]
    idx_exon_aft = np.where((sg.vertices[0, :] == event.exons2[2, 0]) & (sg.vertices[1, :] == event.exons2[2, 1]))[0]

    ### find introns corresponding to event
    intron_pre = (event.exons2[0, 1], event.exons2[1, 0])
    intron_aft = (event.exons2[1, 1], event.exons2[2, 0])
    intron_skip = (event.exons2[0, 1], event.exons2[2, 0])

    ### find segments corresponding to exons
    seg_exon_pre = np.sort(np.where(segs.seg_match[idx_exon_pre, :])[1])
    seg_exon_aft = np.sort(np.where(segs.seg_match[idx_exon_aft, :])[1])
    seg_exon = np.sort(np.where(segs.seg_match[idx_exon, :])[1])

    seg_lens = segs.segments[1, :] - segs.segments[0, :]

    # exon pre cov
    info[1] = np.sum(counts_segments[seg_exon_pre] * seg_lens[seg_exon_pre]) /np.sum(seg_lens[seg_exon_pre])
    # exon aft cov
    info[3] = np.sum(counts_segments[seg_exon_aft] * seg_lens[seg_exon_aft]) /np.sum(seg_lens[seg_exon_aft])
    # exon cov
    info[2] = np.sum(counts_segments[seg_exon] * seg_lens[seg_exon]) /np.sum(seg_lens[seg_exon])

    ### check if coverage of skipped exon is >= than FACTOR times average of pre and after
    if (info[2] >= options.exon_skip['min_skip_rel_cov'] * (info[1] + info[3]) / 2) or \
       (options.use_anno_support and intron_pre in gene.introns_anno and intron_aft in gene.introns_anno): 
        verified[0] = 1

    ### check intron confirmation as sum of valid intron scores
    ### intron score is the number of reads confirming this intron
    # exon_pre_exon_conf
    idx = np.where(counts_edges[:, 0] == np.ravel_multi_index([seg_exon_pre[-1], seg_exon[0]], segs.seg_edges.shape))[0]
    info[4] = counts_edges[idx, 1]
    if (info[4] >= options.exon_skip['min_non_skip_count']) or \
       (options.use_anno_support and intron_pre in gene.introns_anno):
        verified[1] = 1
    # exon_exon_aft_conf
    idx = np.where(counts_edges[:, 0] == np.ravel_multi_index([seg_exon[-1], seg_exon_aft[0]], segs.seg_edges.shape))[0]
    info[5] = counts_edges[idx, 1]
    if (info[5] >= options.exon_skip['min_non_skip_count']) or \
       (options.use_anno_support and intron_aft in gene.introns_anno):
        verified[2] = 1
    # exon_pre_exon_aft_conf
    idx = np.where(counts_edges[:, 0] == np.ravel_multi_index([seg_exon_pre[-1], seg_exon_aft[0]], segs.seg_edges.shape))[0]
    info[6] = counts_edges[idx, 1]
    if (info[6] >= options.exon_skip['min_skip_count']) or \
       (options.use_anno_support and intron_skip in gene.introns_anno):
        verified[3] = 1

    return (verified, info)


def verify_alt_prime(event, gene, counts_segments, counts_edges, options):
    # [verified, info] = verify_exon_skip(event, fn_bam, cfg)

    verified = [0, 0]
    # (0) coverage of diff region is at least FACTOR * coverage constant region
    # (1) both alternative introns are >= threshold 

    info = np.array([1, 0, 0, 0, 0, 0], dtype='float')
    # (0) valid, (1) exon1_const_cov (e1_cov), (2) exon_diff_cov (e2_cov), (3) exon2_const_cov (e3_cov),
    # (4) e1e3_conf, (5) e2_conf

    ### check validity of exon coordinates (>=0)
    if np.any(event.exons1 < 0) or np.any(event.exons2 < 0):
        info[0] = 0 
        return (verified, info)

    ### check validity of intron coordinates (only one side is differing)
    if (event.exons1[0, 1] != event.exons2[0, 1]) and (event.exons1[1, 0] != event.exons2[1, 0]):
        info[0] = 0 
        return (verified, info)

    sg = gene.splicegraph
    segs = gene.segmentgraph

    ### find exons corresponding to event
    idx_exon11 = np.where((sg.vertices[0, :] == event.exons1[0, 0]) & (sg.vertices[1, :] == event.exons1[0, 1]))[0]
    if idx_exon11.shape[0] == 0:
        segs_exon11 = np.where((segs.segments[0, :] >= event.exons1[0, 0]) & (segs.segments[1, :] <= event.exons1[0, 1]))[0]
    else:
        segs_exon11 = np.where(segs.seg_match[idx_exon11, :])[1]
    idx_exon12 = np.where((sg.vertices[0, :] == event.exons1[1, 0]) & (sg.vertices[1, :] == event.exons1[1, 1]))[0]
    if idx_exon12.shape[0] == 0:
        segs_exon12 = np.where((segs.segments[0, :] >= event.exons1[1, 0]) & (segs.segments[1, :] <= event.exons1[1, 1]))[0]
    else:
        segs_exon12 = np.where(segs.seg_match[idx_exon12, :])[1]
    idx_exon21 = np.where((sg.vertices[0, :] == event.exons2[0, 0]) & (sg.vertices[1, :] == event.exons2[0, 1]))[0]
    if idx_exon21.shape[0] == 0:
        segs_exon21 = np.where((segs.segments[0, :] >= event.exons2[0, 0]) & (segs.segments[1, :] <= event.exons2[0, 1]))[0]
    else:
        segs_exon21 = np.where(segs.seg_match[idx_exon21, :])[1]
    idx_exon22 = np.where((sg.vertices[0, :] == event.exons2[1, 0]) & (sg.vertices[1, :] == event.exons2[1, 1]))[0]
    if idx_exon22.shape[0] == 0:
        segs_exon22 = np.where((segs.segments[0, :] >= event.exons2[1, 0]) & (segs.segments[1, :] <= event.exons2[1, 1]))[0]
    else:
        segs_exon22 = np.where(segs.seg_match[idx_exon22, :] > 0)[1]

    assert(segs_exon11.shape[0] > 0)
    assert(segs_exon12.shape[0] > 0)
    assert(segs_exon21.shape[0] > 0)
    assert(segs_exon22.shape[0] > 0)

    if segs_exon11.shape == segs_exon21.shape and np.all(segs_exon11 == segs_exon21):
        seg_const1 = segs_exon11
        seg_diff = np.setdiff1d(segs_exon12, segs_exon22)
        if seg_diff.shape[0] == 0:
            seg_diff = np.setdiff1d(segs_exon22, segs_exon12)
        seg_const2 = np.intersect1d(segs_exon12, segs_exon22)
    elif segs_exon12.shape == segs_exon22.shape and np.all(segs_exon12 == segs_exon22):
        seg_const2 = segs_exon12
        seg_diff = np.setdiff1d(segs_exon11, segs_exon21)
        if seg_diff.shape[0] == 0:
            seg_diff = np.setdiff1d(segs_exon21, segs_exon11)
        seg_const1 = np.intersect1d(segs_exon21, segs_exon11)
    else:
        print("ERROR: both exons differ in alt prime event in verify_alt_prime", file=sys.stderr)
        sys.exit(1)

    ### find introns corresponding to event
    intron1 = (event.exons1[0, 1], event.exons1[1, 0])
    intron2 = (event.exons2[0, 1], event.exons2[1, 0])

    seg_lens = segs.segments[1, :] - segs.segments[0, :]

    # exon_diff_cov
    info[2] = np.sum(counts_segments[seg_diff] * seg_lens[seg_diff]) / np.sum(seg_lens[seg_diff])
    # exon1_const_cov
    info[1] = np.sum(counts_segments[seg_const1] * seg_lens[seg_const1]) / np.sum(seg_lens[seg_const1])
    # exon2_const_cov
    info[3] = np.sum(counts_segments[seg_const2] * seg_lens[seg_const2]) / np.sum(seg_lens[seg_const2])

    if (info[2] >= options.alt_prime['min_diff_rel_cov'] * ((info[1] * np.sum(seg_lens[seg_const1])) + (info[3] * np.sum(seg_lens[seg_const2]))) / np.sum(seg_lens[np.r_[seg_const1, seg_const2]])) or \
       (options.use_anno_support and intron1 in gene.introns_anno and intron2 in gene.introns_anno):
        verified[0] = 1

    ### check intron confirmations as sum of valid intron scores
    ### intron score is the number of reads confirming this intron
    # intron1_conf 
    idx = np.where(counts_edges[:, 0] == np.ravel_multi_index([segs_exon11[-1], segs_exon12[0]], segs.seg_edges.shape))[0]
    assert(idx.shape[0] > 0)
    info[4] = counts_edges[idx, 1]
    # intron2_conf 
    idx = np.where(counts_edges[:, 0] == np.ravel_multi_index([segs_exon21[-1], segs_exon22[0]], segs.seg_edges.shape))[0]
    assert(idx.shape[0] > 0)
    info[5] = counts_edges[idx, 1]

    if (min(info[4], info[5]) >= options.alt_prime['min_intron_count']) or \
       (options.use_anno_support and intron1 in gene.introns_anno and intron2 in gene.introns_anno):
        verified[1] = 1

    return (verified, info)

def verify_mutex_exons(event, gene, counts_segments, counts_edges, options):
    #

    verified = [0, 0, 0, 0]
    # (0) coverage of first alt exon is >= than FACTOR times average of pre and after 
    # (1) coverage of second alt exon is >= than FACTOR times average of pre and after 
    # (2) both introns neighboring first alt exon are confirmed >= threshold
    # (3) both introns neighboring second alt exon are confirmed >= threshold

    info = np.array([1, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float')
    # (0) valid, (1) exon_pre_cov (e1_cov), (2) exon1_cov (e2_cov), (3) exon2_cov (e3_cov), (4) exon_aft_cov (e4_cov), 
    # (5) exon_pre_exon1_conf (e1e2_conf), (6) exon_pre_exon2_conf (e1e3_conf), (7) exon1_exon_aft_conf (e2e4_conf), (8) exon2_exon_aft_conf (e3e4_conf)

    ### check validity of exon coordinates (>=0)
    if np.any(event.exons1 < 0) or np.any(event.exons2 < 0):
        info[0] = 0
        return (verified, info)
    ### check validity of exon coordinates (start < stop && non-overlapping)
    elif np.any(event.exons1[:, 1] - event.exons1[:, 0] < 1) or np.any(event.exons2[:, 1] - event.exons2[:, 0] < 1) or \
         (event.exons1[1, 1] > event.exons2[1, 0] and event.exons1[1, 0] < event.exons2[1, 0]) or \
         (event.exons2[1, 1] > event.exons1[1, 0] and event.exons2[1, 0] < event.exons1[1, 0]):
        info[0] = 0
        return (verified, info)

    sg = gene.splicegraph
    segs = gene.segmentgraph

    ### find exons corresponding to event
    idx_exon_pre  = np.where((sg.vertices[0, :] == event.exons1[0, 0]) & (sg.vertices[1, :] == event.exons1[0, 1]))[0]
    idx_exon_aft  = np.where((sg.vertices[0, :] == event.exons1[-1, 0]) & (sg.vertices[1, :] == event.exons1[-1, 1]))[0]
    idx_exon1  = np.where((sg.vertices[0, :] == event.exons1[1, 0]) & (sg.vertices[1, :] == event.exons1[1, 1]))[0]
    idx_exon2  = np.where((sg.vertices[0, :] == event.exons2[1, 0]) & (sg.vertices[1, :] == event.exons2[1, 1]))[0]
    
    ### find introns and exons corresponding to event
    intron_pre_ex1 = (event.exons1[0, 1], event.exons1[1, 0])
    intron_pre_ex2 = (event.exons2[0, 1], event.exons2[1, 0])
    intron_ex1_aft = (event.exons1[1, 1], event.exons1[2, 0])
    intron_ex2_aft = (event.exons2[1, 1], event.exons2[2, 0])
    exon1 = (event.exons1[1, 0], event.exons1[1, 1])
    exon2 = (event.exons2[1, 0], event.exons2[1, 1])

    ### find segments corresponding to exons
    seg_exon_pre = np.sort(np.where(segs.seg_match[idx_exon_pre, :])[1])
    seg_exon_aft = np.sort(np.where(segs.seg_match[idx_exon_aft, :])[1])
    seg_exon1 = np.sort(np.where(segs.seg_match[idx_exon1, :])[1])
    seg_exon2 = np.sort(np.where(segs.seg_match[idx_exon2, :])[1])

    seg_lens = segs.segments[1, :] - segs.segments[0, :]

    # exon pre cov
    info[1] = np.sum(counts_segments[seg_exon_pre] * seg_lens[seg_exon_pre]) / np.sum(seg_lens[seg_exon_pre])
    # exon1 cov
    info[2] = np.sum(counts_segments[seg_exon1] * seg_lens[seg_exon1]) / np.sum(seg_lens[seg_exon1])
    # exon2 cov
    info[3] = np.sum(counts_segments[seg_exon2] * seg_lens[seg_exon2]) / np.sum(seg_lens[seg_exon2])
    # exon aft cov
    info[4] = np.sum(counts_segments[seg_exon_aft] * seg_lens[seg_exon_aft]) / np.sum(seg_lens[seg_exon_aft])

    ### check if coverage of first exon is >= than FACTOR times average of pre and after
    if (info[2] >= options.mutex_exons['min_skip_rel_cov'] * (info[1] + info[4])/2) or \
       (options.use_anno_support and exon1 in set([(_[0], _[1]) for _ in np.hstack(gene.exons)])):
        verified[0] = 1
    if (info[3] >= options.mutex_exons['min_skip_rel_cov'] * (info[1] + info[4])/2) or \
       (options.use_anno_support and exon2 in set([(_[0], _[1]) for _ in np.hstack(gene.exons)])):
        verified[1] = 1

    ### check intron confirmation as sum of valid intron scores
    ### intron score is the number of reads confirming this intron
    # exon_pre_exon1_conf
    idx = np.where(counts_edges[:, 0] == np.ravel_multi_index([seg_exon_pre[-1], seg_exon1[0]], segs.seg_edges.shape))[0]
    if len(idx.shape) > 0 and idx.shape[0] > 0:
        info[5] = counts_edges[idx[0], 1]
    # exon_pre_exon2_conf
    idx = np.where(counts_edges[:, 0] == np.ravel_multi_index([seg_exon_pre[-1], seg_exon2[0]], segs.seg_edges.shape))[0]
    if len(idx.shape) > 0 and idx.shape[0] > 0:
        info[6] = counts_edges[idx[0], 1]
    # exon1_exon_aft_conf
    idx = np.where(counts_edges[:, 0] == np.ravel_multi_index([seg_exon1[-1], seg_exon_aft[0]], segs.seg_edges.shape))[0]
    if len(idx.shape) > 0 and idx.shape[0] > 0:
        info[7] = counts_edges[idx[0], 1]
    # exon2_exon_aft_conf
    idx = np.where(counts_edges[:, 0] == np.ravel_multi_index([seg_exon2[-1], seg_exon_aft[0]], segs.seg_edges.shape))[0]
    if len(idx.shape) > 0 and idx.shape[0] > 0:
        info[8] = counts_edges[idx[0], 1]

    # set verification flags for intron confirmation
    if (min(info[5], info[7]) >= options.mutex_exons['min_conf_count']) or \
       (options.use_anno_support and intron_pre_ex1 in gene.introns_anno and intron_ex1_aft in gene.introns_anno):
        verified[2] = 1
    if (min(info[6], info[8]) >= options.mutex_exons['min_conf_count']) or \
       (options.use_anno_support and intron_pre_ex2 in gene.introns_anno and intron_ex2_aft in gene.introns_anno):
        verified[3] = 1

    return (verified, info)

def verify_wrapper(ev, genes, gidx_min, gene_ids_edges, gene_ids_segs, edge_idx, sample_idx, event_type, fn_count, options, idx):

    IN = h5py.File(fn_count, 'r')
    counts = []
    verified = []
    for i in range(ev.shape[0]):
        
        #sys.stdout.write('.')
        #if i > 0 and i % 50 == 0:
        #    sys.stdout.write('%i (%i)\n' % (i, ev.shape[0]))
        #sys.stdout.flush()
    
        g_idx = ev[i].gene_idx
        _g_idx = g_idx - gidx_min

        if genes[_g_idx].segmentgraph is None or genes[_g_idx].segmentgraph.segments.shape[1] == 0:
            genes[_g_idx].from_sparse()
            genes[_g_idx].segmentgraph = Segmentgraph(genes[_g_idx])
            genes[_g_idx].to_sparse()

        ### there are no edges present in the event
        if gene_ids_edges.shape[0] == 0:
            ver, info = verify_empty(event_type)
            counts.append(np.array([info]))
            verified.append([ver for _ in sample_idx])
            continue

        gr_idx_segs = np.where(gene_ids_segs == g_idx)[0]
        gr_idx_edges = np.where(gene_ids_edges == g_idx)[0]
        if gr_idx_edges.shape[0] == 0:
            ver, info = verify_empty(event_type)
            counts.append(np.array([info]))
            verified.append([ver for _ in sample_idx])
            continue

        if isinstance(sample_idx, int):
            sample_idx = [sample_idx]

        ### laod relevant count data from HDF5
        segments = np.atleast_2d(IN['segments'][gr_idx_segs, :])[:, sample_idx]
        seg_pos = np.atleast_2d(IN['seg_pos'][gr_idx_segs, :])[:, sample_idx]
        edges = np.atleast_2d(IN['edges'][gr_idx_edges, :])[:, sample_idx]
        curr_edge_idx = edge_idx[gr_idx_edges]

        verified.append([])
        for s_idx in range(len(sample_idx)):
            if event_type == 'exon_skip':
                ver, info = verify_exon_skip(ev[i], genes[_g_idx], segments[:, s_idx].T,  np.c_[curr_edge_idx, edges[:, s_idx]], options)
            elif event_type in ['alt_3prime', 'alt_5prime']:
                ver, info = verify_alt_prime(ev[i], genes[_g_idx], segments[:, s_idx].T,  np.c_[curr_edge_idx, edges[:, s_idx]], options)
            elif event_type == 'intron_retention':
                ver, info = verify_intron_retention(ev[i], genes[_g_idx], segments[:, s_idx].T,  np.c_[curr_edge_idx, edges[:, s_idx]], seg_pos[:, s_idx].T, options)
            elif event_type == 'mult_exon_skip':
                ver, info = verify_mult_exon_skip(ev[i], genes[_g_idx], segments[:, s_idx].T,  np.c_[curr_edge_idx, edges[:, s_idx]], options)
            elif event_type == 'mutex_exons':
                ver, info = verify_mutex_exons(ev[i], genes[_g_idx], segments[:, s_idx].T,  np.c_[curr_edge_idx, edges[:, s_idx]], options)

            verified[-1].append(ver)
            if s_idx == 0:
                counts.append(info[np.newaxis, :])
            else:
                counts[-1] = np.append(counts[-1], info[np.newaxis, :], axis=0) 
        verified[-1] = np.array(verified[-1], dtype='bool')

    IN.close()

    return (np.dstack(counts), np.dstack(verified), idx)


def verify_all_events(ev, sample_idx=None, list_bam=None, event_type=None, options=None, out_fn=None):

    ### verify the events if demanded
    if options.verify_alt_events:

        validate_tag = ''
        if options.validate_sg:
            validate_tag = '.validated'

        if options.merge == 'single':
            (genes, inserted) = pickle.load(open('%s/spladder/genes_graph_conf%i.%s%s.pickle' % (options.outdir, options.confidence, options.samples[sample_idx], validate_tag), 'rb'))
            fn_count = '%s/spladder/genes_graph_conf%i.%s%s.count.hdf5' % (options.outdir, options.confidence, options.samples[sample_idx], validate_tag)
        else:
            (genes, inserted) = pickle.load(open('%s/spladder/genes_graph_conf%i.%s%s.pickle' % (options.outdir, options.confidence, options.merge, validate_tag), 'rb'))
            fn_count = '%s/spladder/genes_graph_conf%i.%s%s.count.hdf5' % (options.outdir, options.confidence, options.merge, validate_tag)

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
        IN.close()

        ### sort events by gene idx
        gene_idx_all = np.array([x.gene_idx for x in ev])
        s_idx = np.argsort(gene_idx_all)
        ev = ev[s_idx]
        gene_idx_all = gene_idx_all[s_idx]
        old_idx = np.argsort(s_idx)

        if options.parallel > 1:
            pool = mp.Pool(processes=options.parallel, initializer=lambda: sig.signal(sig.SIGINT, sig.SIG_IGN))
            binsize = 50
            maxsize = ev.shape[0]
            idx_chunks = [np.arange(x, min(x + binsize, maxsize)) for x in range(0, maxsize, binsize)]
            counts = np.empty((len(idx_chunks), ), dtype='object')
            verified = np.empty((len(idx_chunks), ), dtype='object')
            gidx_min = [gene_idx_all[cidx][0] for cidx in idx_chunks]   
            gidx_max = [gene_idx_all[cidx][-1] + 1 for cidx in idx_chunks]   

            try:
                result = [pool.apply_async(verify_wrapper, args=(ev[cidx], genes[gidx_min[c]:gidx_max[c]], gidx_min[c], gene_ids_edges, gene_ids_segs, edge_idx, sample_idx, event_type, fn_count, options, c)) for c,cidx in enumerate(idx_chunks)]
                res_cnt = 0
                while result:
                    tmp = result.pop(0).get()
                    counts[tmp[2]] = tmp[0]
                    verified[tmp[2]] = tmp[1]
                    if options.verbose:
                        log_progress(res_cnt, len(idx_chunks))
                        res_cnt += 1
                if options.verbose:
                    log_progress(len(idx_chunks), len(idx_chunks))
                    print('')
                pool.terminate()
                pool.join()
            except KeyboardInterrupt:
                print('Keyboard Interrupt - exiting', file=sys.stderr)
                pool.terminate()
                pool.join()
                sys.exit(1)

            ### integrate results in list into a coherent results list
            if len(counts) > 0:
                counts = np.dstack(counts)
                verified = np.dstack(verified)
        else:
            counts, verified = verify_wrapper(ev, genes, 0, gene_ids_edges, gene_ids_segs, edge_idx, sample_idx, event_type, fn_count, options, 0)[:2]

    # re-establish initial sort order
    ev = ev[old_idx]
    counts = counts[:, :, old_idx]
    verified = verified[:, :, old_idx]

    if out_fn is not None:
        pickle.dump((ev, counts, verified), open(out_fn, 'wb'), -1)

    return (ev, counts, verified)

def verify_empty(event_type):
    
    if event_type == 'exon_skip':
        verified = [0, 0, 0, 0]
        info = np.array([1, 0, 0, 0, 0, 0, 0], dtype='float')
    elif event_type in ['alt_3prime', 'alt_5prime']:
        verified = [0, 0]
        info = np.array([1, 0, 0, 0, 0, 0], dtype='float')
    elif event_type == 'intron_retention':
        verified = [0, 0]
        info = np.array([1, 0, 0, 0, 0, 0], dtype='float')
    elif event_type == 'mult_exon_skip':
        verified = [0, 0, 0, 0, 0]
        info = np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float')
    elif event_type == 'mutex_exons':
        verified = [0, 0, 0, 0]
        info = np.array([1, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float')

    return (verified, info)
