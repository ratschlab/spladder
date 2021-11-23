#! /usr/bin/env python
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2021 Andre Kahles, Jonas Behr, Gunnar Raetsch
# Copyright (C) 2009-2011 Max Planck Society
# Copyright (C) 2012-2016 Memorial Sloan-Kettering Cancer Center
# Copyright (C) 2016-2021 ETH Zurich
#
# SplAdder wrapper script to start the interpreter with the correct list of arguments

import sys
import os
import numpy as np
import pickle
import h5py

from . import settings
from .core.spladdercore import spladder_core
from .alt_splice.collect import collect_events
from .alt_splice.analyze import analyze_events
from .count import count_graph_coverage_wrapper, collect_single_quantification_results, compute_gene_expression
from .editgraph import filter_by_edgecount
from .merge import run_merge
from .helpers import *

from .spladder_prep import prep_annotation, prep_sparse_bam_filtered, prep_sparse_bam_full

from .classes.gene import Gene
sys.modules['modules.classes.gene'] = Gene

def _prep_workdir(options):

    ### create out-directory
    if not os.path.exists(options.outdir):
        os.makedirs(options.outdir)

    ### create spladder sub-directory
    if not os.path.exists(os.path.join(options.outdir, 'spladder')):
        os.makedirs(os.path.join(options.outdir, 'spladder'))

    ### create tmp sub-directory
    if not os.path.exists(options.tmpdir):
        os.makedirs(os.path.join(options.tmpdir))


def spladder(options):

    ### parse parameters from options object
    options = settings.parse_args(options)

    ### load confidence level settings
    if not options.no_reset_conf:
        options = settings.set_confidence_level(options)

    ### do not compute components of merged set, if result file already exists
    fn_out_merge = get_filename('fn_out_merge', options)
    fn_out_merge_val = get_filename('fn_out_merge_val', options)

    if not os.path.exists(fn_out_merge):
        ### iterate over files, if merge strategy is single
        if options.merge in ['single', 'merge_graphs']:
            idxs = list(range(len(options.samples)))
        else:
            idxs = [0]
        
        ### preparatory work
        _prep_workdir(options)
        options = prep_annotation(options)

        ### convert input BAMs to sparse arrays - filtered case
        if options.sparse_bam:
            prep_sparse_bam_filtered(options)

        ### build individual graphs
        for idx in idxs:
            CFG_ = dict()
            if options.merge != 'merge_bams':
                CFG_['bam_fnames'] = options.bam_fnames
                CFG_['samples'] = options.samples
                options.bam_fnames = options.bam_fnames[idx]
                options.samples = options.samples[idx]
                options.out_fname = '%s/spladder/genes_graph_conf%i.%s.pickle' % (options.outdir, options.confidence, options.samples)
            else:
                options.out_fname = '%s/spladder/genes_graph_conf%i.%s.pickle' % (options.outdir, options.confidence, options.merge)

            ### assemble output filename to check if we are already done
            fn_out = options.out_fname
            if options.do_prune:
                fn_out = re.sub('.pickle$', '_pruned.pickle', fn_out)
            if options.do_gen_isoforms:
                fn_out = re.sub('.pickle$', '_with_isoforms.pickle', fn_out)
    
            if os.path.exists(fn_out):
                print('%s - All result files already exist.' % fn_out, file=sys.stdout)
            else:
                spladder_core(options)

            for key in CFG_:
                setattr(options, key, CFG_[key].copy())

        ### merge parts if necessary
        if options.merge == 'merge_graphs':
            run_merge(options)

    ### generate validated version of splice graph
    if options.merge == 'merge_graphs' and options.validate_sg and not os.path.exists(fn_out_merge_val):
        (genes, inserted) = pickle.load(open(fn_out_merge, 'rb'))
        genes = filter_by_edgecount(genes, options)
        pickle.dump((genes, inserted), open(fn_out_merge_val, 'wb'), -1)
        del genes

    ### convert input BAMs to sparse arrays - unfiltered case
    if options.sparse_bam:
        prep_sparse_bam_full(options)
    if options.merge == 'single' or options.qmode == 'collect':
        idxs = list(range(len(options.samples)))
    else:
        idxs = [0]

    for idx in idxs:
        ### don't do any counting in chunked merge mode, unless we are in the final chunk level
        if len(options.chunked_merge) > 0:
            curr_level, max_level, chunk_start, chunk_end = [int(_) for _ in options.chunked_merge[0]]
            if curr_level < max_level:
                break

        ### get count output file
        if options.merge == 'single':
            fn_in_count = get_filename('fn_count_in', options, sample_idx=idx)
            fn_out_count = get_filename('fn_count_out', options, sample_idx=idx)
        elif options.merge == 'merge_graphs' and options.qmode == 'single':
            fn_in_count = get_filename('fn_count_in', options)
            fn_out_count = get_filename('fn_count_out', options, sample_idx=0)
        else:
            fn_in_count = get_filename('fn_count_in', options)
            fn_out_count = get_filename('fn_count_out', options)
        fn_out_gene_count = re.sub(r'.count.hdf5$', '', fn_out_count) + '.gene_exp.hdf5'

        ### count segment graph
        if options.quantify_graph:
            if not os.path.exists(fn_out_count):
                ### collect graph counts
                if options.merge == 'single':
                    count_graph_coverage_wrapper(fn_in_count, fn_out_count, options, sample_idx=idx)
                elif options.merge == 'merge_graphs' and options.qmode == 'single':
                    count_graph_coverage_wrapper(fn_in_count, fn_out_count, options, qmode='single')
                elif options.merge == 'merge_graphs' and options.qmode == 'collect':
                    collect_single_quantification_results(fn_out_count, idxs, options)
                else:
                    count_graph_coverage_wrapper(fn_in_count, fn_out_count, options)
                
            ### collect gene expression counts from graph
            if not os.path.exists(fn_out_gene_count):
                if options.merge == 'single':
                    compute_gene_expression(options, fn_in_count, fn_out_count, fn_out_gene_count, sample_idx=[0])
                else:
                    compute_gene_expression(options, fn_in_count, fn_out_count, fn_out_gene_count)

    ### handle alternative splicing part
    if options.extract_as:
        collect_events(options)

        if options.quantify_graph:
            for idx in idxs:
                for e_idx in range(len(options.event_types)):
                    if options.merge == 'single':
                        analyze_events(options, options.event_types[e_idx], sample_idx=idx)
                    else:
                        analyze_events(options, options.event_types[e_idx])


