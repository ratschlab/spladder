#! /usr/bin/env python
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2014 Andre Kahles, Jonas Behr, Gunnar Raetsch
# Copyright (C) 2009-2011 Max Planck Society
# Copyright (C) 2012-2014 Memorial Sloan-Kettering Cancer Center
#
# SplAdder wrapper script to start the interpreter with the correct list of arguments

import sys
import os
import scipy as sp
import pickle
import h5py

from . import settings
from .core.spladdercore import spladder_core
from .alt_splice.collect import collect_events
from .alt_splice.analyze import analyze_events
from .count import count_graph_coverage_wrapper, collect_single_quantification_results
from .editgraph import filter_by_edgecount
from . import init
from . import rproc as rp
from .merge import run_merge
from .helpers import *

from .classes.gene import Gene
sys.modules['modules.classe.gene'] = Gene

def spladder(options):

    ### parse parameters from options object
    options = settings.parse_args(options)

    ### load confidence level settings
    if not options.no_reset_conf:
        options = settings.set_confidence_level(options)

    ### do not compute components of merged set, if result file already exists
    fn_out_merge = get_filename('fn_out_merge', options)
    fn_out_merge_val = get_filename('fn_out_merge_val', options)

    if options.spladderfile == '-' and not os.path.exists(fn_out_merge):
        ### iterate over files, if merge strategy is single
        if options.merge in ['single', 'merge_graphs']:
            idxs = list(range(len(options.samples)))
        else:
            idxs = [0]
        
        ### set parallelization
        if options.pyproc:
            jobinfo = []

        ### create out-directory
        if not os.path.exists(options.outdir):
            os.makedirs(options.outdir)

        ### create spladder sub-directory
        if not os.path.exists(os.path.join(options.outdir, 'spladder')):
            os.makedirs(os.path.join(options.outdir, 'spladder'))

        ### pre-process annotation, if necessary
        if options.annotation.split('.')[-1] != 'pickle':
            if not os.path.exists(options.annotation + '.pickle'):
                if options.annotation.split('.')[-1].lower() in ['gff', 'gff3']:
                    (genes, options) = init.init_genes_gff3(options)
                elif options.annotation.split('.')[-1].lower() in ['gtf']:
                    (genes, options) = init.init_genes_gtf(options)
                else:
                    print('ERROR: Unknown annotation format. File needs to end in gtf or gff/gff3\nCurrent file: %s' % (options.annotation), file=sys.stderr)
                    sys.exit(1)
            options.annotation += '.pickle'

        ### add anotation contigs into lookup table
        if not hasattr(options, 'genes'):
            genes = pickle.load(open(options.annotation, 'rb'), encoding='latin1')
        else:
            genes = options.genes
        options = init.append_chrms(sp.unique(sp.array([x.chr for x in genes], dtype='str')), options)
        del genes

        ### convert input BAMs to sparse arrays - filtered case
        if options.sparse_bam:
            for bfn in options.bam_fnames:
                if bfn.endswith('bam') and not os.path.exists(re.sub(r'.bam$', '', bfn) + '.conf_%i' % options.confidence + '.filt.hdf5'):

                    if not hasattr(options, 'chrm_lookup'):
                        IN = pysam.Samfile(bfn, 'rb')
                        options = append_chrms([x['SN'] for x in parse_header(IN.text)['SQ']], options)
                        IN.close()

                    OUT = h5py.File(re.sub(r'.bam$', '', bfn) + '.conf_%i' % options.confidence + '.filt.hdf5', 'w')
                    if options.parallel > 1:
                        import multiprocessing as mp
                        pool = mp.Pool(processes=options.parallel)
                        result = [pool.apply_async(summarize_chr, args=(bfn, str(chrm), options,), kwds={'filter':options.read_filter}) for chrm in sorted(options.chrm_lookup)]
                        while result:
                            tmp = result.pop(0).get()
                            OUT.create_dataset(name=(tmp[0] + '_reads_row'), data=tmp[1].row.astype('uint8'), compression='gzip')
                            OUT.create_dataset(name=(tmp[0] + '_reads_col'), data=tmp[1].col, compression='gzip')
                            OUT.create_dataset(name=(tmp[0] + '_reads_dat'), data=tmp[1].data, compression='gzip')
                            OUT.create_dataset(name=(tmp[0] + '_reads_shp'), data=tmp[1].shape)
                            OUT.create_dataset(name=(tmp[0] + '_introns_m'), data=tmp[2], compression='gzip')
                            OUT.create_dataset(name=(tmp[0] + '_introns_p'), data=tmp[3], compression='gzip')
                            del tmp
                    else:
                        for chrm in options.chrm_lookup:
                            tmp = summarize_chr(bfn, str(chrm), options, filter=options.read_filter)
                            OUT.create_dataset(name=(chrm + '_reads_row'), data=tmp[1].row.astype('uint8'), compression='gzip')
                            OUT.create_dataset(name=(chrm + '_reads_col'), data=tmp[1].col, compression='gzip')
                            OUT.create_dataset(name=(chrm + '_reads_dat'), data=tmp[1].data, compression='gzip')
                            OUT.create_dataset(name=(chrm + '_reads_shp'), data=tmp[1].shape)
                            OUT.create_dataset(name=(chrm + '_introns_m'), data=tmp[2], compression='gzip')
                            OUT.create_dataset(name=(chrm + '_introns_p'), data=tmp[3], compression='gzip')
                    OUT.close()
                elif options.verbose:
                    print('Filtered sparse BAM representation for %s already exists.' % bfn, file=sys.stdout)

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

            ### assemble out filename to check if we are already done
            fn_out = options.out_fname
            if options.do_prune:
                fn_out = re.sub('.pickle$', '_pruned.pickle', fn_out)
            if options.do_gen_isoforms:
                fn_out = re.sub('.pickle$', '_with_isoforms.pickle', fn_out)
    
            if os.path.exists(fn_out):
                print('%s - All result files already exist.' % fn_out, file=sys.stdout)
            else:
                if options.pyproc:
                    jobinfo.append(rp.rproc('spladder_core', options, 15000, options.options_rproc, 60*60))
                else:
                    spladder_core(options)

            for key in CFG_:
                setattr(options, key, CFG_[key].copy())

        ### collect results after parallelization
        if options.pyproc:
            rp.rproc_wait(jobinfo, 30, 1.0, -1)

        ### merge parts if necessary
        if options.merge == 'merge_graphs':
            run_merge(options)

    if options.spladderfile != '-' and options.merge == 'merge_graphs' and options.validate_sg and not os.path.exists(fn_out_merge_val):
        (genes, inserted) = pickle.load(open(fn_out_merge, 'rb'))
        genes = filter_by_edgecount(genes, options)
        pickle.dump((genes, inserted), open(fn_out_merge_val, 'wb'), -1)
        del genes

    ### convert input BAMs to sparse arrays - unfiltered case
    if options.sparse_bam:
        for bfn in options.bam_fnames:
            if bfn.endswith('bam') and not os.path.exists(re.sub(r'.bam$', '', bfn) + '.hdf5'):
                #cnts = dict()

                if not hasattr(options, 'chrm_lookup'):
                    IN = pysam.Samfile(bfn, 'rb')
                    options = append_chrms([x['SN'] for x in parse_header(IN.text)['SQ']], options)
                    IN.close()

                OUT = h5py.File(re.sub(r'.bam$', '', bfn) + '.hdf5', 'w') 
                if options.parallel > 1:
                    import multiprocessing as mp
                    pool = mp.Pool(processes=options.parallel)
                    result = [pool.apply_async(summarize_chr, args=(bfn, str(chrm), options,)) for chrm in sorted(options.chrm_lookup)]
                    while result:
                        tmp = result.pop(0).get()
                        OUT.create_dataset(name=(tmp[0] + '_reads_row'), data=tmp[1].row.astype('uint8'), compression='gzip')
                        OUT.create_dataset(name=(tmp[0] + '_reads_col'), data=tmp[1].col, compression='gzip')
                        OUT.create_dataset(name=(tmp[0] + '_reads_dat'), data=tmp[1].data, compression='gzip')
                        OUT.create_dataset(name=(tmp[0] + '_reads_shp'), data=tmp[1].shape)
                        OUT.create_dataset(name=(tmp[0] + '_introns_m'), data=tmp[2], compression='gzip')
                        OUT.create_dataset(name=(tmp[0] + '_introns_p'), data=tmp[3], compression='gzip')
                else:
                    for chrm in options.chrm_lookup:
                        tmp = summarize_chr(bfn, str(chrm), options)
                        OUT.create_dataset(name=(chrm + '_reads_row'), data=tmp[1].row.astype('uint8'), compression='gzip')
                        OUT.create_dataset(name=(chrm + '_reads_col'), data=tmp[1].col, compression='gzip')
                        OUT.create_dataset(name=(chrm + '_reads_dat'), data=tmp[1].data, compression='gzip')
                        OUT.create_dataset(name=(chrm + '_reads_shp'), data=tmp[1].shape)
                        OUT.create_dataset(name=(chrm + '_introns_m'), data=tmp[2], compression='gzip')
                        OUT.create_dataset(name=(chrm + '_introns_p'), data=tmp[3], compression='gzip')
                OUT.close()
            elif options.verbose:
                print('Sparse BAM representation for %s already exists.' % bfn, file=sys.stdout)

    if options.merge == 'single' or options.qmode == 'collect':
        idxs = list(range(len(options.samples)))
    else:
        idxs = [0]

    for idx in idxs:
        ### get count output file
        if options.merge == 'single':
            fn_in_count = get_filename('fn_count_in', options, sample_idx=idx)
            fn_out_count = get_filename('fn_count_out', options, sample_idx=idx)
        elif options.merge == 'merge_graphs' and options.qmode == 'single':
            fn_in_count = get_filename('fn_count_in', options)
            fn_out_count = get_filename('fn_count_out', options, sample_idx=0)
        else:
            if options.qmode != 'collect':
                fn_in_count = get_filename('fn_count_in', options)
            fn_out_count = get_filename('fn_count_out', options)

        ### count segment graph
        if options.extract_as or options.quantify_graph:
            if not os.path.exists(fn_out_count):
                if options.merge == 'single':
                    count_graph_coverage_wrapper(fn_in_count, fn_out_count, options, sample_idx=idx)
                elif options.merge == 'merge_graphs' and options.qmode == 'single':
                    count_graph_coverage_wrapper(fn_in_count, fn_out_count, options, qmode='single')
                elif options.merge == 'merge_graphs' and options.qmode == 'collect':
                    collect_single_quantification_results(fn_out_count, idxs, options)
                else:
                    count_graph_coverage_wrapper(fn_in_count, fn_out_count, options)

        ### count intron coverage phenotype
        if options.intron_cov:
            fn_out_intron_count = fn_out_count.replace('pickle', 'introns.pickle')
            count_intron_coverage_wrapper(fn_in_count, fn_out_intron_count, options)

    ### handle alternative splicing part
    if options.extract_as:
        collect_events(options)

        for idx in idxs:
            for e_idx in range(len(options.event_types)):
                if options.merge == 'single':
                    analyze_events(options, options.event_types[e_idx], sample_idx=idx)
                else:
                    analyze_events(options, options.event_types[e_idx])


