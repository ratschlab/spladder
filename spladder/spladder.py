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

def parse_options(argv):

    """Parses options from the command line """

    def _switch(s):
        return bool(s)

    from argparse import ArgumentParser

    parser = ArgumentParser()
    required = parser.add_argument_group('MANDATORY')
    required.add_argument('-b', '--bams', dest='bams', metavar='FILE1,FILE2,...', help='alignment files in BAM format (comma separated list)', default='-')
    required.add_argument('-o', '--outdir', dest='outdir', metavar='DIR', help='output directory', default='-')
    required.add_argument('-a', '--annotation', dest='annotation', metavar='FILE', help='file name for annotation in GTF/GFF3 or format', default='-')
    general = parser.add_argument_group('GENERAL')
    general.add_argument('--parallel', dest='parallel', metavar='<INT>', type=int, help='use INT processors [1]', default=1)
    general.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='use verbose output mode [off]', default=False)
    general.add_argument('-d', '--debug', dest='debug', action='store_true', help='use debug output mode [off]', default=False)
    inputs = parser.add_argument_group('INPUT OPTIONS')
    inputs.add_argument('-n', '--readlen', dest='readlen', metavar='INT', type=int, help='read length (used for automatic confidence levele settings) [36]', default=36)
    inputs.add_argument('--primary-only', dest='primary_only', action='store_true', help='only use primary alignments [on]', default=True)
    inputs.add_argument('--no-primary-only', dest='primary_only', action='store_false', default=True)
    inputs.add_argument('--var-aware', dest='var_aware', action='store_true', help='alignment files are variation aware (presence of XM and XG tags) [off]', default=False)
    inputs.add_argument('--no-var-aware', dest='var_aware', action='store_false', default=False)
    inputs.add_argument('--labels', dest='labels', metavar='STRING', help='use labels instead of bam file names (comma separated list) [-]', default='-')
    #inputs.add_argument('-S', '--ref-strain', dest='refstrain', metavar='STRING', help='reference strain [-]', default='-')
    #inputs.add_argument('-x', '--same-genome', dest='same_genome', metavar='y|n', help='input alignments share the same genome [y]', default='y')
    inputs.add_argument('--input-graph', dest='spladderfile', metavar='FILE', help='use existing SplAdder graph as input (advanced) [-]', default='-')
    inputs.add_argument('--filter-overlap-genes', dest='filter_overlap_genes', action='store_true', help='remove genes from annotation that overlap each other [off]', default=False)
    inputs.add_argument('--filter-overlap-exons', dest='filter_overlap_exons', action='store_true', help='remove exons from annotation that overlap each other [off]', default=False)
    inputs.add_argument('--filter-overlap-transcripts', dest='filter_overlap_transcripts', action='store_true', help='remove transcripts from annotation that overlap each other [off]', default=False)
    inputs.add_argument('--ignore-mismatches', dest='ignore_mismatches', action='store_true', help='does not filter by edit operations - does not require NM in BAM [off]', default=False)
    outputs = parser.add_argument_group('OUTPUT OPTIONS')
    outputs.add_argument('-l', '--logfile', dest='logfile', metavar='FILE', help='log file name [stdout]', default='-')
    outputs.add_argument('--output-txt', dest='output_txt', action='store_true', help='outputs all events in txt format (can be big) [off]', default=False)
    outputs.add_argument('--output-txt-conf', dest='output_confirmed_txt', action='store_true', help='outputs confirmed events in txt format [on]', default=True)
    outputs.add_argument('--no-output-txt-conf', dest='output_confirmed_txt', action='store_false', default=True)
    outputs.add_argument('--output-gff3', dest='output_gff3', action='store_true', help='outputs all events in GFF3 format [off]', default=False)
    outputs.add_argument('--output-gff3-conf', dest='output_confirmed_gff3', action='store_true', help='outputs confirmed events in GFF3 format [on]', default=True)
    outputs.add_argument('--no-output-gff3-conf', dest='output_confirmed_gff3', action='store_false', default=True)
    outputs.add_argument('--output-struc', dest='output_struc', action='store_true', help='outputs all events in structured splicing syntax similar to astalavista [off]', default=False)
    outputs.add_argument('--output-struc-conf', dest='output_confirmed_struc', action='store_true', help='outputs confirmed events in structured splicing syntax similar to astalavista [off]', default=False)
    outputs.add_argument('--output-bed', dest='output_bed', action='store_true', help='output all events in BED format [off]', default=False)
    outputs.add_argument('--output-conf-bed', dest='output_confirmed_bed', action='store_true', help='output confirmed events in BED format [off]', default=False)
    outputs.add_argument('--output-conf-tcga', dest='output_confirmed_tcga', action='store_true', help='output confirmed events in format used for TCGA [off]', default=False)
    outputs.add_argument('--output-conf-icgc', dest='output_confirmed_icgc', action='store_true', help='output confirmed events in format used for ICGC [off]', default=False)
    outputs.add_argument('--sparse-bam', dest='sparse_bam', action='store_true', help='store BAM content as sparse representation for later use [off]', default=False)
    outputs.add_argument('--compress-text', dest='compress_text', action='store_true', help='compress text output [on]', default=True)
    outputs.add_argument('--no-compress-text', dest='compress_text', action='store_false', default=True)
    graph = parser.add_argument_group('GRAPH OPTIONS')
    graph.add_argument('-c', '--confidence', dest='confidence', metavar='INT', type=int, help='confidence level (0 lowest to 3 highest) [3]', default=3)
    graph.add_argument('-I', '--iterations', dest='iterations', metavar='INT', type=int, help='number of iterations to insert new introns into the graph [5]', default=5)
    graph.add_argument('-M', '--merge-strat', dest='merge', metavar='<STRAT>', help='merge strategy, where <STRAT> is one of: single, merge_bams, merge_graphs, merge_all [merge_graphs]', default='merge_graphs')
    graph.add_argument('--insert-ir', dest='insert_ir', action='store_true', help='insert intron retentions [on]', default=True)
    graph.add_argument('--no-insert-ir', dest='insert_ir', action='store_false', default=True)
    graph.add_argument('--insert-es', dest='insert_es', action='store_true', help='insert cassette exons [on]', default=True)
    graph.add_argument('--no-insert-es', dest='insert_es', action='store_false', default=True)
    graph.add_argument('--insert-ni', dest='insert_ni', action='store_true', help='insert new intron edges [on]', default=True)
    graph.add_argument('--no-insert-ni', dest='insert_ni', action='store_false', default=True)
    graph.add_argument('--remove-se', dest='remove_se', action='store_true', help='remove short exons [off]', default=False)
    graph.add_argument('--no-remove-se', dest='remove_se', action='store_false', default=False)
    graph.add_argument('--validate-sg',  dest='validate_sg', action='store_true', help='validate splice graph [off]', default=False)
    graph.add_argument('--no-validate-sg', dest='validate_sg', action='store_false', default=False)
    graph.add_argument('--re-infer-sg', dest='infer_sg', action='store_true', help='re-infer splice graph [off] (advanced)', default=False)
    graph.add_argument('--no-re-infer-sg', dest='infer_sg', action='store_false', default=False)
    splice = parser.add_argument_group('SPLICE OPTIONS')
    splice.add_argument('--event-types', dest='event_types', metavar='STRING', help='list of alternative splicing events to extract [exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mtex_exons]', default='exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons')
    splice.add_argument('--extract-ase', dest='extract_as', action='store_true', help='extract alternative splicing events [on]', default=True)
    splice.add_argument('--ase-edge-limit', dest='detect_edge_limit', metavar='INT', help='max number of edges in the graph to still extract events for a gene [500]', default=500)
    splice.add_argument('--no-extract-ase', dest='extract_as', action='store_false', default=True)
    splice.add_argument('--curate-alt-prime', dest='curate_alt_prime', action='store_true', help='curate alt prime events [on]', default=True)
    splice.add_argument('--no-curate-alt-prime', dest='curate_alt_prime', action='store_false', default=True)
    splice.add_argument('--quantify-graph', dest='quantify_graph', action='store_true', help='quantify graph (implicitly on when --extract-ase is used) [off]', default=False)
    splice.add_argument('--no-quantify-graph', dest='quantify_graph', action='store_false', default=False)
    experimental = parser.add_argument_group('EXPERIMENTAL - BETA STATE')
    experimental.add_argument('--pyproc', dest='pyproc', action='store_true', help='use parallel implementation [off]', default=False)
    #experimental.add_argument('-R', '--replicates', dest='replicates', metavar='1,1,2,2,...', help='replicate structure of files (same number as alignment files) [all 1 - no replicated]', default='-')
    experimental.add_argument('--intron-cov', dest='intron_cov', action='store_true', help='count intron coverage [off]', default=False)
    experimental.add_argument('--qmode', dest='qmode', metavar='STRING', help='quantification mode: single, collect, all [all]', default='all')

    options = parser.parse_args(argv[1:])

    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)

    options.parser = parser
    return options


def spladder(argv):

    ### get command line options
    options = parse_options(argv)

    ### parse parameters from options object
    CFG = settings.parse_args(options)

    ### add dependencies provided in config section
    #if 'paths' in CFG:
    #    for i in CFG['paths']:
    #        eval('import %s'% CFG['paths'][i])

    ### load confidence level settings
    if not CFG['no_reset_conf']:
        CFG = settings.set_confidence_level(CFG)

    ### do not compute components of merged set, if result file already exists
    fn_out_merge = get_filename('fn_out_merge', CFG)
    fn_out_merge_val = get_filename('fn_out_merge_val', CFG)

    if not 'spladder_infile' in CFG and not os.path.exists(fn_out_merge):
        ### iterate over files, if merge strategy is single
        if CFG['merge_strategy'] in ['single', 'merge_graphs']:
            idxs = list(range(len(CFG['samples'])))
        else:
            idxs = [0]
        
        ### set parallelization
        if CFG['rproc']:
            jobinfo = []

        ### create out-directory
        if not os.path.exists(CFG['out_dirname']):
            os.makedirs(CFG['out_dirname'])

        ### create spladder sub-directory
        if not os.path.exists(os.path.join(CFG['out_dirname'], 'spladder')):
            os.makedirs(os.path.join(CFG['out_dirname'], 'spladder'))

        ### pre-process annotation, if necessary
        if CFG['anno_fname'].split('.')[-1] != 'pickle':
            if not os.path.exists(CFG['anno_fname'] + '.pickle'):
                if CFG['anno_fname'].split('.')[-1].lower() in ['gff', 'gff3']:
                    (genes, CFG) = init.init_genes_gff3(CFG['anno_fname'], CFG, CFG['anno_fname'] + '.pickle')
                elif CFG['anno_fname'].split('.')[-1].lower() in ['gtf']:
                    (genes, CFG) = init.init_genes_gtf(CFG['anno_fname'], CFG, CFG['anno_fname'] + '.pickle')
                else:
                    print('ERROR: Unknown annotation format. File needs to end in gtf or gff/gff3\nCurrent file: %s' % CFG['anno_fname'], file=sys.stderr)
                    sys.exit(1)
            CFG['anno_fname'] += '.pickle'

        ### add anotation contigs into lookup table
        if not 'genes' in CFG:
            genes = pickle.load(open(CFG['anno_fname'], 'rb'), encoding='latin1')
        else:
            genes = CFG['genes']
        CFG = init.append_chrms(sp.unique(sp.array([x.chr for x in genes], dtype='str')), CFG)
        del genes

        ### convert input BAMs to sparse arrays - filtered case
        if CFG['bam_to_sparse']:
            for bfn in CFG['bam_fnames']:
                if bfn.endswith('bam') and not os.path.exists(re.sub(r'.bam$', '', bfn) + '.conf_%i' % CFG['confidence_level'] + '.filt.hdf5'):
                    #cnts = dict()

                    if not 'chrm_lookup' in CFG:
                        IN = pysam.Samfile(bfn, 'rb')
                        CFG = append_chrms([x['SN'] for x in parse_header(IN.text)['SQ']], CFG)
                        IN.close()

                    OUT = h5py.File(re.sub(r'.bam$', '', bfn) + '.conf_%i' % CFG['confidence_level'] + '.filt.hdf5', 'w')
                    if CFG['parallel'] > 1:
                        import multiprocessing as mp
                        pool = mp.Pool(processes=CFG['parallel'])
                        result = [pool.apply_async(summarize_chr, args=(bfn, str(chrm), CFG,), kwds={'filter':CFG['read_filter']}) for chrm in sorted(CFG['chrm_lookup'])]
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
                        for chrm in CFG['chrm_lookup']:
                            tmp = summarize_chr(bfn, str(chrm), CFG, filter=CFG['read_filter'])
                            OUT.create_dataset(name=(chrm + '_reads_row'), data=tmp[1].row.astype('uint8'), compression='gzip')
                            OUT.create_dataset(name=(chrm + '_reads_col'), data=tmp[1].col, compression='gzip')
                            OUT.create_dataset(name=(chrm + '_reads_dat'), data=tmp[1].data, compression='gzip')
                            OUT.create_dataset(name=(chrm + '_reads_shp'), data=tmp[1].shape)
                            OUT.create_dataset(name=(chrm + '_introns_m'), data=tmp[2], compression='gzip')
                            OUT.create_dataset(name=(chrm + '_introns_p'), data=tmp[3], compression='gzip')
                    OUT.close()
                elif CFG['verbose']:
                    print('Filtered sparse BAM representation for %s already exists.' % bfn, file=sys.stdout)

        ### build individual graphs
        for idx in idxs:
            CFG_ = dict()
            if CFG['merge_strategy'] != 'merge_bams':
                CFG_['bam_fnames'] = CFG['bam_fnames']
                CFG_['samples'] = CFG['samples']
                CFG['bam_fnames'] = CFG['bam_fnames'][idx]
                CFG['samples'] = CFG['samples'][idx]
                CFG['out_fname'] = '%s/spladder/genes_graph_conf%i.%s.pickle' % (CFG['out_dirname'], CFG['confidence_level'], CFG['samples'])
            else:
                CFG['out_fname'] = '%s/spladder/genes_graph_conf%i.%s.pickle' % (CFG['out_dirname'], CFG['confidence_level'], CFG['merge_strategy'])

            ### assemble out filename to check if we are already done
            fn_out = CFG['out_fname']
            if CFG['do_prune']:
                fn_out = re.sub('.pickle$', '_pruned.pickle', fn_out)
            if CFG['do_gen_isoforms']:
                fn_out = re.sub('.pickle$', '_with_isoforms.pickle', fn_out)
    
            if os.path.exists(fn_out):
                print('%s - All result files already exist.' % fn_out, file=sys.stdout)
            else:
                if CFG['rproc']:
                    jobinfo.append(rp.rproc('spladder_core', CFG, 15000, CFG['options_rproc'], 60*60))
                else:
                    spladder_core(CFG)

            for key in CFG_:
                try:
                    CFG[key] = CFG_[key].copy()
                except AttributeError:
                    CFG[key] = CFG_[key]

        ### collect results after parallelization
        if CFG['rproc']:
            rp.rproc_wait(jobinfo, 30, 1.0, -1)

        ### merge parts if necessary
        if CFG['merge_strategy'] == 'merge_graphs':
            run_merge(CFG)

    if not 'spladder_infile' in CFG and CFG['merge_strategy'] == 'merge_graphs' and CFG['validate_splicegraphs'] and not os.path.exists(fn_out_merge_val):
        (genes, inserted) = pickle.load(open(fn_out_merge, 'rb'))
        genes = filter_by_edgecount(genes, CFG)
        pickle.dump((genes, inserted), open(fn_out_merge_val, 'wb'), -1)
        del genes

    ### convert input BAMs to sparse arrays - unfiltered case
    if CFG['bam_to_sparse']:
        for bfn in CFG['bam_fnames']:
            if bfn.endswith('bam') and not os.path.exists(re.sub(r'.bam$', '', bfn) + '.hdf5'):
                #cnts = dict()

                if not 'chrm_lookup' in CFG:
                    IN = pysam.Samfile(bfn, 'rb')
                    CFG = append_chrms([x['SN'] for x in parse_header(IN.text)['SQ']], CFG)
                    IN.close()

                OUT = h5py.File(re.sub(r'.bam$', '', bfn) + '.hdf5', 'w') 
                if CFG['parallel'] > 1:
                    import multiprocessing as mp
                    pool = mp.Pool(processes=CFG['parallel'])
                    result = [pool.apply_async(summarize_chr, args=(bfn, str(chrm), CFG,)) for chrm in sorted(CFG['chrm_lookup'])]
                    while result:
                        tmp = result.pop(0).get()
                        OUT.create_dataset(name=(tmp[0] + '_reads_row'), data=tmp[1].row.astype('uint8'), compression='gzip')
                        OUT.create_dataset(name=(tmp[0] + '_reads_col'), data=tmp[1].col, compression='gzip')
                        OUT.create_dataset(name=(tmp[0] + '_reads_dat'), data=tmp[1].data, compression='gzip')
                        OUT.create_dataset(name=(tmp[0] + '_reads_shp'), data=tmp[1].shape)
                        OUT.create_dataset(name=(tmp[0] + '_introns_m'), data=tmp[2], compression='gzip')
                        OUT.create_dataset(name=(tmp[0] + '_introns_p'), data=tmp[3], compression='gzip')
                else:
                    for chrm in CFG['chrm_lookup']:
                        tmp = summarize_chr(bfn, str(chrm), CFG)
                        OUT.create_dataset(name=(chrm + '_reads_row'), data=tmp[1].row.astype('uint8'), compression='gzip')
                        OUT.create_dataset(name=(chrm + '_reads_col'), data=tmp[1].col, compression='gzip')
                        OUT.create_dataset(name=(chrm + '_reads_dat'), data=tmp[1].data, compression='gzip')
                        OUT.create_dataset(name=(chrm + '_reads_shp'), data=tmp[1].shape)
                        OUT.create_dataset(name=(chrm + '_introns_m'), data=tmp[2], compression='gzip')
                        OUT.create_dataset(name=(chrm + '_introns_p'), data=tmp[3], compression='gzip')
                OUT.close()
            elif CFG['verbose']:
                print('Sparse BAM representation for %s already exists.' % bfn, file=sys.stdout)

    if CFG['merge_strategy'] == 'single' or CFG['quantification_mode'] == 'collect':
        idxs = list(range(len(CFG['samples'])))
    else:
        idxs = [0]

    for idx in idxs:
        ### get count output file
        if CFG['merge_strategy'] == 'single':
            fn_in_count = get_filename('fn_count_in', CFG, sample_idx=idx)
            fn_out_count = get_filename('fn_count_out', CFG, sample_idx=idx)
        elif CFG['merge_strategy'] == 'merge_graphs' and CFG['quantification_mode'] == 'single':
            fn_in_count = get_filename('fn_count_in', CFG)
            fn_out_count = get_filename('fn_count_out', CFG, sample_idx=0)
        else:
            if CFG['quantification_mode'] != 'collect':
                fn_in_count = get_filename('fn_count_in', CFG)
            fn_out_count = get_filename('fn_count_out', CFG)

        ### count segment graph
        if CFG['run_as_analysis'] or CFG['count_segment_graph']:
            if not os.path.exists(fn_out_count):
                if CFG['merge_strategy'] == 'single':
                    count_graph_coverage_wrapper(fn_in_count, fn_out_count, CFG, sample_idx=idx)
                elif CFG['merge_strategy'] == 'merge_graphs' and CFG['quantification_mode'] == 'single':
                    count_graph_coverage_wrapper(fn_in_count, fn_out_count, CFG, qmode='single')
                elif CFG['merge_strategy'] == 'merge_graphs' and CFG['quantification_mode'] == 'collect':
                    collect_single_quantification_results(fn_out_count, idxs, CFG)
                else:
                    count_graph_coverage_wrapper(fn_in_count, fn_out_count, CFG)

        ### count intron coverage phenotype
        if CFG['count_intron_cov']:
            fn_out_intron_count = fn_out_count.replace('pickle', 'introns.pickle')
            count_intron_coverage_wrapper(fn_in_count, fn_out_intron_count, CFG)

    ### handle alternative splicing part
    if CFG['run_as_analysis']:
        collect_events(CFG)

        for idx in idxs:
            for e_idx in range(len(CFG['event_types'])):
                if CFG['merge_strategy'] == 'single':
                    analyze_events(CFG, CFG['event_types'][e_idx], sample_idx=idx)
                else:
                    analyze_events(CFG, CFG['event_types'][e_idx])

def main(argv):
    spladder(argv)

if __name__ == "__main__":
    main(sys.argv)
