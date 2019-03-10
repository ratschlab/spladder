import pickle
import scipy as sp
import os
import pdb


if __package__ is None:
    __package__ = 'modules.alt_splice'

from .events import *
from .detect import *

from ..init import append_chrms
from ..classes.event import Event


def collect_events(options):

    ### which events do we call
    do_exon_skip = ('exon_skip' in options.event_types)
    do_intron_retention = ('intron_retention' in options.event_types)
    do_mult_exon_skip = ('mult_exon_skip' in options.event_types)
    do_alt_3prime = ('alt_3prime' in options.event_types)
    do_alt_5prime = ('alt_5prime' in options.event_types)
    do_mutex_exons = ('mutex_exons' in options.event_types)

    ### init empty event fields
    intron_reten_pos = []
    exon_skip_pos = []
    mult_exon_skip_pos = []
    alt_end_5prime_pos = [] 
    alt_end_3prime_pos = []
    mutex_exons_pos = []

    validate_tag = ''
    if hasattr(options, 'validate_sg') and options.validate_sg:
        validate_tag = '.validated'

    for i in range(len(options.samples)):
        if i == 1:
            break

        strain = options.strains[i]

        if hasattr(options, 'spladderfile') and options.spladderfile != '-':
            genes_fnames = options.spladderfile
        elif options.merge == 'single':
            genes_fnames = '%s/spladder/genes_graph_conf%i.%s.pickle' % (options.outdir, options.confidence, options.samples[i])
        else:
            genes_fnames = '%s/spladder/genes_graph_conf%i.%s%s.pickle' % (options.outdir, options.confidence, options.merge, validate_tag)

        ### define outfile names
        sample_tag = options.merge
        if options.merge == 'single':
            sample_tag = options.samples[i]
        fn_out_ir = '%s/%s_intron_retention_C%i.pickle' % (options.outdir, sample_tag, options.confidence)
        fn_out_es = '%s/%s_exon_skip_C%i.pickle' % (options.outdir, sample_tag, options.confidence)
        fn_out_mes = '%s/%s_mult_exon_skip_C%i.pickle' % (options.outdir, sample_tag, options.confidence) 
        fn_out_a5 = '%s/%s_alt_5prime_C%i.pickle' % (options.outdir, sample_tag, options.confidence)
        fn_out_a3 = '%s/%s_alt_3prime_C%i.pickle' % (options.outdir, sample_tag, options.confidence)
        fn_out_mex = '%s/%s_mutex_exons_C%i.pickle' % (options.outdir, sample_tag, options.confidence)

        print('\nconfidence %i / sample %i' % (options.confidence, i))

        if os.path.exists(genes_fnames):
            print('Loading gene structure from %s ...' % genes_fnames)
            (genes, inserted) = pickle.load(open(genes_fnames, 'rb'))
            print('... done.')
            
            if not hasattr(options, 'chrm_lookup'):
                options = append_chrms(sp.unique(sp.array([x.chr for x in genes], dtype='str')), options)

            ### detect intron retentions from splicegraph
            if do_intron_retention:
                if not os.path.exists(fn_out_ir):
                    idx_intron_reten, intron_intron_reten = detect_events(genes, 'intron_retention', sp.where([x.is_alt for x in genes])[0], options)
                    for k in range(len(idx_intron_reten)):
                        gene = genes[idx_intron_reten[k]]

                        ### perform liftover between strains if necessary
                        exons = gene.splicegraph.vertices
                        if not hasattr(options, 'reference_strain'):
                            exons_col = exons
                            exons_col_pos = exons
                        else:
                            exons_col = convert_strain_pos_intervals(gene.chr, gene.splicegraph.vertices.T, strain, options.reference_strain).T
                            exons_col_pos = convert_strain_pos(gene.chr, gene.splicegraph.vertices.T, strain, options.reference_strain).T
                        if exons_col.shape != exons_col_pos.shape: 
                            print('skipping non-mappable intron retention event')
                            continue

                        ### build intron retention data structure
                        event = Event('intron_retention', gene.chr, gene.strand)
                        event.strain = sp.array([strain])
                        event.exons1 = sp.c_[exons[:, intron_intron_reten[k][0]], exons[:, intron_intron_reten[k][1]]].T
                        event.exons2 = sp.array([exons[:, intron_intron_reten[k][0]][0], exons[:, intron_intron_reten[k][1]][1]])
                        #event.exons2 = exons[:, intron_intron_reten[k][2]]
                        event.exons1_col = sp.c_[exons_col[:, intron_intron_reten[k][0]], exons_col[:, intron_intron_reten[k][1]]]
                        event.exons2_col = sp.array([exons_col[:, intron_intron_reten[k][0]][0], exons_col[:, intron_intron_reten[k][1]][1]])
                        #event.exons2_col = exons_col[:, intron_intron_reten[k][2]]
                        event.gene_name = sp.array([gene.name])
                        event.gene_idx = idx_intron_reten[k]
                        #event.transcript_type = sp.array([gene.transcript_type])
                        intron_reten_pos.append(event)
                else:
                    print('%s already exists' % fn_out_ir)

            ### detect exon_skips from splicegraph
            if do_exon_skip:
                if not os.path.exists(fn_out_es):
                    idx_exon_skip, exon_exon_skip = detect_events(genes, 'exon_skip', sp.where([x.is_alt for x in genes])[0], options)
                    for k in range(len(idx_exon_skip)):
                        gene = genes[idx_exon_skip[k]]

                        ### perform liftover between strains if necessary
                        exons = gene.splicegraph.vertices
                        if not hasattr(options, 'reference_strain'):
                            exons_col = exons
                            exons_col_pos = exons
                        else:
                            exons_col = convert_strain_pos_intervals(gene.chr, gene.splicegraph.vertices.T, strain, options.reference_strain).T
                            exons_col_pos = convert_strain_pos(gene.chr, gene.splicegraph.vertices.T, strain, options.reference_strain).T
                        if exons_col.shape != exons_col_pos.shape: 
                            print('skipping non-mappable exon_skip event')
                            continue

                        ### build exon skip data structure
                        event = Event('exon_skip', gene.chr, gene.strand)
                        event.strain = sp.array([strain])
                        event.exons1 = sp.c_[exons[:, exon_exon_skip[k][0]], exons[:, exon_exon_skip[k][2]]].T
                        event.exons2 = sp.c_[exons[:, exon_exon_skip[k][0]], exons[:, exon_exon_skip[k][1]], exons[:, exon_exon_skip[k][2]]].T
                        event.exons1_col = sp.c_[exons_col[:, exon_exon_skip[k][0]], exons_col[:, exon_exon_skip[k][2]]].T
                        event.exons2_col = sp.c_[exons_col[:, exon_exon_skip[k][0]], exons_col[:, exon_exon_skip[k][1]], exons_col[:, exon_exon_skip[k][2]]].T
                        event.gene_name = sp.array([gene.name])
                        event.gene_idx = idx_exon_skip[k]
                        #event.transcript_type = sp.array([gene.transcript_type])
                        exon_skip_pos.append(event)
                else:
                    print('%s already exists' % fn_out_es)

            ### detect alternative intron_ends from splicegraph
            if do_alt_3prime or do_alt_5prime:
                if not os.path.exists(fn_out_a5) or not os.path.exists(fn_out_a3):
                    idx_alt_end_5prime, exon_alt_end_5prime, idx_alt_end_3prime, exon_alt_end_3prime = detect_events(genes, 'alt_prime', sp.where([x.is_alt for x in genes])[0], options)
                    ### handle 5 prime events
                    for k in range(len(idx_alt_end_5prime)):
                        gene = genes[idx_alt_end_5prime[k]]

                        ### perform liftover between strains if necessary
                        exons = gene.splicegraph.vertices
                        if not hasattr(options, 'reference_strain'):
                            exons_col = exons
                            exons_col_pos = exons
                        else:
                            exons_col = convert_strain_pos_intervals(gene.chr, gene.splicegraph.vertices.T, strain, options.reference_strain).T
                            exons_col_pos = convert_strain_pos(gene.chr, gene.splicegraph.vertices.T, strain, options.reference_strain).T
                        if exons_col.shape != exons_col_pos.shape: 
                            print('skipping non-mappable alt 5-prime event')
                            continue
                        
                        for k1 in range(len(exon_alt_end_5prime[k]['fiveprimesites']) - 1):
                            for k2 in range(k1 + 1, len(exon_alt_end_5prime[k]['fiveprimesites'])):

                                exon_alt1_col = exons_col[:, exon_alt_end_5prime[k]['fiveprimesites'][k1]].T
                                exon_alt2_col = exons_col[:, exon_alt_end_5prime[k]['fiveprimesites'][k2]].T

                                ### check if exons overlap
                                if (exon_alt1_col[0] >= exon_alt2_col[1]) or (exon_alt1_col[1] <= exon_alt2_col[0]):
                                    continue

                                event = Event('alt_5prime', gene.chr, gene.strand)
                                event.strain = sp.array([strain])
                                if gene.strand == '+':
                                    event.exons1 = sp.c_[exons[:, exon_alt_end_5prime[k]['fiveprimesites'][k1]], exons[:, exon_alt_end_5prime[k]['threeprimesite']]].T
                                    event.exons2 = sp.c_[exons[:, exon_alt_end_5prime[k]['fiveprimesites'][k2]], exons[:, exon_alt_end_5prime[k]['threeprimesite']]].T
                                    event.exons1_col = sp.c_[exons_col[:, exon_alt_end_5prime[k]['fiveprimesites'][k1]], exons_col[:, exon_alt_end_5prime[k]['threeprimesite']]].T
                                    event.exons2_col = sp.c_[exons_col[:, exon_alt_end_5prime[k]['fiveprimesites'][k2]], exons_col[:, exon_alt_end_5prime[k]['threeprimesite']]].T
                                else:
                                    event.exons1 = sp.c_[exons[:, exon_alt_end_5prime[k]['threeprimesite']], exons[:, exon_alt_end_5prime[k]['fiveprimesites'][k1]]].T
                                    event.exons2 = sp.c_[exons[:, exon_alt_end_5prime[k]['threeprimesite']], exons[:, exon_alt_end_5prime[k]['fiveprimesites'][k2]]].T
                                    event.exons1_col = sp.c_[exons_col[:, exon_alt_end_5prime[k]['threeprimesite']], exons_col[:, exon_alt_end_5prime[k]['fiveprimesites'][k1]]].T
                                    event.exons2_col = sp.c_[exons_col[:, exon_alt_end_5prime[k]['threeprimesite']], exons_col[:, exon_alt_end_5prime[k]['fiveprimesites'][k2]]].T
                                event.gene_name = sp.array([gene.name])
                                event.gene_idx = idx_alt_end_5prime[k]

                                ### assert that first isoform is always the shorter one
                                if sp.sum(event.exons1[:, 1] - event.exons1[:, 0]) > sp.sum(event.exons2[:, 1] - event.exons2[:, 0]):
                                    _tmp = event.exons1.copy()
                                    event.exons1 = event.exons2.copy()
                                    event.exons2 = _tmp
                                #event.transcript_type = sp.array([gene.transcript_type])
                                if do_alt_5prime:
                                    alt_end_5prime_pos.append(event)

                    ### handle 3 prime events
                    for k in range(len(idx_alt_end_3prime)):
                        gene = genes[idx_alt_end_3prime[k]]

                        ### perform liftover between strains if necessary
                        exons = gene.splicegraph.vertices
                        if not hasattr(options, 'reference_strain'):
                            exons_col = exons
                            exons_col_pos = exons
                        else:
                            exons_col = convert_strain_pos_intervals(gene.chr, gene.splicegraph.vertices.T, strain, options.reference_strain).T
                            exons_col_pos = convert_strain_pos(gene.chr, gene.splicegraph.vertices.T, strain, options.reference_strain).T
                        if exons_col.shape != exons_col_pos.shape: 
                            print('skipping non-mappable alt 3-prime event')
                            continue

                        for k1 in range(len(exon_alt_end_3prime[k]['threeprimesites']) - 1):
                            for k2 in range(k1 + 1, len(exon_alt_end_3prime[k]['threeprimesites'])):

                                exon_alt1_col = exons_col[:, exon_alt_end_3prime[k]['threeprimesites'][k1]].T
                                exon_alt2_col = exons_col[:, exon_alt_end_3prime[k]['threeprimesites'][k2]].T

                                ### check if exons overlap
                                if (exon_alt1_col[0] >= exon_alt2_col[1]) or (exon_alt1_col[1] <= exon_alt2_col[0]):
                                    continue

                                event = Event('alt_3prime', gene.chr, gene.strand)
                                event.strain = sp.array([strain])
                                if gene.strand == '+':
                                    event.exons1 = sp.c_[exons[:, exon_alt_end_3prime[k]['threeprimesites'][k1]], exons[:, exon_alt_end_3prime[k]['fiveprimesite']]].T
                                    event.exons2 = sp.c_[exons[:, exon_alt_end_3prime[k]['threeprimesites'][k2]], exons[:, exon_alt_end_3prime[k]['fiveprimesite']]].T
                                    event.exons1_col = sp.c_[exons_col[:, exon_alt_end_3prime[k]['threeprimesites'][k1]], exons_col[:, exon_alt_end_3prime[k]['fiveprimesite']]].T
                                    event.exons2_col = sp.c_[exons_col[:, exon_alt_end_3prime[k]['threeprimesites'][k2]], exons_col[:, exon_alt_end_3prime[k]['fiveprimesite']]].T
                                else:
                                    event.exons1 = sp.c_[exons[:, exon_alt_end_3prime[k]['fiveprimesite']], exons[:, exon_alt_end_3prime[k]['threeprimesites'][k1]]].T
                                    event.exons2 = sp.c_[exons[:, exon_alt_end_3prime[k]['fiveprimesite']], exons[:, exon_alt_end_3prime[k]['threeprimesites'][k2]]].T
                                    event.exons1_col = sp.c_[exons_col[:, exon_alt_end_3prime[k]['fiveprimesite']], exons_col[:, exon_alt_end_3prime[k]['threeprimesites'][k1]]].T
                                    event.exons2_col = sp.c_[exons_col[:, exon_alt_end_3prime[k]['fiveprimesite']], exons_col[:, exon_alt_end_3prime[k]['threeprimesites'][k2]]].T
                                event.gene_name = sp.array([gene.name])
                                event.gene_idx = idx_alt_end_3prime[k]

                                ### assert that first isoform is always the shorter one
                                if sp.sum(event.exons1[:, 1] - event.exons1[:, 0]) > sp.sum(event.exons2[:, 1] - event.exons2[:, 0]):
                                    _tmp = event.exons1.copy()
                                    event.exons1 = event.exons2.copy()
                                    event.exons2 = _tmp

                                #event.transcript_type = sp.array([gene.transcript_type])
                                if do_alt_3prime:
                                    alt_end_3prime_pos.append(event)
                else:
                    print('%s and %s already exists' % (fn_out_a5, fn_out_a3))

            ### detect multiple_exon_skips from splicegraph
            if do_mult_exon_skip:
                if not os.path.exists(fn_out_mes):
                    idx_mult_exon_skip, exon_mult_exon_skip = detect_events(genes, 'mult_exon_skip', sp.where([x.is_alt for x in genes])[0], options)
                    for k, gidx in enumerate(idx_mult_exon_skip):
                        gene = genes[gidx] 

                        ### perform liftover between strains if necessary
                        exons = gene.splicegraph.vertices
                        if not hasattr(options, 'reference_strain'):
                            exons_col = exons
                            exons_col_pos = exons
                        else:
                            exons_col = convert_strain_pos_intervals(gene.chr, gene.splicegraph.vertices.T, strain, options.reference_strain).T
                            exons_col_pos = convert_strain_pos(gene.chr, gene.splicegraph.vertices.T, strain, options.reference_strain).T
                        if exons_col.shape != exons_col_pos.shape: 
                            print('skipping non-mappable multiple exon skip event')
                            continue

                        ### build multiple exon skip data structure
                        event = Event('mult_exon_skip', gene.chr, gene.strand)
                        event.strain = sp.array([strain])
                        event.exons1 = sp.c_[exons[:, exon_mult_exon_skip[k][0]], exons[:, exon_mult_exon_skip[k][2]]].T
                        event.exons2 = sp.c_[exons[:, exon_mult_exon_skip[k][0]], exons[:, exon_mult_exon_skip[k][1]], exons[:, exon_mult_exon_skip[k][2]]].T
                        event.exons1_col = sp.c_[exons_col[:, exon_mult_exon_skip[k][0]], exons_col[:, exon_mult_exon_skip[k][2]]].T
                        event.exons2_col = sp.c_[exons_col[:, exon_mult_exon_skip[k][0]], exons_col[:, exon_mult_exon_skip[k][1]], exons_col[:, exon_mult_exon_skip[k][2]]].T
                        event.gene_name = sp.array([gene.name])
                        event.gene_idx = gidx
                        #event.transcript_type = sp.array([gene.transcript_type])
                        mult_exon_skip_pos.append(event)
                else:
                    print('%s already exists' % fn_out_mes)

            ### detect mutually exclusive exons from splicegraph
            if do_mutex_exons:
                if not os.path.exists(fn_out_mex):
                    idx_mutex_exons, exon_mutex_exons = detect_events(genes, 'mutex_exons', sp.where([x.is_alt for x in genes])[0], options)
                    if len(idx_mutex_exons) > 0:
                        for k in range(len(exon_mutex_exons)):
                            gene = genes[idx_mutex_exons[k]]

                            ### perform liftover between strains if necessary
                            exons = gene.splicegraph.vertices
                            if not hasattr(options, 'reference_strain'):
                                exons_col = exons
                                exons_col_pos = exons
                            else:
                                exons_col = convert_strain_pos_intervals(gene.chr, gene.splicegraph.vertices.T, strain, options.reference_strain).T
                                exons_col_pos = convert_strain_pos(gene.chr, gene.splicegraph.vertices.T, strain, options.reference_strain).T

                            if exons_col.shape != exons_col_pos.shape: 
                                print('skipping non-mappable mutex exons event')
                                continue

                            ### build data structure for mutually exclusive exons
                            event = Event('mutex_exons', gene.chr, gene.strand)
                            event.strain = sp.array([strain])
                            event.exons1 = sp.c_[exons[:, exon_mutex_exons[k][0]], exons[:, exon_mutex_exons[k][1]], exons[:, exon_mutex_exons[k][3]]].T
                            event.exons2 = sp.c_[exons[:, exon_mutex_exons[k][0]], exons[:, exon_mutex_exons[k][2]], exons[:, exon_mutex_exons[k][3]]].T
                            event.exons1_col = sp.c_[exons_col[:, exon_mutex_exons[k][0]], exons_col[:, exon_mutex_exons[k][1]], exons_col[:, exon_mutex_exons[k][3]]].T
                            event.exons2_col = sp.c_[exons_col[:, exon_mutex_exons[k][0]], exons_col[:, exon_mutex_exons[k][2]], exons_col[:, exon_mutex_exons[k][3]]].T
                            event.gene_name = sp.array([gene.name])
                            event.gene_idx = idx_mutex_exons[k]
                            #event.transcript_type = sp.array([gene.transcript_type])
                            mutex_exons_pos.append(event)
                else:
                    print('%s already exists' % fn_out_mex)

        ### genes file does not exist
        else:
            print('result file not found: %s' % genes_fnames)

    ################################################%
    ### COMBINE INTRON RETENTIONS
    ################################################%
    if do_intron_retention:
        if not os.path.exists(fn_out_ir):

            ### post process event structure by sorting and making events unique
            events_all = post_process_event_struct(sp.array(intron_reten_pos), options)

            ### store intron retentions
            print('saving intron retentions to %s' % fn_out_ir)
            pickle.dump(events_all, open(fn_out_ir, 'wb'), -1)
        else:
            print('%s already exists' % fn_out_ir)
    
    ################################################%
    ### COMBINE EXON SKIPS
    ################################################%
    if do_exon_skip:
        if not os.path.exists(fn_out_es):

            ### post process event structure by sorting and making events unique
            events_all = post_process_event_struct(sp.array(exon_skip_pos), options)

            ### store exon skip events
            print('saving exon skips to %s' % fn_out_es)
            pickle.dump(events_all, open(fn_out_es, 'wb'), -1)
        else:
            print('%s already exists' % fn_out_es)

    ################################################%
    ### COMBINE MULTIPLE EXON SKIPS
    ################################################%
    if do_mult_exon_skip:
        if not os.path.exists(fn_out_mes):

            ### post process event structure by sorting and making events unique
            events_all = post_process_event_struct(sp.array(mult_exon_skip_pos), options)

            ### store multiple exon skip events
            print('saving multiple exon skips to %s' % fn_out_mes)
            pickle.dump(events_all, open(fn_out_mes, 'wb'), -1)
        else:
            print('%s already exists' % fn_out_mes)

    ################################################%
    ### COMBINE ALT FIVE PRIME EVENTS
    ################################################%
    if do_alt_5prime:
        if not os.path.exists(fn_out_a5):
          
            ### post process event structure by sorting and making events unique
            events_all = post_process_event_struct(sp.array(alt_end_5prime_pos), options)

            ### curate alt prime events
            ### cut to min len, if alt exon lengths differ
            ### remove, if alt exons do not overlap
            if options.curate_alt_prime:
                events_all = curate_alt_prime(events_all, options)

            ### store alt 5 prime events
            print('saving alt 5 prime events to %s' % fn_out_a5)
            pickle.dump(events_all, open(fn_out_a5, 'wb'), -1)
        else:
            print('%s already exists' % fn_out_a5)

    ################################################%
    ### COMBINE ALT THREE PRIME EVENTS
    ################################################%
    if do_alt_3prime:
        if not os.path.exists(fn_out_a3):

            ### post process event structure by sorting and making events unique
            events_all = post_process_event_struct(sp.array(alt_end_3prime_pos), options)

            ### curate alt prime events
            ### cut to min len, if alt exon lengths differ
            ### remove, if alt exons do not overlap
            if options.curate_alt_prime:
                events_all = curate_alt_prime(events_all, options)

            ### store alt 3 prime events
            print('saving alt 3 prime events to %s' % fn_out_a3)
            pickle.dump(events_all, open(fn_out_a3, 'wb'), -1)
        else:
            print('%s already exists' % fn_out_a3)

    ################################################%
    ### COMBINE MUTUALLY EXCLUSIVE EXONS
    ################################################%
    if do_mutex_exons:
        if not os.path.exists(fn_out_mex):

            ### post process event structure by sorting and making events unique
            events_all = post_process_event_struct(sp.array(mutex_exons_pos), options)

            ### store multiple exon skip events
            print('saving mutually exclusive exons to %s' % fn_out_mex)
            pickle.dump(events_all, open(fn_out_mex, 'wb'), -1)
        else:
            print('%s already exists' % fn_out_mex)


