import cPickle
import scipy as sp
import os
import pdb

from .events import *
from .detect import *

from ..classes.event import Event


def collect_events(CFG):

    ### which events do we call
    do_exon_skip = ('exon_skip' in CFG['event_types'])
    do_intron_retention = ('intron_retention' in CFG['event_types'])
    do_mult_exon_skip = ('mult_exon_skip' in CFG['event_types'])
    do_alt_3prime = ('alt_3prime' in CFG['event_types'])
    do_alt_5prime = ('alt_5prime' in CFG['event_types'])

    ### init empty event fields
    if do_intron_retention:
        intron_reten_pos = sp.zeros((len(CFG['replicate_idxs']), 1), dtype = 'object')
    if do_exon_skip:
        exon_skip_pos = sp.zeros((len(CFG['replicate_idxs']), 1), dtype = 'object')
    if do_alt_3prime or do_alt_5prime:
        alt_end_5prime_pos = sp.zeros((len(CFG['replicate_idxs']), 1), dtype = 'object')
        alt_end_3prime_pos = sp.zeros((len(CFG['replicate_idxs']), 1), dtype = 'object')
    if do_mult_exon_skip:
        mult_exon_skip_pos = sp.zeros((len(CFG['replicate_idxs']), 1), dtype = 'object')

    validate_tag = ''
    if 'validate_splicegraphs' in CFG and CFG['validate_splicegraphs']:
        validate_tag = '.validated'

    ### set offset that can be used for half open intervals
    if CFG['is_half_open']:
        ho_offset = 1
    else:
        ho_offset = 0

    for i in range(len(CFG['samples'])):
        if CFG['same_genestruct_for_all_samples'] == 1 and i == 1:
            break

        if i > 0:
            if do_intron_retention:
                intron_reten_pos = sp.c_[intron_reten_pos, sp.zeros((len(CFG['replicate_idxs']), 1), dtype = 'object')]
            if do_exon_skip:
                exon_skip_pos = sp.c_[exon_skip_pos, sp.zeros((len(CFG['replicate_idxs']), 1), dtype = 'object')]
            if do_alt_3prime or do_alt_5prime:
                alt_end_5prime_pos = sp.c_[alt_end_5prime_pos, sp.zeros((len(CFG['replicate_idxs']), 1), dtype = 'object')]
                alt_end_3prime_pos = sp.c_[alt_end_3prime_pos, sp.zeros((len(CFG['replicate_idxs']), 1), dtype = 'object')]
            if do_mult_exon_skip:
                mult_exon_skip_pos = sp.c_[mult_exon_skip_pos, sp.zeros((len(CFG['replicate_idxs']), 1), dtype = 'object')]

        strain = CFG['strains'][i]

        for ridx in CFG['replicate_idxs']:
            if len(CFG['replicate_idxs']) > 1:
                rep_tag = '_R%i' % ridx
            else:
                rep_tag = ''

            if 'spladder_infile' in CFG:
                genes_fnames = CFG['spladder_infile']
            elif CFG['merge_strategy'] == 'single':
                genes_fnames = '%s/spladder/genes_graph_conf%i%s.%s.pickle' % (CFG['out_dirname'], CFG['confidence_level'], rep_tag, CFG['samples'][i])
            else:
                genes_fnames = '%s/spladder/genes_graph_conf%i%s.%s%s.pickle' % (CFG['out_dirname'], CFG['confidence_level'], rep_tag, CFG['merge_strategy'], validate_tag)

            ### define outfile names
            fn_out_ir = '%s/%s_intron_retention%s_C%i.pickle' % (CFG['out_dirname'], CFG['merge_strategy'], rep_tag, CFG['confidence_level'])
            fn_out_es = '%s/%s_exon_skip%s_C%i.pickle' % (CFG['out_dirname'], CFG['merge_strategy'], rep_tag, CFG['confidence_level'])
            fn_out_mes = '%s/%s_mult_exon_skip%s_C%i.pickle' % (CFG['out_dirname'], CFG['merge_strategy'], rep_tag, CFG['confidence_level']) 
            fn_out_a5 = '%s/%s_alt_5prime%s_C%i.pickle' % (CFG['out_dirname'], CFG['merge_strategy'], rep_tag, CFG['confidence_level'])
            fn_out_a3 = '%s/%s_alt_3prime%s_C%i.pickle' % (CFG['out_dirname'], CFG['merge_strategy'], rep_tag, CFG['confidence_level'])

            intron_reten_pos[ridx, i] = []
            exon_skip_pos[ridx, i] = []
            mult_exon_skip_pos[ridx, i] = []
            alt_end_5prime_pos[ridx, i] = []
            alt_end_3prime_pos[ridx, i] = []

            print '\nconfidence %i / sample %i / replicate %i' % (CFG['confidence_level'], i, ridx)

            if os.path.exists(genes_fnames):
                print 'Loading gene structure from %s ...' % genes_fnames
                (genes, inserted) = cPickle.load(open(genes_fnames, 'r'))
                print '... done.'

                ### detect intron retentions from splicegraph
                if do_intron_retention:
                    if not os.path.exists(fn_out_ir):
                        idx_intron_reten, intron_intron_reten = detect_intronreten(genes, sp.where([x.is_alt for x in genes])[0])
                        for k in range(len(idx_intron_reten)):
                            gene = genes[idx_intron_reten[k]]

                            ### perform liftover between strains if necessary
                            exons = gene.splicegraph.vertices
                            if not 'reference_strain' in CFG:
                                exons_col = exons
                                exons_col_pos = exons
                            else:
                                exons_col = convert_strain_pos_intervals(gene.chr_num, gene.splicegraph.vertices.T, strain, CFG['reference_strain']).T
                                exons_col_pos = convert_strain_pos(gene.chr_num, gene.splicegraph.vertices.T, strain, CFG['reference_strain']).T
                            if exons_col.shape != exons_col_pos.shape: 
                                print 'skipping non-mappable intron retention event'
                                continue

                            ### build intron retention data structure
                            event = Event('intron_retention', gene.chr, gene.chr_num, gene.strand)
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
                            intron_reten_pos[ridx, i].append(event)
                    else:
                        print '%s already exists' % fn_out_ir

                ### detect exon_skips from splicegraph
                if do_exon_skip:
                    if not os.path.exists(fn_out_es):
                        idx_exon_skip, exon_exon_skip = detect_exonskips(genes, sp.where([x.is_alt for x in genes])[0])
                        for k in range(len(idx_exon_skip)):
                            gene = genes[idx_exon_skip[k]]

                            ### perform liftover between strains if necessary
                            exons = gene.splicegraph.vertices
                            if not 'reference_strain' in CFG:
                                exons_col = exons
                                exons_col_pos = exons
                            else:
                                exons_col = convert_strain_pos_intervals(gene.chr_num, gene.splicegraph.vertices.T, strain, CFG['reference_strain']).T
                                exons_col_pos = convert_strain_pos(gene.chr_num, gene.splicegraph.vertices.T, strain, CFG['reference_strain']).T
                            if exons_col.shape != exons_col_pos.shape: 
                                print 'skipping non-mappable exon_skip event'
                                continue

                            ### build exon skip data structure
                            event = Event('exon_skip', gene.chr, gene.chr_num, gene.strand)
                            event.strain = sp.array([strain])
                            event.exons1 = sp.c_[exons[:, exon_exon_skip[k][0]], exons[:, exon_exon_skip[k][2]]].T
                            event.exons2 = sp.c_[exons[:, exon_exon_skip[k][0]], exons[:, exon_exon_skip[k][1]], exons[:, exon_exon_skip[k][2]]].T
                            event.exons1_col = sp.c_[exons_col[:, exon_exon_skip[k][0]], exons_col[:, exon_exon_skip[k][2]]].T
                            event.exons2_col = sp.c_[exons_col[:, exon_exon_skip[k][0]], exons_col[:, exon_exon_skip[k][1]], exons_col[:, exon_exon_skip[k][2]]].T
                            event.gene_name = sp.array([gene.name])
                            event.gene_idx = idx_exon_skip[k]
                            #event.transcript_type = sp.array([gene.transcript_type])
                            exon_skip_pos[ridx, i].append(event)
                    else:
                        print '%s already exists' % fn_out_es

                ### detect alternative intron_ends from splicegraph
                if do_alt_3prime or do_alt_5prime:
                    if not os.path.exists(fn_out_a5) or not os.path.exists(fn_out_a3):
                        idx_alt_end_5prime, exon_alt_end_5prime, idx_alt_end_3prime, exon_alt_end_3prime = detect_altprime(genes, sp.where([x.is_alt for x in genes])[0])
                        ### handle 5 prime events
                        for k in range(len(idx_alt_end_5prime)):
                            gene = genes[idx_alt_end_5prime[k]]

                            ### perform liftover between strains if necessary
                            exons = gene.splicegraph.vertices
                            if not 'reference_strain' in CFG:
                                exons_col = exons
                                exons_col_pos = exons
                            else:
                                exons_col = convert_strain_pos_intervals(gene.chr_num, gene.splicegraph.vertices.T, strain, CFG['reference_strain']).T
                                exons_col_pos = convert_strain_pos(gene.chr_num, gene.splicegraph.vertices.T, strain, CFG['reference_strain']).T
                            if exons_col.shape != exons_col_pos.shape: 
                                print 'skipping non-mappable alt 5-prime event'
                                continue
                            
                            for k1 in range(len(exon_alt_end_5prime[k]['fiveprimesites']) - 1):
                                for k2 in range(k1 + 1, len(exon_alt_end_5prime[k]['fiveprimesites'])):

                                    exon_alt1_col = exons_col[:, exon_alt_end_5prime[k]['fiveprimesites'][k1]].T
                                    exon_alt2_col = exons_col[:, exon_alt_end_5prime[k]['fiveprimesites'][k2]].T

                                    ### check if exons overlap
                                    if (exon_alt1_col[0] > exon_alt2_col[1]) or (exon_alt1_col[1] < exon_alt2_col[0]):
                                        continue

                                    event = Event('alt_5prime', gene.chr, gene.chr_num, gene.strand)
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
                                    alt_end_5prime_pos[ridx, i].append(event)

                        ### handle 3 prime events
                        for k in range(len(idx_alt_end_3prime)):
                            gene = genes[idx_alt_end_3prime[k]]

                            ### perform liftover between strains if necessary
                            exons = gene.splicegraph.vertices
                            if not 'reference_strain' in CFG:
                                exons_col = exons
                                exons_col_pos = exons
                            else:
                                exons_col = convert_strain_pos_intervals(gene.chr_num, gene.splicegraph.vertices.T, strain, CFG['reference_strain']).T
                                exons_col_pos = convert_strain_pos(gene.chr_num, gene.splicegraph.vertices.T, strain, CFG['reference_strain']).T
                            if exons_col.shape != exons_col_pos.shape: 
                                print 'skipping non-mappable alt 3-prime event'
                                continue

                            for k1 in range(len(exon_alt_end_3prime[k]['threeprimesites']) - 1):
                                for k2 in range(k1 + 1, len(exon_alt_end_3prime[k]['threeprimesites'])):

                                    exon_alt1_col = exons_col[:, exon_alt_end_3prime[k]['threeprimesites'][k1]].T
                                    exon_alt2_col = exons_col[:, exon_alt_end_3prime[k]['threeprimesites'][k2]].T

                                    ### check if exons overlap
                                    if (exon_alt1_col[0] > exon_alt2_col[1]) or (exon_alt1_col[1] < exon_alt2_col[0]):
                                        continue

                                    event = Event('alt_3prime', gene.chr, gene.chr_num, gene.strand)
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
                                    alt_end_3prime_pos[ridx, i].append(event)
                    else:
                        print '%s and %s already exists' % (fn_out_a5, fn_out_a3)

                ### detect multiple_exon_skips from splicegraph
                if do_mult_exon_skip:
                    if not os.path.exists(fn_out_mes):
                        idx_mult_exon_skip, exon_mult_exon_skip, id_mult_exon_skip = detect_multipleskips(genes, sp.where([x.is_alt for x in genes])[0])
                        if len(id_mult_exon_skip) > 0:
                            for k in id_mult_exon_skip:
                                #k_ = sp.where(sp.array(id_mult_exon_skip) == k)[0]
                                #gene = genes[idx_mult_exon_skip[k_[0]]]
                                gene = genes[idx_mult_exon_skip[k]]

                                ### perform liftover between strains if necessary
                                exons = gene.splicegraph.vertices
                                if not 'reference_strain' in CFG:
                                    exons_col = exons
                                    exons_col_pos = exons
                                else:
                                    exons_col = convert_strain_pos_intervals(gene.chr_num, gene.splicegraph.vertices.T, strain, CFG['reference_strain']).T
                                    exons_col_pos = convert_strain_pos(gene.chr_num, gene.splicegraph.vertices.T, strain, CFG['reference_strain']).T
                                if exons_col.shape != exons_col_pos.shape: 
                                    print 'skipping non-mappable multiple exon skip event'
                                    continue

                                ### build multiple exon skip data structure
                                event = Event('mult_exon_skip', gene.chr, gene.chr_num, gene.strand)
                                event.strain = sp.array([strain])
                                event.exons1 = sp.c_[exons[:, exon_mult_exon_skip[k][0]], exons[:, exon_mult_exon_skip[k][2]]].T
                                event.exons2 = sp.c_[exons[:, exon_mult_exon_skip[k][0]], exons[:, exon_mult_exon_skip[k][1]], exons[:, exon_mult_exon_skip[k][2]]].T
                                event.exons1_col = sp.c_[exons_col[:, exon_mult_exon_skip[k][0]], exons_col[:, exon_mult_exon_skip[k][2]]].T
                                event.exons2_col = sp.c_[exons_col[:, exon_mult_exon_skip[k][0]], exons_col[:, exon_mult_exon_skip[k][1]], exons_col[:, exon_mult_exon_skip[k][2]]].T
                                event.gene_name = sp.array([gene.name])
                                event.gene_idx = idx_mult_exon_skip[k]
                                #event.transcript_type = sp.array([gene.transcript_type])
                                mult_exon_skip_pos[ridx, i].append(event)
                    else:
                        print '%s already exists' % fn_out_mes
            ### genes file does not exist
            else:
                print 'result file not found: %s' % genes_fnames

    ### combine events for all samples
    for ridx in CFG['replicate_idxs']:

        ################################################%
        ### COMBINE INTRON RETENTIONS
        ################################################%
        if do_intron_retention:
            if not os.path.exists(fn_out_ir):
                intron_reten_pos_all = sp.array([item for sublist in intron_reten_pos[ridx, :] for item in sublist])

                ### post process event structure by sorting and making events unique
                events_all = post_process_event_struct(intron_reten_pos_all)

                ### store intron retentions
                print 'saving intron retentions to %s' % fn_out_ir
                cPickle.dump(events_all, open(fn_out_ir, 'w'), -1)
            else:
                print '%s already exists' % fn_out_ir
        
        ################################################%
        ### COMBINE EXON SKIPS
        ################################################%
        if do_exon_skip:
            if not os.path.exists(fn_out_es):
                exon_skip_pos_all = sp.array([item for sublist in exon_skip_pos[ridx, :] for item in sublist])

                ### post process event structure by sorting and making events unique
                events_all = post_process_event_struct(exon_skip_pos_all)

                ### store exon skip events
                print 'saving exon skips to %s' % fn_out_es
                cPickle.dump(events_all, open(fn_out_es, 'w'), -1)
            else:
                print '%s already exists' % fn_out_es

        ################################################%
        ### COMBINE MULTIPLE EXON SKIPS
        ################################################%
        if do_mult_exon_skip:
            if not os.path.exists(fn_out_mes):
                mult_exon_skip_pos_all = sp.array([item for sublist in mult_exon_skip_pos[ridx, :] for item in sublist])

                ### post process event structure by sorting and making events unique
                events_all = post_process_event_struct(mult_exon_skip_pos_all)

                ### store multiple exon skip events
                print 'saving multiple exon skips to %s' % fn_out_mes
                cPickle.dump(events_all, open(fn_out_mes, 'w'), -1)
            else:
                print '%s already exists' % fn_out_mes

        ################################################%
        ### COMBINE ALT FIVE PRIME EVENTS
        ################################################%
        if do_alt_5prime:
            if not os.path.exists(fn_out_a5):
                alt_end_5prime_pos_all = sp.array([item for sublist in alt_end_5prime_pos[ridx, :] for item in sublist])
              
                ### post process event structure by sorting and making events unique
                events_all = post_process_event_struct(alt_end_5prime_pos_all)

                ### curate alt prime events
                ### cut to min len, if alt exon lengths differ
                ### remove, if alt exons do not overlap
                if CFG['curate_alt_prime_events']:
                    events_all = curate_alt_prime(events_all, CFG)

                ### store alt 5 prime events
                print 'saving alt 5 prime events to %s' % fn_out_a5
                cPickle.dump(events_all, open(fn_out_a5, 'w'), -1)
            else:
                print '%s already exists' % fn_out_a5

        ################################################%
        ### COMBINE ALT THREE PRIME EVENTS
        ################################################%
        if do_alt_3prime:
            if not os.path.exists(fn_out_a3):
                alt_end_3prime_pos_all = sp.array([item for sublist in alt_end_3prime_pos[ridx, :] for item in sublist])
               
                ### post process event structure by sorting and making events unique
                events_all = post_process_event_struct(alt_end_3prime_pos_all)

                ### curate alt prime events
                ### cut to min len, if alt exon lengths differ
                ### remove, if alt exons do not overlap
                if CFG['curate_alt_prime_events']:
                    events_all = curate_alt_prime(events_all, CFG)

                ### store alt 3 prime events
                print 'saving alt 3 prime events to %s' % fn_out_a3
                cPickle.dump(events_all, open(fn_out_a3, 'w'), -1)
            else:
                print '%s already exists' % fn_out_a3
