import sys
import os
import scipy as sp
import cPickle
import h5py

### local imports
from verify import *
from write import *

def _prepare_count_hdf5(CFG, OUT, events, event_features):
    
    ### load gene info
    if 'spladder_infile' in CFG and os.path.exists(CFG['spladder_infile']):
        genes = cPickle.load(open(CFG['spladder_infile']))
    else:
        prune_tag = ''
        if CFG['do_prune']:
            prune_tag = '_pruned'
        validate_tag = ''
        if CFG['validate_splicegraphs']:
            validate_tag = '.validated'
        genes = cPickle.load(open('%s/spladder/genes_graph_conf%i.%s%s%s.pickle' % (CFG['out_dirname'], CFG['confidence_level'], CFG['merge_strategy'], validate_tag, prune_tag)))

    ### write strain and gene indices to hdf5
    OUT.create_dataset(name='strains', data=CFG['strains'])
    feat = OUT.create_group(name='event_features')
    for f in event_features:
        feat.create_dataset(name=f, data=event_features[f])
    OUT.create_dataset(name='gene_names', data=sp.array([x.name for x in genes], dtype='str'))
    OUT.create_dataset(name='gene_chr', data=sp.array([x.chr for x in genes], dtype='str'))
    OUT.create_dataset(name='gene_strand', data=sp.array([x.strand for x in genes], dtype='str'))
    OUT.create_dataset(name='gene_pos', data=sp.array([[x.start, x.stop] for x in genes], dtype='int'))


def analyze_events(CFG, event_type):

    if CFG['rproc']:
        os.makedirs('%s/event_count_chunks' % CFG['out_dirname'])

    for replicate in CFG['replicate_idxs']:
        
        print 'confidence %i / replicate %i' % (CFG['confidence_level'], replicate)

        if len(CFG['replicate_idxs']) > 1:
            rep_tag = '_R%i' % r_idx
        else:
            rep_tag = ''

        fn_out = '%s/%s_%s%s_C%i.pickle' % (CFG['out_dirname'], CFG['merge_strategy'], event_type, rep_tag, CFG['confidence_level'])
        fn_out_conf = fn_out.replace('.pickle', '.confirmed.pickle')
        fn_out_count = fn_out.replace('.pickle', '.counts.hdf5')

        ### define result files
        fn_out_txt = fn_out.replace('.pickle', '.txt')
        fn_out_conf_txt = fn_out_conf.replace('.pickle', '.txt')
        fn_out_conf_tcga = fn_out_conf.replace('.pickle', '.tcga.txt')
        fn_out_conf_gff3 = fn_out_conf.replace('.pickle', '.gff3')

        ### check if there is anything to do
        if os.path.exists(fn_out_txt) and os.path.exists(fn_out_conf_txt) and os.path.exists(fn_out_conf_tcga) and os.path.exists(fn_out_conf_gff3):
            print 'All output files for %s exist.\n' % event_type
            continue

        event_features = {'exon_skip':7, 'alt_3prime':5, 'alt_5prime':5, 'intron_retention':6, 'mult_exon_skip':9}

        ### check, if confirmed version exists
        if not os.path.exists(fn_out_count):

            events_all_ = cPickle.load(open(fn_out, 'r'))
            if isinstance(events_all_, tuple):
                events_all = events_all_[0]
                events_all_strains = events_all_[1]
            else:
                events_all = events_all_
                events_all_strains = None

            ### add strain information, so we can do two way chunking!
            if events_all_strains is None:
                events_all_strains = CFG['strains']
                
            ### handle case where we did not find any event of this type
            if sp.sum([x.event_type == event_type for x in events_all]) == 0:
                OUT = h5py.File(fn_out_count, 'w')
                OUT.create_dataset(name='event_counts', data=[0])
                _prepare_count_hdf5(CFG, OUT, events_all, event_features)
                OUT.close()
                confirmed_idx = sp.array([], dtype='int')
            else:
                if not CFG['rproc']:
                    OUT = h5py.File(fn_out_count, 'w')
                    #events_all = verify_all_events(events_all, range(len(CFG['strains'])), CFG['bam_fnames'][replicate, :], event_type, CFG)
                    # TODO handle replicate setting
                    (events_all, counts) = verify_all_events(events_all, range(len(CFG['strains'])), CFG['bam_fnames'], event_type, CFG)
                    OUT.create_dataset(name='event_counts', data=counts, compression='gzip')
                    _prepare_count_hdf5(CFG, OUT, events_all, event_features)
                else:
                    jobinfo = rproc_empty()
                    chunk_size_events = 1000
                    chunk_size_strains = 500
                    job_nr = 1
                    for i in range(0, events_all.shape[0], chunk_size_events):
                        idx_events = sp.arange(i, min(i + chunk_size_events - 1, events_all.shape[0]))
                        for j in range(0, len(CFG['strains']), chunk_size_strains):
                            idx_strains = sp.arang(j, min(j + chunk_size_strains - 1, len(CFG['strains'])))
                            PAR['ev'] = events_all[idx_events]
                            PAR['strain_idx'] = idx_strains
                            #PAR['list_bam'] = CFG['bam_fnames'][replicate, :]
                            # TODO handle replicate setting
                            PAR['list_bam'] = CFG['bam_fnames']
                            PAR['out_fn'] = '%s/event_count_chunks/%s_%i_%i_R%i_C%i.pickle' % (CFG['out_dirname'], event_type, i, j, replicate, CFG['confidence_level'])
                            PAR['event_type'] = event_type
                            PAR['CFG'] = CFG
                            if os.path.exists(PAR['out_fn']):
                                print 'Chunk event %i, strain %i already completed' % (i, j)
                            else:
                                print 'Submitting job %i, event chunk %i, strain chunk %i' % (job_nr, i, j)
                                jobinfo[job_nr] = rproc('verify_all_events', PAR, 8000, CFG['options_rproc'], 60)
                                job_nr += 1
                    
                    jobinfo, nr_crashed = rproc_wait(jobinfo, 20, 1, 1)
                    
                    events_all_ = []
                    gene_idx_ = []
                    print 'Collecting results from chunks ...'
                    for i in range(0, events_all.shape[0], chunk_size_events):
                        idx_events = sp.arange(i, min(i + chunk_size_events - 1, events_all.shape[0]))
                        ev_ = []
                        for j in range(0, len(CFG['strains']), chunk_size_strains):
                            idx_strains = sp.arange(j, min(j + chunk_size_strains - 1, len(CFG['strains'])))
                            print '\r%i (%i), %i (%i)' % (i, events_all.shape[0], j, len(CFG['strains']))
                            out_fn = '%s/event_count_chunks/%s_%i_%i_R%i_C%i.pickle' % (CFG['out_dirname'], event_type, i, j, replicate, CFG['confidence_level'])
                            if not os.path.exists(out_fn):
                                print >> sys.stderr, 'ERROR: not finished %s' % out_fn
                                sys.exit(1)
                            ev, counts = cPickle.load(open(out_fn, 'r'))
                            if j == 0:
                                ev_ = ev
                            else:
                                for jj in range(len(ev_)):
                                    ev_[jj].verified[idx_strains, :] = ev[jj].verified
                                    
                        if i == 0:
                            OUT.create_dataset(name='event_counts', data=counts, maxshape=(len(CFG['strains']), len(event_features), None))
                        else:
                            tmp = OUT['event_counts'].shape
                            OUT['event_counts'].resize((tmp[0], tmp[1], tmp[2] + len(ev_)))
                            OUT['event_counts'][:, :, tmp[2]:] = counts
                        events_all_ = sp.r_[events_all_, ev_]
                        gene_idx_ = sp.r_[gene_idx_, [x.gene_idx for x in ev_]]

                    assert(events_all.shape[0] == events_all_.shape[0])
                    assert(sp.all([sp.all(events_all[e].exons1 == events_all_[e].exons1) for e in events_all.shape[0]]))
                    OUT.create_dataset(name='gene_idx', data=gene_idx_)
                    events_all = events_all_
                
                ### write more event infos to hdf5
                if event_type == 'exon_skip':
                    event_pos = sp.array([x.exons2.ravel() for x in events_all])
                elif event_type == 'intron_retention':
                    event_pos = sp.array([x.exons2.ravel() for x in events_all])
                elif event_type in ['alt_3prime', 'alt_5prime']:
                    event_pos = sp.array([unique_rows(sp.c_[x.exons1, x.exons2]).ravel() for x in events_all])
                elif event_type == 'mult_exon_skip':
                    event_pos = sp.array([x.exons2[[0, 1, -2, -1], :].ravel() for x in events_all])

                OUT.create_dataset(name='event_pos', data=event_pos)

                for i in range(events_all.shape[0]):
                    events_all[i].num_verified = sp.sum(events_all[i].verified, axis=0)
                    events_all[i].confirmed = sp.array(events_all[i].num_verified).min()
                
                num_verified = sp.array([x.num_verified for x in events_all])

                #verified_count = []
                #for min_verified = 1:length(CFG.strains),
                #    verified_count(min_verified) = sum([events_all.confirmed] >= min_verified) ;
                
                confirmed_idx = sp.where([x.confirmed for x in events_all] >= 1)[0]
                
                OUT.create_dataset(name='conf_idx', data=confirmed_idx)
                OUT.create_dataset(name='verified', data=num_verified)

                ### close HDF5
                OUT.close()

            ### save events
            cPickle.dump((events_all, events_all_strains), open(fn_out, 'w'), -1)
            cPickle.dump(confirmed_idx, open(fn_out_conf, 'w'), -1)

        else:
            (events_all, events_all_strains) = cPickle.load(open(fn_out, 'r'))
            confirmed_idx = cPickle.load(open(fn_out_conf, 'r'))

        if events_all.shape[0] == 0:
            print '\nNo %s event could be found. - Nothing to report' % event_type
            continue
        else:
            print '\nReporting complete %s events:' % event_type

        if os.path.exists(fn_out_txt):
            print '%s already exists' % fn_out_txt
        else:
            write_events_txt(fn_out_txt, CFG['strains'], events_all)

        if events_confirmed.shape[0] == 0:
            print '\nNo %s event could be confirmed. - Nothing to report.' % event_type
            continue
        else:
            print '\nReporting confirmed %s events:' % event_type

        if os.path.exists(fn_out_conf_gff3):
            print '%s already exists' % fn_out_conf_gff3
        else:
            write_events_gff3(fn_out_conf_gff3, events_confirmed)

        if os.path.exists(fn_out_conf_txt):
            print '%s already exists' % fn_out_conf_txt
        else:
            write_events_txt(fn_out_conf_txt, CFG['strains'], events_confirmed)

        if os.path.exists(fn_out_conf_tcga):
            print '%s already exists' % fn_out_conf_tcga
        else:
            write_events_tcga(fn_out_conf_tcga, CFG['strains'], events_confirmed)

        fn_out_conf_txt = fn_out_conf.replace('.pickle', '.filt0.05.txt')
        if os.path.exists(fn_out_conf_txt):
            print '%s already exists' % fn_out_conf_txt
        else:
            print '\nWriting filtered events (sample freq 0.05):'
            cf_idx = sp.where([x.confirmed for x in events_confirmed] >= (0.05 * len(events_confirmed[0].detected)))[0]
            write_events_txt(fn_out_conf_txt, CFG['strains'], events_confirmed[cf_idx])

        fn_out_conf_txt = fn_out_conf.replace('.pickle', '.filt0.1.txt')
        if os.path.exists(fn_out_conf_txt):
            print '%s already exists' %  fn_out_conf_txt
        else:
            print '\nWriting filtered events (sample freq 0.01):'
            cf_idx = sp.where([x.confirmed for x in events_confirmed] >= (0.1 * len(events_confirmed[0].detected)))[0]
            write_events_txt(fn_out_conf_txt, CFG['strains'], events_confirmed[cf_idx])
