def analyze_events(CFG, event_type):

    if CFG['rproc']:
        os.makedirs('%s/event_count_chunks' % CFG['out_dirname'])

    for replicate in CFG['replicate_idxs']:
        
        print, 'confidence %i / replicate %i' % (CFG['confidence_level'], replicate)

        if len(CFG['replicate_idxs') > 1:
            rep_tag = '_R%i' % r_idx
        else:
            rep_tag = ''

        fn_out = '%s/%s_%s%s_C%i.mat' % (CFG['out_dirname'], CFG['merge_strategy'], event_type, rep_tag, CFG['confidence_level'])
        fn_out_conf = fn_out.replace('.mat', '.confirmed.mat')
        fn_out_count = fn_out.replace('.mat', '.counts.mat')
        fn_out_info = fn_out.replace('.mat', '.info.mat')

        ### define result files
        fn_out_txt = fn_out.replace('.mat', '.txt')
        fn_out_conf_txt = fn_out_conf.replace('.mat', '.txt')
        fn_out_conf_tcga = fn_out_conf.replace('.mat', '.tcga.txt')
        fn_out_conf_gff3 = fn_out_conf.replace('.mat', '.gff3')

        ### check if there is anything to do
        if os.path.exists(fn_out_txt) and os.path.exists(fn_out_conf_txt) and os.path.exists(fn_out_conf_tcga) and os.path.exists(fn_out_conf_gff3):
            print 'All output files for %s exist.\n' % event_type
            continue

        ### check, if confirmed version exists
        if not os.path.exists(fn_out_info):

            events_all_ = cPickle.load(open(fn_out, 'r'))
            if isinstance(events_all_, tuple):
                events_all = events_all_[0]
                events_all_strains = events_all_[1]
            else:
                events_all = events_all_

            ### handle case where we did not find any event of this type
            if len([1 for x in events_all if x.event_type == event_type]) == 0:
                cPickle.dump(events_all, open(fn_out_count.replace('.mat', '.chunk1.mat'), 'w'), -1)
                events_confirmed = events_all
                cPickle.dump(events_confirmed, open(fn_out_conf.replace('.mat', '.chunk1.mat'), 'w'), -1)
                max_len_all = 0
                max_len_confirmed = 0
                chunksize = 200
                cPickle.dump((chunksize, max_len_all, max_len_confirmed), open(fn_out_info, 'w'), -1)
            else:
                ### add strain information, so we can do two way chunking!
                if events_all_strains is None:
                    events_all_strains = CFG['strains']
                    cPickle.dump((events_all, events_all_strains), open(fn_out, 'w'), -1)
                
                if not CFG['rproc']:
                    events_all = verify_all_events(events_all, range(len(CFG['strains'])), CFG['bam_fnames'][replicate, :], event_type, CFG)
                else:
                    jobinfo = rproc_empty()
                    chunk_size_events = 25
                    chunk_size_strains = 50
                    job_nr = 1
                    for i in range(0, events_all.shape[0], chunk_size_events):
                        idx_events = sp.arange(i, min(i + chunk_size_events - 1, events_all.shape[0]))
                        for j in range(0, len(CFG['strains']), chunk_size_strains):
                            idx_strains = sp.arang(j, min(j + chunk_size_strains - 1, len(CFG['strains'])))
                            PAR['ev'] = events_all[idx_events]
                            PAR['strain_idx'] = idx_strains
                            PAR['list_bam'] = CFG['bam_fnames'][replicate, :]
                            PAR['out_fn'] = '%s/event_count_chunks/%s_%i_%i_R%i_C%i.mat' % (CFG['out_dirname'], event_type, i, j, replicate, CFG['confidence_level'])
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
                    print 'Collecting results from chunks ...'
                    for i in range(0, events_all.shape[0], chunk_size_events):
                        idx_events = sp.arange(i, min(i + chunk_size_events - 1, events_all.shape[0]))
                        ev_ = []
                        for j in range(0, len(CFG['strains'], chunk_size_strains):
                            idx_strains = sp.arange(j, min(j + chunk_size_strains - 1, len(CFG['strains'])))
                            print '\r%i (%i), %i (%i)', i, events_all.shape[0], j, len(CFG['strains']))
                            out_fn = '%s/event_count_chunks/%s_%i_%i_R%i_C%i.mat' % (CFG['out_dirname'], event_type, i, j, replicate, CFG['confidence_level'])
                            if not os.path.exists(out_fn):
                                print >> sys.stderr, 'ERROR: not finished %s' % out_fn
                                sys.exit(1)
                            ev = cPickle.load(open(out_fn, 'r'))
                            if j == 0:
                                ev_ = ev
                            else:
                                for jj in range(len(ev_)):
                                    ev_[jj].verified[idx_strains, :] = ev[jj].verified
                                    ev_[jj].info[idx_strains] = ev[jj].info
                        events_all_ = sp.r_[events_all_, ev_]
                    assert(events_all.shape[0] == events_all_.shape[0])
                    assert(sp.all([sp.all(events_all[e].exons1 == events_all_[e].exons1) for e in events_all.shape[0]]))
                    events_all = events_all_
                
                for i in range(events_all.shape[0]):
                    events_all[i].num_verified = sp.sum(events_all[i].verified, axis=0)
                    events_all[i].confirmed = sp.array(events_all[i].num_verified).min()
                
                #verified_count = []
                #for min_verified = 1:length(CFG.strains),
                #    verified_count(min_verified) = sum([events_all.confirmed] >= min_verified) ;
                
                min_verified = 1
                v_idx = sp.where([x.confirmed for x in events_all] >= min_verified)[0]
                events_confirmed = events_all[v_idx]
                
                ### save events per chromosome 
                chunksize = 200
                max_len_all = events_all.shape[0]
                max_len_confirmed = events_confirmed.shape[0]

                for c1_idx in range(0, max_len_all, chunksize):
                    fn_out_count_chunk = fn_out_count.replace('.mat', '.chunk%i.mat' % c1_idx)
                    if not os.path.exists(fn_out_count_chunk):
                        tmp_idx = sp.arange(c1_idx, min(c1_idx + chunksize - 1, max_len_all))
                        print 'saving counted %s events to %s' % (event_type, fn_out_count_chunk)
                        cPickle.dump(events_all[tmp_idx], open(fn_out_count_chunk, 'w'), -1)
                    else:
                        print '%s exists' % fn_out_count_chunk

                for c1_idx in range(0, max_len_confirmed, chunksize):
                    fn_out_conf_chunk = fn_out_conf.replace('.mat', '.chunk%i.mat' % c1_idx)
                    if not os.path.exists(fn_out_conf_chunk):
                        tmp_idx = sp.arange(c1_idx, min(c1_idx + chunksize - 1, max_len_confirmed))
                        print 'saving confirmed %s events to %s' % (event_type, fn_out_conf_chunk) 
                        cPickle.dump(events_confirmed[tmp_idx], open(fn_out_conf_chunk, 'w'))
                    else:
                        print '%s exists' % fn_out_conf_chunk

                ### save summary file
                cPickle.dump((chunksize, max_len_all, max_len_confirmed), open(fn_out_info, 'w'), -1)
        else:
            print '%s exists - loading chunk-wise!\n' % fn_out_info
            (chunksize, max_len_all, max_len_confirmed) = cPickle.load(open(fn_out_info, 'r'))
            for c1_idx in range(0, max_len_confirmed, chunksize):
                print '...loading chunk %i' % c1_idx
                fn_out_conf_chunk = fn_out_conf.replace('.mat', '.chunk%i.mat' % c1_idx)
                events_confirmed_ = cPickle.load(open(fn_out_conf_chunk, 'r'))
                if c1_idx == 0:
                    events_confirmed = events_confirmed_.copy()
                else:
                    events_confirmed = sp.r_[events_confirmed, events_confirmed_]

            for c1_idx in range(0, max_len_all, chunksize):
                print '...loading chunk %i' % c1_idx
                fn_out_count_chunk = fn_out_count.replace('.mat', '.chunk%i.mat' % c1_idx)
                events_all_ = cPcikle.load(open(fn_out_count_chunk, 'r'))
                if c1_idx == 0:
                    events_all = events_all_.copy()
                else:
                    events_all = sp.r_[events_all, events_all_]

        if isempty(events_all) || isempty([events_all.event_type]),
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

        fn_out_conf_txt = fn_out_conf.replace('.mat', '.filt0.05.txt')
        if os.path.exists(fn_out_conf_txt):
            print '%s already exists' % fn_out_conf_txt
        else:
            print '\nWriting filtered events (sample freq 0.05):'
            cf_idx = sp.where([x.confirmed for x in events_confirmed] >= (0.05 * len(events_confirmed[0].detected)))[0]
            write_events_txt(fn_out_conf_txt, CFG['strains'], events_confirmed[cf_idx])

        fn_out_conf_txt = fn_out_conf.replace('.mat', '.filt0.1.txt')
        if os.path.exists(fn_out_conf_txt):
            print '%s already exists' %  fn_out_conf_txt
        else:
            print '\nWriting filtered events (sample freq 0.01):'
            cf_idx = sp.where([x.confirmed for x in events_confirmed] >= (0.1 * len(events_confirmed[0].detected)))[0]
            write_events_txt(fn_out_conf_txt, CFG['strains'], events_confirmed[cf_idx])
