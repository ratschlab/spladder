import sys
import os
import scipy as sp
import pickle
import h5py

if __name__ == "__main__":
    __package__ = "modules.alt_splice"

### local imports
from .verify import *
from .write import *
from ..rproc import rproc, rproc_wait
from ..helpers import compute_psi, codeUTF8

def _prepare_count_hdf5(options, OUT, event_features, sample_idx=None):
    
    ### load gene info
    if hasattr(options, 'spladderfile') and os.path.exists(options.spladderfile):
        (genes, inserted) = pickle.load(open(options.spladderfile), 'rb')
    else:
        prune_tag = ''
        if options.do_prune:
            prune_tag = '_pruned'
        validate_tag = ''
        if options.validate_sg:
            validate_tag = '.validated'
        if not sample_idx is None:
            (genes, inserted) = pickle.load(open('%s/spladder/genes_graph_conf%i.%s%s%s.pickle' % (options.outdir, options.confidence, options.samples[sample_idx], validate_tag, prune_tag), 'rb'))
        else:
            (genes, inserted) = pickle.load(open('%s/spladder/genes_graph_conf%i.%s%s%s.pickle' % (options.outdir, options.confidence, options.merge, validate_tag, prune_tag), 'rb'))

    ### write strain and gene indices to hdf5
    OUT.create_dataset(name='strains', data=codeUTF8(options.strains))
    feat = OUT.create_group(name='event_features')
    for f in event_features:
        feat.create_dataset(name=f, data=codeUTF8(sp.array(event_features[f], dtype='str')))
    OUT.create_dataset(name='gene_names', data=codeUTF8(sp.array([x.name for x in genes], dtype='str')))
    OUT.create_dataset(name='gene_chr', data=codeUTF8(sp.array([x.chr for x in genes], dtype='str')))
    OUT.create_dataset(name='gene_strand', data=codeUTF8(sp.array([x.strand for x in genes], dtype='str')))
    OUT.create_dataset(name='gene_pos', data=sp.array([[x.start, x.stop] for x in genes], dtype='int'))


def analyze_events(options, event_type, sample_idx=None):

    if options.pyproc and not os.path.exists('%s/event_count_chunks' % options.outdir):
        os.makedirs('%s/event_count_chunks' % options.outdir)

    print('analyzing events with confidence %i' % options.confidence)

    if options.merge == 'single':
        fn_out = '%s/%s_%s_C%i.pickle' % (options.outdir, options.samples[sample_idx], event_type, options.confidence)
    else:
        fn_out = '%s/%s_%s_C%i.pickle' % (options.outdir, options.merge, event_type, options.confidence)
    fn_out_conf = fn_out.replace('.pickle', '.confirmed.pickle')
    fn_out_count = fn_out.replace('.pickle', '.counts.hdf5')

    ### define result files
    if options.compress_text:
        gz_suffix = '.gz'
    else:
        gz_suffix = ''
    fn_out_txt = fn_out.replace('.pickle', '.txt' + gz_suffix)
    fn_out_struc = fn_out.replace('.pickle', '.struc.txt' + gz_suffix)
    fn_out_bed = fn_out.replace('.pickle', '.bed')
    fn_out_gff3 = fn_out.replace('.pickle', '.gff3')
    fn_out_conf_txt = fn_out_conf.replace('.pickle', '.txt' + gz_suffix)
    fn_out_conf_struc = fn_out_conf.replace('.pickle', '.struc.txt' + gz_suffix)
    fn_out_conf_bed = fn_out_conf.replace('.pickle', '.bed')
    fn_out_conf_gff3 = fn_out_conf.replace('.pickle', '.gff3')
    fn_out_conf_tcga = fn_out_conf.replace('.pickle', '.tcga.txt' + gz_suffix)
    fn_out_conf_icgc = fn_out_conf.replace('.pickle', '.icgc.txt' + gz_suffix)

    ### check if there is anything to do
    if os.path.exists(fn_out_txt) and os.path.exists(fn_out_conf_txt) and os.path.exists(fn_out_conf_tcga) and os.path.exists(fn_out_conf_icgc) and os.path.exists(fn_out_conf_gff3):
        print('All output files for %s exist.\n' % event_type)
        return

    event_features = {'mult_exon_skip': ['valid', 'exon_pre_cov', 'exons_cov', 'exon_aft_cov', 'exon_pre_exon_conf', 'exon_exon_aft_conf', 'exon_pre_exon_aft_conf', 'sum_inner_exon_conf', 'num_inner_exon', 'len_inner_exon'],
                      'intron_retention': ['valid', 'intron_cov', 'exon1_cov', 'exon2_cov', 'intron_conf', 'intron_cov_region'],
                      'exon_skip': ['valid', 'exon_cov', 'exon_pre_cov', 'exon_aft_cov', 'exon_pre_exon_conf', 'exon_exon_aft_conf', 'exon_pre_exon_aft_conf'],
                      'mutex_exons': ['valid', 'exon_pre_cov', 'exon1_cov', 'exon2_cov', 'exon_aft_cov', 'exon_pre_exon1_conf', 'exon_pre_exon2_conf', 'exon1_exon_aft_conf', 'exon2_exon_aft_conf'],
                      'alt_3prime': ['valid', 'exon_diff_cov', 'exon_const_cov', 'intron1_conf', 'intron2_conf'],
                      'alt_5prime': ['valid', 'exon_diff_cov', 'exon_const_cov', 'intron1_conf', 'intron2_conf']}

    ### check, if confirmed version exists
    if not os.path.exists(fn_out_count):

        events_all = pickle.load(open(fn_out, 'rb'))
        events_all_strains = options.strains

        ### handle case where we did not find any event of this type
        if sp.sum([x.event_type == event_type for x in events_all]) == 0:
            OUT = h5py.File(fn_out_count, 'w')
            OUT.create_dataset(name='event_counts', data=[0])
            _prepare_count_hdf5(options, OUT, event_features, sample_idx=sample_idx)
            OUT.close()
            confirmed_idx = sp.array([], dtype='int')
        else:
            if not options.pyproc:
                if options.merge == 'single':
                    (events_all, counts) = verify_all_events(events_all, sample_idx, options.bam_fnames, event_type, options)
                else:
                    (events_all, counts) = verify_all_events(events_all, sp.arange(len(options.strains)), options.bam_fnames, event_type, options)
                verified = sp.array([x.verified for x in events_all], dtype='bool')
                for ev in events_all:
                    ev.verified = []

                psi = sp.empty((counts.shape[0], counts.shape[2]), dtype='float')
                iso1 = sp.empty((counts.shape[0], counts.shape[2]), dtype='int32')
                iso2 = sp.empty((counts.shape[0], counts.shape[2]), dtype='int32')
                for i in range(counts.shape[2]):
                    (psi[:, i], iso1[:, i], iso2[:, i])  = compute_psi(counts[:, :, i], event_type, options)

                OUT = h5py.File(fn_out_count, 'w')
                OUT.create_dataset(name='event_counts', data=counts, compression='gzip')
                OUT.create_dataset(name='psi', data=psi, compression='gzip')
                OUT.create_dataset(name='iso1', data=iso1, compression='gzip')
                OUT.create_dataset(name='iso2', data=iso2, compression='gzip')
                OUT.create_dataset(name='gene_idx', data=sp.array([x.gene_idx for x in events_all], dtype='int'), compression='gzip')
                OUT.create_dataset(name='verified', data=verified, compression='gzip')
                _prepare_count_hdf5(options, OUT, event_features, sample_idx=sample_idx)
            else:
                jobinfo = []
                PAR = dict()
                chunk_size_events = 5000
                chunk_size_strains = 500
                for i in range(0, events_all.shape[0], chunk_size_events):
                    idx_events = sp.arange(i, min(i + chunk_size_events, events_all.shape[0]))
                    for j in range(0, len(options.strains), chunk_size_strains):
                        idx_strains = sp.arange(j, min(j + chunk_size_strains, len(options.strains)))
                        PAR['ev'] = events_all[idx_events].copy()
                        PAR['strain_idx'] = idx_strains
                        PAR['list_bam'] = options.bam_fnames
                        PAR['out_fn'] = '%s/event_count_chunks/%s_%i_%i_C%i.pickle' % (options.outdir, event_type, i, j, options.confidence)
                        PAR['event_type'] = event_type
                        PAR['options'] = options
                        if os.path.exists(PAR['out_fn']):
                            print('Chunk event %i, strain %i already completed' % (i, j))
                        else:
                            print('Submitting job %i, event chunk %i/%i, strain chunk %i' % (len(jobinfo) + 1, i, events_all.shape[0], j))
                            jobinfo.append(rproc('verify_all_events', PAR, 10000, options.options_rproc, 60 * 12))
                            #verify_all_events(PAR)
                
                rproc_wait(jobinfo, 20, 1.0, 1)
                
                gene_idx_ = []
                verified = []
                collect_ids = []

                print('Collecting results from chunks ...')
                OUT = h5py.File(fn_out_count, 'w')
                for i in range(0, events_all.shape[0], chunk_size_events):
                    idx_events = sp.arange(i, min(i + chunk_size_events, events_all.shape[0]))
                    for j in range(0, len(options.strains), chunk_size_strains):
                        idx_strains = sp.arange(j, min(j + chunk_size_strains, len(options.strains)))
                        print('\r%i (%i), %i (%i)' % (i, events_all.shape[0], j, len(options.strains)))
                        out_fn = '%s/event_count_chunks/%s_%i_%i_C%i.pickle' % (options.outdir, event_type, i, j, options.confidence)
                        if not os.path.exists(out_fn):
                            print('ERROR: not finished %s' % out_fn, file=sys.stderr)
                            sys.exit(1)
                        ev_, counts_ = pickle.load(open(out_fn, 'rb'))
                        if j == 0:
                            ev = ev_
                            counts = counts_
                            verified_ = [x.verified.astype('bool') for x in ev]
                            collect_ids_ = [x.id for x in ev]
                        else:
                            counts = sp.r_[counts, counts_]
                            for jj in range(len(ev_)):
                                verified_[jj] = sp.r_[verified_[jj], ev_[jj].verified]
                            del counts_
                                
                    psi = sp.empty((counts.shape[0], counts.shape[2]), dtype='float')
                    iso1 = sp.empty((counts.shape[0], counts.shape[2]), dtype='int32')
                    iso2 = sp.empty((counts.shape[0], counts.shape[2]), dtype='int32')
                    for j in range(counts.shape[2]):
                        (psi[:, j], iso1[:, j], iso2[:, j]) = compute_psi(counts[:, :, j], event_type, options) 

                    if i == 0:
                        OUT.create_dataset(name='event_counts', data=counts, maxshape=(len(options.strains), len(event_features[event_type]), None), compression='gzip')
                        OUT.create_dataset(name='psi', data=sp.atleast_2d(psi), maxshape=(psi.shape[0], None), compression='gzip')
                        OUT.create_dataset(name='iso1', data=sp.atleast_2d(iso1), maxshape=(iso1.shape[0], None), compression='gzip')
                        OUT.create_dataset(name='iso2', data=sp.atleast_2d(iso2), maxshape=(iso2.shape[0], None), compression='gzip')
                    else:
                        tmp = OUT['event_counts'].shape
                        OUT['event_counts'].resize((tmp[0], tmp[1], tmp[2] + len(ev)))
                        OUT['event_counts'][:, :, tmp[2]:] = counts
                        tmp = OUT['psi'].shape
                        OUT['psi'].resize((tmp[0], tmp[1] + len(ev)))
                        OUT['psi'][:, tmp[1]:] = psi
                        tmp = OUT['iso1'].shape
                        OUT['iso1'].resize((tmp[0], tmp[1] + len(ev)))
                        OUT['iso1'][:, tmp[1]:] = iso1
                        tmp = OUT['iso2'].shape
                        OUT['iso2'].resize((tmp[0], tmp[1] + len(ev)))
                        OUT['iso2'][:, tmp[1]:] = iso2
                    verified.extend(verified_)
                    collect_ids.extend(collect_ids_)
                    gene_idx_ = sp.r_[gene_idx_, [x.gene_idx for x in ev]]
                    del iso1, iso2, psi, counts, ev, ev_

                verified = sp.array(verified, dtype='bool')

                assert(events_all.shape[0] == verified.shape[0])
                assert(sp.all([events_all[e].id for e in range(events_all.shape[0])] == collect_ids))

                OUT.create_dataset(name='verified', data=verified, dtype='bool', compression='gzip')
                OUT.create_dataset(name='gene_idx', data=gene_idx_)
                _prepare_count_hdf5(options, OUT, event_features, sample_idx=sample_idx)
            
            ### write more event infos to hdf5
            if event_type == 'exon_skip':
                event_pos = sp.array([x.exons2.ravel() for x in events_all])
            elif event_type == 'intron_retention':
                event_pos = sp.array([x.exons1.ravel() for x in events_all])
            elif event_type in ['alt_3prime', 'alt_5prime']:
                event_pos = sp.array([unique_rows(sp.c_[x.exons1, x.exons2]).ravel() for x in events_all])
            elif event_type == 'mult_exon_skip':
                event_pos = sp.array([x.exons2[[0, 1, -2, -1], :].ravel() for x in events_all])
            elif event_type == 'mutex_exons':
                event_pos = sp.array([sp.c_[x.exons1[0, :], x.exons1[1, :], x.exons2[1, :], x.exons2[2, :]] for x in events_all])

            OUT.create_dataset(name='event_pos', data=event_pos)

            num_verified = sp.sum(verified, axis=1)
            confirmed = num_verified.min(axis=1)
            OUT.create_dataset(name='num_verified', data=num_verified)
            OUT.create_dataset(name='confirmed', data=confirmed)

            #verified_count = []
            #for min_verified = 1:length(options.strains),
            #    verified_count(min_verified) = sum([events_all.confirmed] >= min_verified) ;
            
            confirmed_idx = sp.where(confirmed >= 1)[0]
            if confirmed_idx.shape[0] > 0:
                OUT.create_dataset(name='conf_idx', data=confirmed_idx)

            ### close HDF5
            OUT.close()

        ### save events
        #cPickle.dump((events_all_info, events_all_strains), open(fn_out_info, 'w'), -1)
        pickle.dump(confirmed_idx, open(fn_out_conf, 'wb'), -1)

    else:
        print('\nLoading event data from %s' % fn_out)
        events_all = pickle.load(open(fn_out, 'rb'))
        confirmed_idx = pickle.load(open(fn_out_conf, 'rb'))

    if events_all.shape[0] == 0:
        print('\nNo %s event could be found. - Nothing to report' % event_type)
        return
    else:
        print('\nReporting complete %s events:' % event_type)

    if options.output_txt:
        if os.path.exists(fn_out_txt):
            print('%s already exists' % fn_out_txt)
        else:
            write_events_txt(fn_out_txt, options.strains, events_all, fn_out_count)

    if options.output_struc:
        if os.path.exists(fn_out_struc):
            print('%s already exists' % fn_out_struc)
        else:
            write_events_structured(fn_out_struc, events_all, fn_out_count)

    if confirmed_idx.shape[0] == 0:
        print('\nNo %s event could be confirmed. - Nothing to report.' % event_type)
        return 
    else:
        print('\nReporting confirmed %s events:' % event_type)

    if isinstance(sample_idx, int):
        sample_idx = [sample_idx]
    elif sample_idx is None:
        sample_idx = sp.arange(options.strains.shape[0])

    if options.output_gff3:
        if os.path.exists(fn_out_gff3):
            print('%s already exists' % fn_out_gff3)
        else:
            write_events_gff3(fn_out_gff3, events_all)

    if options.output_confirmed_gff3:
        if os.path.exists(fn_out_conf_gff3):
            print('%s already exists' % fn_out_conf_gff3)
        else:
            write_events_gff3(fn_out_conf_gff3, events_all, confirmed_idx)

    if options.output_confirmed_txt:
        if os.path.exists(fn_out_conf_txt):
            print('%s already exists' % fn_out_conf_txt)
        else:
            write_events_txt(fn_out_conf_txt, options.strains[sample_idx], events_all, fn_out_count, event_idx=confirmed_idx)

    if options.output_confirmed_bed:
        if os.path.exists(fn_out_conf_bed):
            print('%s already exists' % fn_out_conf_bed)
        else:
            write_events_bed(fn_out_conf_bed, events_all, idx=confirmed_idx)

    if options.output_bed:
        if os.path.exists(fn_out_bed):
            print('%s already exists' % fn_out_bed)
        else:
            write_events_bed(fn_out_bed, events_all)

    if options.output_confirmed_struc:
        if os.path.exists(fn_out_conf_struc):
            print('%s already exists' % fn_out_conf_struc)
        else:
            write_events_structured(fn_out_conf_struc, events_all, fn_out_count, confirmed_idx)

    if options.output_confirmed_tcga:
        if os.path.exists(fn_out_conf_tcga):
            print('%s already exists' % fn_out_conf_tcga)
        else:
            write_events_tcga(fn_out_conf_tcga, options.strains[sample_idx], events_all, fn_out_count, event_idx=confirmed_idx)

    if options.output_confirmed_icgc:
        if os.path.exists(fn_out_conf_icgc):
            print('%s already exists' % fn_out_conf_icgc)
        else:
            write_events_icgc(fn_out_conf_icgc, options.strains[sample_idx], events_all, fn_out_count, event_idx=confirmed_idx)

    if options.output_filtered_txt:
        fn_out_conf_txt = fn_out_conf.replace('.pickle', '.filt0.05.txt')
        if os.path.exists(fn_out_conf_txt):
            print('%s already exists' % fn_out_conf_txt)
        else:
            print('\nWriting filtered events (sample freq 0.05):')
            cf_idx = sp.where([x.confirmed for x in events_all[confirmed_idx]] >= (0.05 * options.strains.shape[0]))[0]
            write_events_txt(fn_out_conf_txt, options.strains[sample_idx], events_all, fn_out_count, event_idx=confirmed_idx[cf_idx])

        fn_out_conf_txt = fn_out_conf.replace('.pickle', '.filt0.1.txt')
        if os.path.exists(fn_out_conf_txt):
            print('%s already exists' %  fn_out_conf_txt)
        else:
            print('\nWriting filtered events (sample freq 0.01):')
            cf_idx = sp.where([x.confirmed for x in events_all[confirmed_idx]] >= (0.01 * options.strains.shape[0]))[0]
            write_events_txt(fn_out_conf_txt, options.strains[sample_idx], events_all, fn_out_count, event_idx=confirmed_idx[cf_idx])
