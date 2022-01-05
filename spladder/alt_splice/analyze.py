import sys
import os
import numpy as np
import pickle
import h5py

if __name__ == "__main__":
    __package__ = "modules.alt_splice"

### local imports
from .verify import *
from .write import *
from ..helpers import compute_psi, codeUTF8

def _prepare_count_hdf5(options, OUT, event_features, sample_idx=None):
    
    ### load gene info
    if hasattr(options, 'spladderfile') and os.path.exists(options.spladderfile):
        (genes, inserted) = pickle.load(open(options.spladderfile), 'rb')
    else:
        validate_tag = ''
        if options.validate_sg:
            validate_tag = '.validated'
        if not sample_idx is None:
            (genes, inserted) = pickle.load(open('%s/spladder/genes_graph_conf%i.%s%s.pickle' % (options.outdir, options.confidence, options.samples[sample_idx], validate_tag), 'rb'))
        else:
            (genes, inserted) = pickle.load(open('%s/spladder/genes_graph_conf%i.%s%s.pickle' % (options.outdir, options.confidence, options.merge, validate_tag), 'rb'))

    ### write sample and gene indices to hdf5
    OUT.create_dataset(name='samples', data=codeUTF8(options.samples))
    OUT['strains'] = h5py.SoftLink('/samples')
    OUT.create_dataset(name='event_features', data=codeUTF8(np.array(event_features, dtype='str')))
    OUT.create_dataset(name='gene_names', data=codeUTF8(np.array([x.name for x in genes], dtype='str')))
    OUT.create_dataset(name='gene_chr', data=codeUTF8(np.array([x.chr for x in genes], dtype='str')))
    OUT.create_dataset(name='gene_strand', data=codeUTF8(np.array([x.strand for x in genes], dtype='str')))
    OUT.create_dataset(name='gene_pos', data=np.array([[x.start, x.stop] for x in genes], dtype='int'))


def analyze_events(event_type, bam_fnames, options, sample_idx=None):

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

    event_features = {'mult_exon_skip': ['valid', 'e1_cov', 'e2_cov', 'e3_cov', 'e1e2_conf', 'e2e3_conf', 'e1e3_conf', 'sum_e2_conf', 'num_e2', 'len_e2'],
                      'intron_retention': ['valid', 'e1_cov', 'e2_cov', 'e3_cov', 'e1e3_conf', 'e2_cov_region'],
                      'exon_skip': ['valid', 'e1_cov', 'e2_cov', 'e3_cov', 'e1e2_conf', 'e2e3_conf', 'e1e3_conf'],
                      'mutex_exons': ['valid', 'e1_cov', 'e2_cov', 'e3_cov', 'e4_cov', 'e1e2_conf', 'e1e3_conf', 'e2e4_conf', 'e3e4_conf'],
                      'alt_3prime': ['valid', 'e1_cov', 'e2_cov', 'e3_cov', 'e1e3_conf', 'e2_conf'],
                      'alt_5prime': ['valid', 'e1_cov', 'e2_cov', 'e3_cov', 'e1e3_conf', 'e2_conf']}

    ### check, if confirmed version exists
    if not os.path.exists(fn_out_count):

        events_all = pickle.load(open(fn_out, 'rb'))
        events_all_samples = options.samples

        ### handle case where we did not find any event of this type
        if np.sum([x.event_type == event_type for x in events_all]) == 0:
            OUT = h5py.File(fn_out_count, 'w')
            OUT.create_dataset(name='event_counts', data=[0])
            _prepare_count_hdf5(options, OUT, event_features[event_type], sample_idx=sample_idx)
            OUT.close()
            confirmed_idx = np.array([], dtype='int')
        else:
            if options.merge == 'single':
                (events_all, counts, verified) = verify_all_events(events_all, sample_idx, bam_fnames, event_type, options)
            else:
                (events_all, counts, verified) = verify_all_events(events_all, np.arange(len(options.samples)), bam_fnames, event_type, options)

            psi = np.empty((counts.shape[0], counts.shape[2]), dtype='float')
            iso1 = np.empty((counts.shape[0], counts.shape[2]), dtype='int32')
            iso2 = np.empty((counts.shape[0], counts.shape[2]), dtype='int32')
            for i in range(counts.shape[2]):
                (psi[:, i], iso2[:, i], iso1[:, i])  = compute_psi(counts[:, :, i], event_type, options)

            OUT = h5py.File(fn_out_count, 'w')
            OUT.create_dataset(name='event_counts', data=counts, compression='gzip')
            OUT.create_dataset(name='psi', data=psi, compression='gzip')
            OUT.create_dataset(name='iso1', data=iso1, compression='gzip')
            OUT.create_dataset(name='iso2', data=iso2, compression='gzip')
            OUT.create_dataset(name='gene_idx', data=np.array([x.gene_idx for x in events_all], dtype='int'), compression='gzip')
            OUT.create_dataset(name='verified', data=verified, compression='gzip')
            _prepare_count_hdf5(options, OUT, event_features[event_type], sample_idx=sample_idx)
           
            ### write more event infos to hdf5
            if event_type == 'exon_skip':
                event_pos = np.array([x.exons2.ravel() for x in events_all])
            elif event_type == 'intron_retention':
                event_pos = np.array([x.exons1.ravel() for x in events_all])
            elif event_type in ['alt_3prime', 'alt_5prime']:
                event_pos = np.array([unique_rows(np.c_[x.exons1, x.exons2]).ravel() for x in events_all])
            elif event_type == 'mult_exon_skip':
                event_pos = np.array([x.exons2[[0, 1, -2, -1], :].ravel() for x in events_all])
            elif event_type == 'mutex_exons':
                event_pos = np.array([np.c_[x.exons1[0, :], x.exons1[1, :], x.exons2[1, :], x.exons2[2, :]] for x in events_all])

            OUT.create_dataset(name='event_pos', data=event_pos)

            num_verified = np.sum(verified, axis=0)
            confirmed = num_verified.min(axis=0)
            OUT.create_dataset(name='num_verified', data=num_verified)
            OUT.create_dataset(name='confirmed', data=confirmed)

            confirmed_idx = np.where(confirmed >= 1)[0]
            if confirmed_idx.shape[0] > 0:
                OUT.create_dataset(name='conf_idx', data=confirmed_idx)

            ### close HDF5
            OUT.close()

        ### save events
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
            write_events_txt(fn_out_txt, options.samples, events_all, fn_out_count)

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
        sample_idx = np.arange(options.samples.shape[0])

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
            write_events_txt(fn_out_conf_txt, options.samples[sample_idx], events_all, fn_out_count, event_idx=confirmed_idx)

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
            write_events_tcga(fn_out_conf_tcga, options.samples[sample_idx], events_all, fn_out_count, event_idx=confirmed_idx)

    if options.output_confirmed_icgc:
        if os.path.exists(fn_out_conf_icgc):
            print('%s already exists' % fn_out_conf_icgc)
        else:
            write_events_icgc(fn_out_conf_icgc, options.samples[sample_idx], events_all, fn_out_count, event_idx=confirmed_idx)

    if options.output_filtered_txt:
        fn_out_conf_txt = fn_out_conf.replace('.pickle', '.filt0.05.txt')
        if os.path.exists(fn_out_conf_txt):
            print('%s already exists' % fn_out_conf_txt)
        else:
            print('\nWriting filtered events (sample freq 0.05):')
            cf_idx = np.where([x.confirmed for x in events_all[confirmed_idx]] >= (0.05 * options.samples.shape[0]))[0]
            write_events_txt(fn_out_conf_txt, options.samples[sample_idx], events_all, fn_out_count, event_idx=confirmed_idx[cf_idx])

        fn_out_conf_txt = fn_out_conf.replace('.pickle', '.filt0.1.txt')
        if os.path.exists(fn_out_conf_txt):
            print('%s already exists' %  fn_out_conf_txt)
        else:
            print('\nWriting filtered events (sample freq 0.01):')
            cf_idx = np.where([x.confirmed for x in events_all[confirmed_idx]] >= (0.01 * options.samples.shape[0]))[0]
            write_events_txt(fn_out_conf_txt, options.samples[sample_idx], events_all, fn_out_count, event_idx=confirmed_idx[cf_idx])
