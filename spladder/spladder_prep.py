import sys
import os
import numpy as np
import pickle
import h5py
import re

from . import settings
from . import init
from .reads import summarize_chr

def prep_annotation(options):

    ### load confidence level settings
    if not options.no_reset_conf:
        options = settings.set_confidence_level(options)

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
    options = init.append_chrms(np.unique(np.array([x.chr for x in genes], dtype='str')), options)
    del genes

    return options

def prep_sparse_bam_filtered(options):

    for bfn in options.bam_fnames:
        if bfn.lower().split('.')[-1] in ['bam', 'cram'] and not os.path.exists(re.sub(r'\.[bB][aA][mM]|\.[cC][rR][aA][mM]$', '', bfn) + '.conf_%i' % options.confidence + '.filt.hdf5'):

            if not hasattr(options, 'chrm_lookup'):
                IN = pysam.Samfile(bfn, 'rb')
                options = append_chrms([x['SN'] for x in parse_header(IN.text)['SQ']], options)
                IN.close()

            OUT = h5py.File(re.sub(r'\.[bB][aA][mM]|\.[cC][rR][aA][mM]$', '', bfn) + '.conf_%i' % options.confidence + '.filt.hdf5', 'w')
            if options.parallel > 1:
                import multiprocessing as mp
                pool = mp.Pool(processes=options.parallel)
                result = [pool.apply_async(summarize_chr, args=(bfn, str(chrm), options,), kwds={'filter':options.read_filter, 'usetmp':True}) for chrm in sorted(options.chrm_lookup)]
                while result:
                    tmp = result.pop(0).get()
                    P_IN = h5py.File(tmp[1], 'r')
                    for k in P_IN.keys():
                        OUT.create_dataset(name=k, data=P_IN[k][:], compression='gzip')
                    P_IN.close()
                    os.remove(tmp[1])
                    del tmp
                pool.close()
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

def prep_sparse_bam_full(options):

    for bfn in options.bam_fnames:
        if bfn.lower().split('.')[-1] in ['bam', 'cram'] and not os.path.exists(re.sub(r'\.[bB][aA][mM]|\.[cC][rR][aA][mM]$', '', bfn) + '.hdf5'):
            #cnts = dict()

            if not hasattr(options, 'chrm_lookup'):
                IN = pysam.Samfile(bfn, 'rb')
                options = append_chrms([x['SN'] for x in parse_header(IN.text)['SQ']], options)
                IN.close()

            OUT = h5py.File(re.sub(r'\.[bB][aA][mM]|\.[cC][rR][aA][mM]$', '', bfn) + '.hdf5', 'w') 
            if options.parallel > 1:
                import multiprocessing as mp
                pool = mp.Pool(processes=options.parallel)
                result = [pool.apply_async(summarize_chr, args=(bfn, str(chrm), options,), kwds={'usetmp':True}) for chrm in sorted(options.chrm_lookup)]
                while result:
                    tmp = result.pop(0).get()
                    P_IN = h5py.File(tmp[1], 'r')
                    for k in P_IN.keys():
                        OUT.create_dataset(name=k, data=P_IN[k][:], compression='gzip')
                    P_IN.close()
                    os.remove(tmp[1])
                    del tmp
                pool.close()
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

def spladder_prep(options):

    ### parse parameters from options object
    options = settings.parse_args(options)

    ### handle annotation prepping
    if options.annotation != '-':
        options = prep_annotation(options)

    ### convert input BAMs to sparse arrays - filtered case
    if options.sparse_bam and options.bams != '-':
        prep_sparse_bam_filtered(options)
        prep_sparse_bam_full(options)


