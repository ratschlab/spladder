import os
import sys
import scipy as sp
import math
import re

def default_settings(options):

    ### import necessary paths 
    options.paths = []
    if 'SPLADDER_SRC_PATH' in os.environ:
        options.paths.append(os.environ['SPLADDER_SRC_PATH'])
        options.paths.append('%s/utils' % os.environ['SPLADDER_SRC_PATH'])
        options.paths.append('%s/alt_splice' % os.environ['SPLADDER_SRC_PATH'])
    else:
        ### infer src path from current script path
        SPLADDER_SRC_PATH = os.getcwd()
        options.paths.append(SPLADDER_SRC_PATH)
        options.paths.append('%s/utils' % SPLADDER_SRC_PATH)
        options.paths.append('%s/alt_splice' % SPLADDER_SRC_PATH)

    options.bam_fnames = []
    # TODO
    #if len(CFG['paths']) > 0:
    #    for p in CFG['paths']:
    #        eval('import %s' % p )

    ### settings for adding new intron edges
    options.intron_edges = dict()
    options.intron_edges['min_exon_len'] = 50
    options.intron_edges['min_exon_len_remove'] = 8
    options.intron_edges['vicinity_region'] = 40
    options.intron_edges['insert_intron_retention'] = 1
    options.intron_edges['gene_merges'] = 0
    options.intron_edges['append_new_terminal_exons'] = 1
    options.intron_edges['append_new_terminal_exons_len'] = 200

    options.intron_window = 5000

    ### settings for cassette exons
    options.cassette_exons = dict()
    options.cassette_exons['min_cassette_cov'] = 5
    options.cassette_exons['min_cassette_region'] = 0.9
    options.cassette_exons['min_cassette_rel_diff'] = 0.5

    ### settings for short exon removal
    options.remove_exons = dict()
    options.remove_exons['terminal_short_extend'] = 40
    options.remove_exons['terminal_short_len'] = 10
    options.remove_exons['min_exon_len'] = 50
    options.remove_exons['min_exon_len_remove'] = 10

    options.verify_alt_events = True

    ### settings for verifying exon skips
    options.exon_skip = dict()
    options.exon_skip['min_non_skip_count'] = 3
    options.exon_skip['min_skip_count'] = 3
    options.exon_skip['min_skip_rel_cov'] = 0.05
   #options.exon_skip['max_exon_fold_diff'] = 4
   #options.exon_skip['max_skip_rel_cov'] = 1.5
    options.exon_skip['intron_tolerance'] = 0

    ### settings for verifying mutually exclusive exons
    options.mutex_exons = dict()
    options.mutex_exons['min_skip_rel_cov'] = 0.05
    options.mutex_exons['min_conf_count'] = 2

    ### settings for verifying multiple exon skips
    options.mult_exon_skip = dict()
    options.mult_exon_skip['min_non_skip_count'] = 3
    options.mult_exon_skip['min_skip_count'] = 3
    options.mult_exon_skip['min_skip_rel_cov'] = 0.05
   #options.mult_exon_skip['max_exon_fold_diff'] = 4
   #options.mult_exon_skip['max_skip_rel_cov'] = 1.5

    ### settings for verifying alt prime events
    options.alt_prime = dict()
    options.alt_prime['min_diff_rel_cov'] = 0.05
    options.alt_prime['min_intron_count'] = 3

    ### settings for verifying intron retention events
    options.intron_retention = dict()
    options.intron_retention['min_retention_cov'] = 3
    options.intron_retention['min_retention_region'] = 0.75
    options.intron_retention['min_retention_rel_cov'] = 0.05
    options.intron_retention['min_non_retention_count'] = 3
   #options.intron_retention['max_retention_rel_cov'] = 1.5
    options.intron_retention['min_retention_max_exon_fold_diff']  = 4

    options.sg_min_edge_count = 10
    options.no_reset_conf = False

    options.do_prune = False
    options.do_gen_isoforms = False
    options.do_merge_all = False

    ### define which output files are written
    options.output_filtered_txt = False

    ### settings for truncation detection mode
    options.detect_trunc = False
    options.min_truncation_cov = 5

    ### min number of reads to compute PSI
    options.psi_min_reads = 10

    ### treat introns as unstranded
    options.introns_unstranded = False

    return options


def parse_args(options, identity='main'):

    ### load all default settings
    options = default_settings(options)

    ref_tag = ''
    
    options.event_types = options.event_types.strip(',').split(',')

    if not os.path.exists(options.outdir):
        print('WARNING: Output directory %s does not exist - will be created\n\n' % options.outdir, file=sys.stderr)
        try:
            os.makedirs(options.outdir)
        except OSError:
            print('ERROR: Output directory %s can not be created.\n\n' % options.outdir, file=sys.stderr)
            sys.exit(2)

    ### options specific for main program
    if identity == 'main':

        ### open log file, if specified
        if options.logfile != '-':
            options.log_fname = options.logfile
        else:
            options.log_fname = 'stdout'

        options.bam_fnames = options.bams.strip(',').split(',')

        if not os.path.isfile(options.annotation):
            print('ERROR: Annotation file %s can not be found\n\n' % options.annotation, file=sys.stderr)
            sys.exit(2)
        
        #if options.refstrain != '-':
        #    CFG['reference_strain'] = options.refstrain
        #    ref_tag = '%s:' % options.refstrain

        ### pyproc options
        if options.pyproc:
            options.options_rproc = dict()
            options.options_rproc['mem_req_resubmit']  = [30000, 60000, 80000]
            options.options_rproc['time_req_resubmit'] = [60*60, 80*60, 90*60]
            options.options_rproc['resubmit'] = 3
            options.options_rproc['priority'] = 100
            options.options_rproc['addpaths'] = options.paths
            if options.environment:
                options.options_rproc['environment'] = options.environment

    if identity == 'viz':

        if options.bams != '-':
            options.bam_fnames = options.bams.strip(':').split(':')
            for g, group in enumerate(options.bam_fnames):
                options.bam_fnames[g] = group.strip(',').split(',')
                ### check existence of files
                for fname in options.bam_fnames[g]:
                    if not os.path.isfile(fname):
                        print('ERROR: Input file %s can not be found\n\n' % fname, file=sys.stderr)
                        sys.exit(2)

    if identity == 'test':
        if options.conditionA == '-':
            print('ERROR: At least one sample for condition A required', file=sys.stderr)
            sys.exit(1)
        if options.conditionB == '-':
            print('ERROR: At least one sample for condition B required', file=sys.stderr)
            sys.exit(1)

        options.conditionA = options.conditionA.strip(',').split(',')
        options.conditionB = options.conditionB.strip(',').split(',')
        if len(options.conditionA) > 0 and options.conditionA[0].lower().endswith('txt'):
            options.conditionA = [str(x) for x in sp.loadtxt(options.conditionA[0], dtype='str')]
        if len(options.conditionB) > 0 and options.conditionB[0].lower().endswith('txt'):
            options.conditionB = [str(x) for x in sp.loadtxt(options.conditionB[0], dtype='str')]
        options.conditionA = [re.sub(r'(.[bB][aA][mM]|.[hH][dD][fF]5)|.[nN][pP][zZ]$', '', x) for x in options.conditionA]
        options.conditionB = [re.sub(r'(.[bB][aA][mM]|.[hH][dD][fF]5)|.[nN][pP][zZ]$', '', x) for x in options.conditionB]
        options.conditionA = [os.path.basename(x) for x in options.conditionA]
        options.conditionB = [os.path.basename(x) for x in options.conditionB]

    ### check if we got a list of bam files in a text file instead of a comma separated list
    if len(options.bam_fnames) > 0:
        if identity == 'main' and options.bam_fnames[0].split('.')[-1] == 'txt':
            options.bam_fnames = [str(x) for x in sp.atleast_1d(sp.loadtxt(options.bam_fnames[0], dtype='str'))]
        elif identity == 'viz':
            for g, group in enumerate(options.bam_fnames):
                if group[0].split('.')[-1] == 'txt':
                    options.bam_fnames[g] = [str(x) for x in sp.atleast_1d(sp.loadtxt(group[0], dtype='str'))]

        ### check existence of alignment files
        for fname in options.bam_fnames:
            if not os.path.isfile(fname):
                print('ERROR: Input file %s can not be found\n\n' % fname, file=sys.stderr)
                sys.exit(2)
            if not os.path.isfile(fname + '.bai'):
                print('ERROR: Input file %s is not indexed. %s.bai can not be found\n\n' % (fname, fname), file=sys.stderr)
                sys.exit(2)

    ### assemble strain list
    options.samples = []
    options.strains = []
    if identity in ['viz']:
        for g, group in enumerate(options.bam_fnames):
            options.strains.append([]) 
            options.samples.append([])
            for i in range(len(group)):
                options.samples[-1].append(re.sub(r'(.[bB][aA][mM]|.[hH][dD][fF]5)$', '', group[i].split('/')[-1]))
                options.strains[-1].append('%s%s' % (ref_tag, options.samples[-1][-1]))
            options.strains[-1] = sp.array(options.strains[-1])
        if options.labels != '':
            options.labels = options.labels.split(',')
            assert (len(options.labels) == len(options.bam_fnames)), "Number of labels (%i given) and file names (%i given) needs to match!" % (len(options.labels), len(options.bam_fnames))
    else:
        for i in range(len(options.bam_fnames)):
            options.samples.append(re.sub(r'(.[bB][aA][mM]|.[hH][dD][fF]5)$', '', options.bam_fnames[i].split('/')[-1]))
            options.strains.append('%s%s' % (ref_tag, options.samples[-1]))
        options.strains = sp.array(options.strains)

    ### adapt graph validation requirement to max number of samples
    options.sg_min_edge_count = min(options.sg_min_edge_count, len(options.samples))

    return options


def set_confidence_level(options):

    ### settings for accepted introns
    options.read_filter = dict()
    if options.confidence == 0:
        options.read_filter['intron'] = 350000 
        options.read_filter['exon_len'] = math.ceil(options.readlen * 0.10)
        options.read_filter['mismatch'] = max(2, math.floor(options.readlen * 0.03)) 
        options.read_filter['mincount'] = 1
    elif options.confidence == 1:
        options.read_filter['intron'] = 350000 
        options.read_filter['exon_len'] = math.ceil(options.readlen * 0.15)
        options.read_filter['mismatch'] = max(1, math.floor(options.readlen * 0.02))
        options.read_filter['mincount'] = 2
    elif options.confidence == 2:
        options.read_filter['intron'] = 350000 
        options.read_filter['exon_len'] = math.ceil(options.readlen * 0.20)
        options.read_filter['mismatch'] = max(1, math.floor(options.readlen * 0.01)) 
        options.read_filter['mincount'] = 2
    elif options.confidence == 3:
        options.read_filter['intron'] = 350000 
        options.read_filter['exon_len'] = math.ceil(options.readlen * 0.25)
        options.read_filter['mismatch'] = 0
        options.read_filter['mincount'] = 2

    ### settings for accepted cassette exons
    options.cassette_exon = dict()
    options.cassette_exon['min_cassette_cov'] = 5 
    options.cassette_exon['min_cassette_region'] = 0.9
    options.cassette_exon['min_cassette_rel_diff'] = 0.5 

    ### settings for accepted intron retentions
    if not hasattr(options, 'intron_retention'):
        options.intron_retention = dict()
    if options.confidence == 0:
      options.intron_retention['min_retention_cov'] = 1
      options.intron_retention['min_retention_region'] = 0.75 
      options.intron_retention['min_retention_rel_cov'] = 0.1
      options.intron_retention['max_retention_rel_cov'] = 2 
      options.intron_retention['min_retention_max_exon_fold_diff'] = 4
    elif options.confidence == 1:
      options.intron_retention['min_retention_cov'] = 2
      options.intron_retention['min_retention_region'] = 0.75
      options.intron_retention['min_retention_rel_cov'] = 0.1
      options.intron_retention['max_retention_rel_cov'] = 1.2
      options.intron_retention['min_retention_max_exon_fold_diff'] = 4 
    elif options.confidence == 2:
      options.intron_retention['min_retention_cov'] = 5 
      options.intron_retention['min_retention_region'] = 0.9 
      options.intron_retention['min_retention_rel_cov'] = 0.2 
      options.intron_retention['max_retention_rel_cov'] = 1.2 
      options.intron_retention['min_retention_max_exon_fold_diff'] = 4 
    elif options.confidence == 3:
      options.intron_retention['min_retention_cov'] = 10 
      options.intron_retention['min_retention_region'] = 0.9  
      options.intron_retention['min_retention_rel_cov'] = 0.2 
      options.intron_retention['max_retention_rel_cov'] = 1.2 
      options.intron_retention['min_retention_max_exon_fold_diff'] = 4 

    options.intron_retention['read_filter'] = options.read_filter

    return options

