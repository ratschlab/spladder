import os
import sys
import scipy as sp
import math
import re

def default_settings():

    CFG = dict()

    ### import necessary paths 
    CFG['paths'] = []
    if 'SPLADDER_SRC_PATH' in os.environ:
        CFG['paths'].append(os.environ['SPLADDER_SRC_PATH'])
        CFG['paths'].append('%s/utils' % os.environ['SPLADDER_SRC_PATH'])
        CFG['paths'].append('%s/alt_splice' % os.environ['SPLADDER_SRC_PATH'])
    else:
        ### infer src path from current script path
        #SPLADDER_SRC_PATH = fileparts(which('default_settings'))
        SPLADDER_SRC_PATH = os.getcwd()
        CFG['paths'].append(SPLADDER_SRC_PATH)
        CFG['paths'].append('%s/utils' % SPLADDER_SRC_PATH)
        CFG['paths'].append('%s/alt_splice' % SPLADDER_SRC_PATH)

    if 'SPLADDER_PATH' in os.environ:
        CFG['paths'].append('%s/mex' % os.environ['SPLADDER_PATH'])
    else:
        #SPLADDER_SRC_PATH = fileparts(which('default_settings'))
        SPLADDER_SRC_PATH = os.getcwd()
        CFG['paths'].append('%s/../mex' % SPLADDER_SRC_PATH)
    
    CFG['bam_fnames'] = []
    # TODO
    #if len(CFG['paths']) > 0:
    #    for p in CFG['paths']:
    #        eval('import %s' % p )

    ### settings for adding new intron edges
    CFG['intron_edges'] = dict()
    CFG['intron_edges']['min_exon_len'] = 50
    CFG['intron_edges']['min_exon_len_remove'] = 8
    CFG['intron_edges']['vicinity_region'] = 40
    CFG['intron_edges']['insert_intron_retention'] = 1
    CFG['intron_edges']['gene_merges'] = 0
    CFG['intron_edges']['append_new_terminal_exons'] = 1
    CFG['intron_edges']['append_new_terminal_exons_len'] = 200

    CFG['intron_window'] = 5000

    ### settings for cassette exons
    CFG['cassette_exons'] = dict()
    CFG['cassette_exons']['min_cassette_cov'] = 5
    CFG['cassette_exons']['min_cassette_region'] = 0.9
    CFG['cassette_exons']['min_cassette_rel_diff'] = 0.5

    ### settings for short exon removal
    CFG['remove_exons'] = dict()
    CFG['remove_exons']['terminal_short_extend'] = 40
    CFG['remove_exons']['terminal_short_len'] = 10
    CFG['remove_exons']['min_exon_len'] = 50
    CFG['remove_exons']['min_exon_len_remove'] = 10

    ### settings for splice graph augmentations
    CFG['do_insert_intron_retentions'] = 1
    CFG['do_insert_cassette_exons'] = 1
    CFG['do_insert_intron_edges'] = 1
    CFG['do_remove_short_exons'] = 0
    CFG['do_infer_splice_graph'] = 0

    CFG['insert_intron_iterations'] = 5
    CFG['merge_strategy'] = 'merge_graphs' ### alternatives are: merge_graphs, single, merge_all
    CFG['confidence_level'] = 3
    CFG['validate_splicegraphs'] = 0
    CFG['same_genestruct_for_all_samples'] = 1
    CFG['curate_alt_prime_events'] = 1
    CFG['replicate_idxs'] = [0]
    CFG['verify_alt_events'] = 1

    ### settings for verifying exon skips
    CFG['exon_skip'] = dict()
    CFG['exon_skip']['min_non_skip_count'] = 3
    CFG['exon_skip']['min_skip_count'] = 3
    CFG['exon_skip']['min_skip_rel_cov'] = 0.05
   #CFG['exon_skip']['max_exon_fold_diff'] = 4
   #CFG['exon_skip']['max_skip_rel_cov'] = 1.5
    CFG['exon_skip']['intron_tolerance'] = 0

    ### settings for verifying mutually exclusive exons
    CFG['mutex_exons'] = dict()
    CFG['mutex_exons']['min_skip_rel_cov'] = 0.05
    CFG['mutex_exons']['min_conf_count'] = 2

    ### settings for verifying multiple exon skips
    CFG['mult_exon_skip'] = dict()
    CFG['mult_exon_skip']['min_non_skip_count'] = 3
    CFG['mult_exon_skip']['min_skip_count'] = 3
    CFG['mult_exon_skip']['min_skip_rel_cov'] = 0.05
   #CFG['mult_exon_skip']['max_exon_fold_diff'] = 4
   #CFG['mult_exon_skip']['max_skip_rel_cov'] = 1.5

    ### settings for verifying alt prime events
    CFG['alt_prime'] = dict()
    CFG['alt_prime']['min_diff_rel_cov'] = 0.05
    CFG['alt_prime']['min_intron_count'] = 3

    ### settings for verifying intron retention events
    CFG['intron_retention'] = dict()
    CFG['intron_retention']['min_retention_cov'] = 3
    CFG['intron_retention']['min_retention_region'] = 0.75
    CFG['intron_retention']['min_retention_rel_cov'] = 0.05
    CFG['intron_retention']['min_non_retention_count'] = 3
   #CFG['intron_retention']['max_retention_rel_cov'] = 1.5
    CFG['intron_retention']['min_retention_max_exon_fold_diff']  = 4

    CFG['count_segment_graph'] = 0

    ### set I/O and verbosity
    CFG['verbose'] = 1
    CFG['debug'] = 0
    CFG['fd_log'] = 1

    CFG['sg_min_edge_count'] = 10
    CFG['no_reset_conf'] = 0

    CFG['do_prune'] = 0
    CFG['do_gen_isoforms'] = 0
    CFG['do_merge_all'] = 0

    CFG['is_half_open'] = 1

    CFG['run_as_analysis'] = 1
    CFG['event_types'] = ['exon_skip', 'intron_retention', 'alt_3prime', 'alt_5prime', 'mult_exon_skip']

    CFG['read_length'] = 36
    CFG['var_aware'] = 0
    CFG['primary_only'] = False

    CFG['rproc'] = 0
    CFG['parallel'] = 1

    CFG['bam_to_sparse'] = 0
    CFG['ignore_mismatch_tag'] = False

    ### define which output files are written
    CFG['output_txt'] = False
    CFG['output_bed'] = False
    CFG['output_struc'] = False
    CFG['output_confirmed_gff3'] = True
    CFG['output_confirmed_txt'] = True
    CFG['output_confirmed_bed'] = False
    CFG['output_confirmed_struc'] = False
    CFG['output_filtered_txt'] = False
    CFG['output_confirmed_tcga'] = False
    CFG['output_confirmed_icgc'] = True
    CFG['compress_text'] = True

    ### settings for truncation detection mode
    CFG['detect_trunc'] = False
    CFG['count_intron_cov'] = False
    CFG['min_truncation_cov'] = 5

    CFG['psi_min_reads'] = 10
    CFG['diagnose_plots'] = False

    CFG['filter_overlapping_genes'] = False
    CFG['filter_overlapping_exons'] = False
    CFG['filter_overlapping_transcripts'] = False

    ### edge limit for detecting events from a graph, if a graph has
    ### more edges, the gene will be ignored for event detection
    CFG['detect_edge_limit'] = 500 #2000

    CFG['introns_unstranded'] = False

    return CFG


def parse_args(options, identity='main'):

    ### load all default settings
    CFG = default_settings()

    ref_tag = ''
    
    ### general options
    if options.verbose in ['n', 'y']:
        CFG['verbose'] = (options.verbose == 'y')
    else:
        print >> sys.stderr, 'ERROR: option verbose should have value y or n, but has %s' % options.verbose
        sys.exit(1)

    if options.debug in ['n', 'y']:
        CFG['debug'] = (options.debug == 'y')
    else:
        print >> sys.stderr, 'ERROR: option debug should have value y or n, but has %s' % options.debug
        sys.exit(1)

    CFG['event_types'] = options.event_types.strip(',').split(',')

    if options.outdir == '-':
        print >> sys.stderr, 'ERROR: please provide the mandatory parameter: out directory\n\n'
        options.parser.print_help()
        sys.exit(2)
    else:
        if not os.path.exists(options.outdir):
            print >> sys.stderr, 'WARNING: Output directory %s does not exist - will be created\n\n' % options.outdir
            try:
                os.makedirs(options.outdir)
            except OSError:
                print >> sys.stderr, 'ERROR: Output directory %s can not be created.\n\n' % options.outdir
                sys.exit(2)
        CFG['out_dirname'] = options.outdir

    ### options specific for main program
    if identity == 'main':
        if options.insert_ir in ['n', 'y']:
            CFG['do_insert_intron_retentions'] = (options.insert_ir == 'y')
        else:
            print >> sys.stderr, 'ERROR: option insert_ir should have value y or n, but has %s' % options.insert_ir
            sys.exit(1)

        if options.insert_es in ['n', 'y']:
            CFG['do_insert_cassette_exons'] = (options.insert_es == 'y')
        else:
            print >> sys.stderr, 'ERROR: option insert_es should have value y or n, but has %s' % options.insert_es
            sys.exit(1)

        if options.insert_ni in ['n', 'y']:
            CFG['do_insert_intron_edges'] = (options.insert_ni == 'y')
        else:
            print >> sys.stderr, 'ERROR: option insert_ni should have value y or n, but has %s' % options.insert_ni
            sys.exit(1)

        if options.remove_se in ['n', 'y']:
            CFG['do_remove_short_exons'] = (options.remove_se == 'y')
        else:
            print >> sys.stderr, 'ERROR: option remove_se should have value y or n, but has %s' % options.remove_se
            sys.exit(1)

        if options.infer_sg in ['n', 'y']:
            CFG['do_infer_splice_graph'] = (options.infer_sg == 'y')
        else:
            print >> sys.stderr, 'ERROR: option infer_sg should have value y or n, but has %s' % options.infer_sg
            sys.exit(1)

        if options.var_aware in ['n', 'y']:
            CFG['var_aware'] = (options.var_aware == 'y')
        else:
            print >> sys.stderr, 'ERROR: option var_aware should have value y or n, but has %s' % options.var_aware
            sys.exit(1)

        if options.primary_only in ['n', 'y']:
            CFG['primary_only'] = (options.primary_only == 'y')
        else:
            print >> sys.stderr, 'ERROR: option primary_only should have value y or n, but has %s' % options.primary_only
            sys.exit(1)

        if options.intron_cov in ['n', 'y']:
            CFG['count_intron_cov'] = (options.intron_cov == 'y')
        else:
            print >> sys.stderr, 'ERROR: option intron_cov should have value y or n, but has %s' % options.intron_cov

        if options.quantify_graph in ['n', 'y']:
            CFG['count_segment_graph'] = (options.quantify_graph == 'y')
        else:
            print >> sys.stderr, 'ERROR: option quantify_graph should have value y or n, but has %s' % options.quantify_graph

        if options.ignore_mismatches in ['n', 'y']:
            CFG['ignore_mismatch_tag'] = (options.ignore_mismatches == 'y')
        else:
            print >> sys.stderr, 'ERROR: option ignore mismatches bam should have value y or n, but has %s' % options.ignore_mismatches
    
        if options.output_struc in ['n', 'y']:
            CFG['output_struc'] = (options.output_struc == 'y')
            CFG['output_confirmed_struc'] = (options.output_struc == 'y')
        else:
            print >> sys.stderr, 'ERROR: option output struc value y or n, but has %s' % options.output_struc

        if options.filter_overlap_genes in ['n', 'y']:
            CFG['filter_overlapping_genes'] = (options.filter_overlap_genes == 'y')
        else:
            print >> sys.stderr, 'ERROR: option filter overlap genes should have value y or n, but has %s' % options.filter_overlap_genes
    
        if options.compress_text in ['n', 'y']:
            CFG['compress_text'] = (options.compress_text == 'y')
        else:
            print >> sys.stderr, 'ERROR: option compress text should have value y or n, but has %s' % options.compress_text
    
        ### option to store sparse BAM representation
        if options.sparse_bam in ['n', 'y']:
            CFG['bam_to_sparse'] = (options.sparse_bam == 'y')
        else:
            print >> sys.stderr, 'ERROR: option sparse_bam should have value y or n, but has %s' % options.sparse_bam

        CFG['insert_intron_iterations'] = options.iterations
        if options.spladderfile != '-':
            CFG['spladder_infile'] = options.spladderfile

        ### settings for the alt splice part
        CFG['same_genestruct_for_all_samples'] = (options.same_genome == 'y')
        if options.replicates != '-':
            CFG['replicate_idxs'] = [int(x) for x in options.replicates.split(',')]
        CFG['curate_alt_prime_events'] = (options.curate_alt_prime == 'y')

        ### open log file, if specified
        if options.logfile != '-':
            CFG['log_fname'] = options.logfile
            CFG['fd_log'] = open(options.logfile, 'w')
        else:
            CFG['log_fname'] = 'stdout'
            CFG['fd_log'] = sys.stdout

        #if options.user != '-':
        #    CFG['user_settings'] = options.user

        ### alt splice analysis
        CFG['run_as_analysis'] = (options.extract_as == 'y')
        
        ### mandatory parameters for main spladder
        if options.bams == '-':
            print >> sys.stderr, 'ERROR: please provide the mandatory parameter: bam files\n\n'
            options.parser.print_help()
            sys.exit(2)
        else:
            CFG['bam_fnames'] = options.bams.strip(',').split(',')
            ### check existence of files
            for fname in CFG['bam_fnames']:
                if not os.path.isfile(fname):
                    print >> sys.stderr, 'ERROR: Input file %s can not be found\n\n' % fname
                    sys.exit(2)

        if options.annotation == '-':
            print >> sys.stderr, 'ERROR: please provide the mandatory parameter: annotation\n\n'
            options.parser.print_help()
            sys.exit(2)
        elif not os.path.isfile(options.annotation):
            print >> sys.stderr, 'ERROR: Annotation file %s can not be found\n\n' % options.annotation
            sys.exit(2)
        else:
            CFG['anno_fname'] = options.annotation
        
        if options.refstrain != '-':
            CFG['reference_strain'] = options.refstrain
            ref_tag = '%s:' % options.refstrain

        ### rproc options
        if options.pyproc == 'y':
            CFG['rproc'] = (options.pyproc == 'y')
            CFG['options_rproc'] = dict()
            CFG['options_rproc']['mem_req_resubmit']  = [30000, 60000, 80000]
            CFG['options_rproc']['time_req_resubmit'] = [60*60, 80*60, 90*60]
            CFG['options_rproc']['resubmit'] = 3
            CFG['options_rproc']['priority'] = 100
            CFG['options_rproc']['addpaths'] = CFG['paths']

        CFG['detect_edge_limit'] = options.detect_edge_limit


    if identity in ['main', 'test']:
        ### parallel processing
        CFG['parallel'] = options.parallel

        CFG['merge_strategy'] = options.merge
        CFG['read_length'] = options.readlen
        CFG['confidence_level'] = options.confidence

    if identity in ['main', 'test', 'viz']:
        if options.validate_sg in ['n', 'y']:
            CFG['validate_splicegraphs'] = (options.validate_sg == 'y')
        else:
            print >> sys.stderr, 'ERROR: validate_sg should have value y or n, but has %s' % options.validate_sg
            sys.exit(1)
    
    if identity == 'viz':
        CFG['event_id'] = options.event_id

        if options.transcripts in ['n', 'y']:
            CFG['plot_transcripts'] = (options.transcripts == 'y')
        else:
            print >> sys.stderr, 'ERROR: transcripts should have value y or n, but has %s' % options.transcripts
            sys.exit(1)

        if options.splicegraph in ['n', 'y']:
            CFG['plot_splicegraph'] = (options.splicegraph == 'y')
        else:
            print >> sys.stderr, 'ERROR: splicegraph should have value y or n, but has %s' % options.splicegraph
            sys.exit(1)

        CFG['plot_labels'] = options.labels

        if options.bams != '-':
            CFG['bam_fnames'] = options.bams.strip(':').split(':')
            for g, group in enumerate(CFG['bam_fnames']):
                CFG['bam_fnames'][g] = group.strip(',').split(',')
                ### check existence of files
                for fname in CFG['bam_fnames'][g]:
                    if not os.path.isfile(fname):
                        print >> sys.stderr, 'ERROR: Input file %s can not be found\n\n' % fname
                        sys.exit(2)


    if identity == 'test':
        CFG['multiTest'] = options.correction
        CFG['max_0_frac'] = options.max_0_frac
        CFG['min_count'] = options.min_count
        
        if options.non_alt_norm in ['n', 'y']:
            CFG['non_alt_norm'] = (options.non_alt_norm == 'y')
        else:
            print >> sys.stderr, 'ERROR: option non_alt_norm should have value y or n, but has %s' % options.non_alt_norm
            sys.exit(1)

        if options.low_memory in ['n', 'y']:
            CFG['low_memory'] = (options.low_memory == 'y')
        else:
            print >> sys.stderr, 'ERROR: option low_memory should have value y or n, but has %s' % options.low_memory
            sys.exit(1)

        if options.matlab in ['n', 'y']:
            CFG['is_matlab'] = (options.matlab == 'y')
        else:
            print >> sys.stderr, 'ERROR: option matlab should have value y or n, but has %s' % options.matlab
            sys.exit(1)

        if options.cap_exp_outliers in ['n', 'y']:
            CFG['cap_exp_outliers'] = (options.cap_exp_outliers == 'y')
        else:
            print >> sys.stderr, 'ERROR: option cap_exp_outliers should have value y or n, but has %s' % options.cap_exp_outliers
            sys.exit(1)

        if options.cap_outliers in ['n', 'y']:
            CFG['cap_outliers'] = (options.cap_outliers == 'y')
        else:
            print >> sys.stderr, 'ERROR: option cap_outliers should have value y or n, but has %s' % options.cap_outliers
            sys.exit(1)

        if options.conditionA == '-':
            print >> sys.stderr, 'ERROR: At least one sample for condition A required'
            sys.exit(1)
        if options.conditionB == '-':
            print >> sys.stderr, 'ERROR: At least one sample for condition B required'
            sys.exit(1)

        if options.diagnose_plots in ['n', 'y']:
            CFG['diagnose_plots'] = (options.diagnose_plots == 'y')
        else:
            print >> sys.stderr, 'ERROR: option diagnose_plots should have value y or n, but has %s' % options.diagnose_plots
            sys.exit(1)

        CFG['conditionA'] = options.conditionA.strip(',').split(',')
        CFG['conditionB'] = options.conditionB.strip(',').split(',')
        if len(CFG['conditionA']) > 0 and CFG['conditionA'][0].lower().endswith('txt'):
            CFG['conditionA'] = [str(x) for x in sp.loadtxt(CFG['conditionA'][0], dtype='str')]
        if len(CFG['conditionB']) > 0 and CFG['conditionB'][0].lower().endswith('txt'):
            CFG['conditionB'] = [str(x) for x in sp.loadtxt(CFG['conditionB'][0], dtype='str')]
        CFG['conditionA'] = [re.sub(r'(.[bB][aA][mM]|.[hH][dD][fF]5)|.[nN][pP][zZ]$', '', x) for x in CFG['conditionA']]
        CFG['conditionB'] = [re.sub(r'(.[bB][aA][mM]|.[hH][dD][fF]5)|.[nN][pP][zZ]$', '', x) for x in CFG['conditionB']]
        CFG['conditionA'] = [os.path.basename(x) for x in CFG['conditionA']]
        CFG['conditionB'] = [os.path.basename(x) for x in CFG['conditionB']]

    ### check if we got a list of bam files in a text file instead of a comma separated list
    if len(CFG['bam_fnames']) > 0:
        if identity == 'main' and CFG['bam_fnames'][0].split('.')[-1] == 'txt':
            CFG['bam_fnames'] = [str(x) for x in sp.atleast_1d(sp.loadtxt(CFG['bam_fnames'][0], dtype='str'))]
        elif identity == 'viz':
            for g, group in enumerate(CFG['bam_fnames']):
                if group[0].split('.')[-1] == 'txt':
                    CFG['bam_fnames'][g] = [str(x) for x in sp.atleast_1d(sp.loadtxt(group[0], dtype='str'))]

    ### assemble strain list
    CFG['samples'] = []
    CFG['strains'] = []
    if identity in ['viz']:
        for g, group in enumerate(CFG['bam_fnames']):
            CFG['strains'].append([]) 
            CFG['samples'].append([])
            for i in range(len(group)):
                CFG['samples'][-1].append(re.sub(r'(.[bB][aA][mM]|.[hH][dD][fF]5)$', '', group[i].split('/')[-1]))
                CFG['strains'][-1].append('%s%s' % (ref_tag, CFG['samples'][-1][-1]))
            CFG['strains'][-1] = sp.array(CFG['strains'][-1])
    else:
        for i in range(len(CFG['bam_fnames'])):
            if options.label != '-':
                CFG['samples'].append('%s_%s' % (options.label, re.sub(r'(.bam|.hdf5)$', '', CFG['bam_fnames'][i].split('/')[-1])))
            else:
                CFG['samples'].append(re.sub(r'(.[bB][aA][mM]|.[hH][dD][fF]5)$', '', CFG['bam_fnames'][i].split('/')[-1]))
            CFG['strains'].append('%s%s' % (ref_tag, CFG['samples'][-1]))
        CFG['strains'] = sp.array(CFG['strains'])

    ### adapt graph validation requirement to max number of samples
    CFG['sg_min_edge_count'] = min(CFG['sg_min_edge_count'], len(CFG['samples']))

    return CFG


def set_confidence_level(CFG):
# CFG = set_confidence_level(CFG),

    ### settings for accepted introns
    CFG['read_filter'] = dict()
    if CFG['confidence_level'] == 0:
        CFG['read_filter']['intron'] = 350000 
        CFG['read_filter']['exon_len'] = math.ceil(CFG['read_length'] * 0.10)
        CFG['read_filter']['mismatch'] = max(2, math.floor(CFG['read_length'] * 0.03)) 
        CFG['read_filter']['mincount'] = 1
    elif CFG['confidence_level'] == 1:
        CFG['read_filter']['intron'] = 350000 
        CFG['read_filter']['exon_len'] = math.ceil(CFG['read_length'] * 0.15)
        CFG['read_filter']['mismatch'] = max(1, math.floor(CFG['read_length'] * 0.02))
        CFG['read_filter']['mincount'] = 2
    elif CFG['confidence_level'] == 2:
        CFG['read_filter']['intron'] = 350000 
        CFG['read_filter']['exon_len'] = math.ceil(CFG['read_length'] * 0.20)
        CFG['read_filter']['mismatch'] = max(1, math.floor(CFG['read_length'] * 0.01)) 
        CFG['read_filter']['mincount'] = 2
    elif CFG['confidence_level'] == 3:
        CFG['read_filter']['intron'] = 350000 
        CFG['read_filter']['exon_len'] = math.ceil(CFG['read_length'] * 0.25)
        CFG['read_filter']['mismatch'] = 0
        CFG['read_filter']['mincount'] = 2

    ### settings for accepted cassette exons
    CFG['cassette_exon'] = dict()
    CFG['cassette_exon']['min_cassette_cov'] = 5 
    CFG['cassette_exon']['min_cassette_region'] = 0.9
    CFG['cassette_exon']['min_cassette_rel_diff'] = 0.5 

    ### settings for accepted intron retentions
    if not 'intron_retention' in CFG:
        CFG['intron_retention'] = dict()
    if CFG['confidence_level'] == 0:
      CFG['intron_retention']['min_retention_cov'] = 1
      CFG['intron_retention']['min_retention_region'] = 0.75 
      CFG['intron_retention']['min_retention_rel_cov'] = 0.1
      CFG['intron_retention']['max_retention_rel_cov'] = 2 
      CFG['intron_retention']['min_retention_max_exon_fold_diff'] = 4
    elif CFG['confidence_level'] == 1:
      CFG['intron_retention']['min_retention_cov'] = 2
      CFG['intron_retention']['min_retention_region'] = 0.75
      CFG['intron_retention']['min_retention_rel_cov'] = 0.1
      CFG['intron_retention']['max_retention_rel_cov'] = 1.2
      CFG['intron_retention']['min_retention_max_exon_fold_diff'] = 4 
    elif CFG['confidence_level'] == 2:
      CFG['intron_retention']['min_retention_cov'] = 5 
      CFG['intron_retention']['min_retention_region'] = 0.9 
      CFG['intron_retention']['min_retention_rel_cov'] = 0.2 
      CFG['intron_retention']['max_retention_rel_cov'] = 1.2 
      CFG['intron_retention']['min_retention_max_exon_fold_diff'] = 4 
    elif CFG['confidence_level'] == 3:
      CFG['intron_retention']['min_retention_cov'] = 10 
      CFG['intron_retention']['min_retention_region'] = 0.9  
      CFG['intron_retention']['min_retention_rel_cov'] = 0.2 
      CFG['intron_retention']['max_retention_rel_cov'] = 1.2 
      CFG['intron_retention']['min_retention_max_exon_fold_diff'] = 4 

    CFG['intron_retention']['read_filter'] = CFG['read_filter'] 

    return CFG

