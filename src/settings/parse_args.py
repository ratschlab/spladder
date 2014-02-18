import sys
import scipy as sp

def parse_args(options):

    ### load all default settings
    CFG = default_settings()

    ### get switches
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

    if options.var_aware in ['n', 'y']:
        CFG['var_aware'] = (options.var_aware == 'y')
    else:
        print >> sys.stderr, 'ERROR: option var_aware should have value y or n, but has %s' % options.var_aware
        sys.exit(1)

    CFG['insert_intron_iterations'] = options.iterations
    CFG['confidence_level'] = options.confidence
    CFG['spladder_infile'] = options.spladderfile

    ### settings for the alt splice part
    CFG['merge_strategy'] = options.merge
    CFG['validate_splicegraphs'] = (options.validate_sg == 'y')
    CFG['same_genestruct_for_all_samples'] = (options.same_genome == 'y')
    CFG['replicate_idxs'] = options.replicates.split(',')
    CFG['curate_alt_prime_events'] = (options.curate_alt_prime == 'y')

    ### open log file, if specified
    if options.logfile != '-':
        CFG['log_fname'] = options.logfile
        CFG['fd_log'] = open(options.logfile, 'w')
    else:
        CFG['log_fname'] = 'stdout'
        CFG['fd_log'] = sys.stdout

    if options.user != '-':
        CFG['user_settings'] = options.user

    CFG['read_length'] = options.readlen

    ### alt splice analysis
    CFG['run_as_analysis'] = (options.extract_as == 'y')
    CFG['event_types'] = options.event_types.split(',')

    ### mandatory parameters
    CFG['bam_fnames'] = options.bams.split(',')
    CFG['anno_fname'] = options.annotation
    CFG['out_dirname'] = options.outdir

    ### check if we got a list of bam files in a text file instead of a comma separated list
    if CFG['bam_fnames'][0].split('.')[-1] == 'txt':
        CFG['bam_fnames'] = list(sp.loadtxt(CFG.bam_fnames[0], dtype='str'))

    if options.refstrain != '-':
        CFG['reference_strain'] = options.refstrain
        ref_tag = '%s:' % options.refstrain
    else:
        ref_tag = ''

    ### assemble strain list
    CFG.list_config = {};
    CFG.samples = {};
    CFG.strains = {};
    for i in range(len(CFG['bam_fnames'])):
        if options.label != '-':
            CFG['samples'].append('%s_%s' % (options.label, CFG['bam_fnames'][i].split('/')[-1].replace('.bam', '')))
        else:
            CFG['samples'].append(CFG['bam_fnames'][i].split('/')[-1].replace('.bam', ''))
        CFG['strains'].append('%s%s' % (ref_tag, CFG['samples'][-1]))

    ### rproc options
    if options.parallel == 'y':
        CFG['rproc'] = (options.parallel == 'y')
        CFG['options_rproc'] = dict()
        CFG['options_rproc']['mem_req_resubmit']  = [30000, 60000, 80000]
        CFG['options_rproc']['time_req_resubmit'] = [60*60, 80*60, 90*60]
        CFG['options_rproc']['resubmit'] = 3
        CFG['options_rproc']['priority'] = 100
        CFG['options_rproc']['addpaths'] = CFG['paths']

    return CFG
