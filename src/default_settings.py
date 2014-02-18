
def default_setings():

    ### import necessary paths 
    CFG['paths'] = []
    if 'SPLADDER_SRC_PATH' in os.environ:
        CFG['paths'].append(os.environ['SPLADDER_SRC_PATH'])
        CFG['paths'].append('%s/rproc' % os.environ['SPLADDER_SRC_PATH'])
        CFG['paths'].append('%s/utils' % os.environ['SPLADDER_SRC_PATH'])
        CFG['paths'].append('%s/alt_splice' % os.environ['SPLADDER_SRC_PATH')]
    else:
        ### infer src path from current script path
        #SPLADDER_SRC_PATH = fileparts(which('default_settings'))
        SPLADDER_SRC_PATH = os.getcwd()
        CFG['paths'].append(SPLADDER_SRC_PATH)
        CFG['paths'].append('%s/rproc' % SPLADDER_SRC_PATH)
        CFG['paths'].append('%s/utils' % SPLADDER_SRC_PATH)
        CFG['paths'].append('%s/alt_splice' % SPLADDER_SRC_PATH)

    if 'SPLADDER_PATH' in os.environ:
        CFG['paths'].append('%s/mex' % os.environ['SPLADDER_PATH'])
    else:
        #SPLADDER_SRC_PATH = fileparts(which('default_settings'))
        SPLADDER_SRC_PATH = os.getcwd()
        CFG['paths'].append('%s/../mex' % SPLADDER_SRC_PATH)
    
    if len(CFG['paths']) > 0:
        for p in paths:
            eval('import %s' % p )

    ### settings for adding new intron edges
    CFG['intron_edges'] = []
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
    CFG['merge_strategy'] = 'merge_bams' ### alternatives are: merge_graphs, single, merge_all
    CFG['confidence_level'] = 3
    CFG['validate_splicegraphs'] = 0
    CFG['same_genestruct_for_all_samples'] = 1
    CFG['curate_alt_prime_events'] = 1
    CFG['replicate_idxs'] = 1
    CFG['verify_alt_events'] = 1

    ### settings for verifying exon skips
    CFG['exon_skip'] = dict()
    CFG['exon_skip']['min_non_skip_count'] = 3
    CFG['exon_skip']['min_skip_count'] = 3
    CFG['exon_skip']['min_skip_rel_cov'] = 0.05
   #CFG['exon_skip']['max_exon_fold_diff'] = 4
   #CFG['exon_skip']['max_skip_rel_cov'] = 1.5
    CFG['exon_skip']['intron_tolerance'] = 0

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

    ### set I/O and verbosity
    CFG['verbose'] = 1
    CFG['debug'] = 0
    CFG['fd_log'] = 1

    CFG['sg_min_edge_count'] = 1
    CFG['no_reset_conf'] = 0

    CFG['do_prune'] = 0
    CFG['do_gen_isoforms'] = 0
    CFG['do_merge_all'] = 0

    CFG['is_half_open'] = 1

    CFG['run_as_analysis'] = 1
    CFG['event_types'] = ['exon_skip', 'intron_retention', 'alt_3prime', 'alt_5prime', 'mult_exon_skip']

    CFG['read_length'] = 36

    CFG['rproc'] = 0
    CFG['var_aware'] = 0

    return CFG
