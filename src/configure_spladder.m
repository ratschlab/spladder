function configure_spladder()

    %%% settings for adding new intron edges
    CFG.intron_edges = [] ;
    CFG.intron_edges.min_exon_len = 50;
    CFG.intron_edges.min_exon_len_remove = 1;
    CFG.intron_edges.vicinity_region = 40 ;
    CFG.intron_edges.insert_intron_retention = 1 ;
    CFG.intron_edges.gene_merges = 0 ;
    CFG.intron_edges.append_new_terminal_exons = 1;
    CFG.intron_edges.append_new_terminal_exons_len = 200;

    CFG.intron_window = 5000 ;

    CFG = set_confidene_level(CFG, 3);

    %%% settings for cassette exons
    CFG.cassette_exons = [];
    CFG.cassette_exons.min_cassette_cov = 5;
    CFG.cassette_exons.min_cassette_region = 0.9;
    CFG.cassette_exons.min_cassette_rel_diff = 0.5;

    %%% settings for short exon removal
    CFG.remove_exons = [];
	CFG.remove_exons.terminal_short_extend = 40 ;
	CFG.remove_exons.terminal_short_len = 10 ;
	CFG.remove_exons.min_exon_len = 50 ;
	CFG.remove_exons.min_exon_len_remove = 10 ;

    %%% settings for splice graph augmentations
    CFG.do_insert_intron_retentions = 1 ;
    CFG.do_insert_cassette_exons = 1;
    CFG.do_insert_intron_edges = 1;
    CFG.do_remove_short_exons = 0 ;
    CFG.do_infer_splice_graph = 0;

    CFG.insert_intron_iterations = 5;

    %%% set I/O and verbosity
    CFG.verbose = 1;
    CFG.debug = 0;
    CFG.fd_log = 1;

