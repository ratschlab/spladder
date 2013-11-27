function CFG = set_confidence_level_level(CFG),
%function CFG = set_confidence_level_level(CFG),

    %%% settings for accepted introns
    CFG.intron_filter=[] ;
    if CFG.confidence_level == 0,
        CFG.intron_filter.intron = 20000 ; 
        CFG.intron_filter.exon_len = 8;
        CFG.intron_filter.mismatch = 2 ;
        CFG.intron_filter.mincount = 1 ;
    elseif CFG.confidence_level == 1,
        CFG.intron_filter.intron = 20000 ; 
        CFG.intron_filter.exon_len = 12;
        CFG.intron_filter.mismatch = 1 ;
        CFG.intron_filter.mincount = 2 ;
    elseif CFG.confidence_level == 2
        CFG.intron_filter.intron = 20000 ; 
        CFG.intron_filter.exon_len = 20;
        CFG.intron_filter.mismatch = 1 ;
        CFG.intron_filter.mincount = 5 ;
    elseif CFG.confidence_level == 3
        CFG.intron_filter.intron = 20000 ; 
        CFG.intron_filter.exon_len = 25;
        CFG.intron_filter.mismatch = 0 ;
        CFG.intron_filter.mincount = 10 ;
    end ;

    %%% settings for accepted cassette exons
    CFG.cassette_exon = [];
    CFG.cassette_exon.min_cassette_cov = 5 ; 
    CFG.cassette_exon.min_cassette_region = 0.9; 
    CFG.cassette_exon.min_cassette_rel_diff = 0.5; 

    %%% settings for accepted intron retentions
    CFG.intron_retention=[] ;
    if CFG.confidence_level == 0,
      CFG.intron_retention.min_retention_cov = 1;
      CFG.intron_retention.min_retention_region = 0.75; 
      CFG.intron_retention.min_retention_rel_cov = 0.1;
      CFG.intron_retention.max_retention_rel_cov = 2 ;
      CFG.intron_retention.min_retention_max_exon_fold_diff = 4; 
    elseif CFG.confidence_level == 1,
      CFG.intron_retention.min_retention_cov = 2;
      CFG.intron_retention.min_retention_region = 0.75; 
      CFG.intron_retention.min_retention_rel_cov = 0.1;
      CFG.intron_retention.max_retention_rel_cov = 1.2;
      CFG.intron_retention.min_retention_max_exon_fold_diff = 4; 
    elseif CFG.confidence_level == 2,
      CFG.intron_retention.min_retention_cov = 5 ;
      CFG.intron_retention.min_retention_region = 0.9 ; 
      CFG.intron_retention.min_retention_rel_cov = 0.2 ; 
      CFG.intron_retention.max_retention_rel_cov = 1.2 ; 
      CFG.intron_retention.min_retention_max_exon_fold_diff = 4 ; 
    elseif CFG.confidence_level == 3,
      CFG.intron_retention.min_retention_cov = 10 ;
      CFG.intron_retention.min_retention_region = 0.9 ; 
      CFG.intron_retention.min_retention_rel_cov = 0.2 ; 
      CFG.intron_retention.max_retention_rel_cov = 1.2 ; 
      CFG.intron_retention.min_retention_max_exon_fold_diff = 4 ; 
    end;

    CFG.intron_retention.read_filter = CFG.intron_filter ;

return

