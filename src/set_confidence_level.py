import math

def set_confidence_level(CFG):
# CFG = set_confidence_level(CFG),

    ### settings for accepted introns
    CFG['read_filter'] = dict()
    if CFG['confidence_level'] == 0:
        CFG['read_filter']['intron'] = 20000 
        CFG['read_filter']['exon_len'] = math.ceil(CFG['read_length'] * 0.10)
        CFG['read_filter']['mismatch'] = math.floor(CFG['read_length'] * 0.03) 
        CFG['read_filter']['mincount'] = 1
    elif CFG['confidence_level'] == 1:
        CFG['read_filter']['intron'] = 20000 
        CFG['read_filter']['exon_len'] = math.ceil(CFG['read_length'] * 0.15)
        CFG['read_filter']['mismatch'] = math.floor(CFG['read_length'] * 0.02)
        CFG['read_filter']['mincount'] = 2
    elif CFG['confidence_level'] == 2:
        CFG['read_filter']['intron'] = 20000 
        CFG['read_filter']['exon_len'] = math.ceil(CFG['read_length'] * 0.20)
        CFG['read_filter']['mismatch'] = math.floor(CFG['read_length'] * 0.01) 
        CFG['read_filter']['mincount'] = 3
    elif CFG['confidence_level'] == 3:
        CFG['read_filter']['intron'] = 20000 
        CFG['read_filter']['exon_len'] = math.ceil(CFG['read_length'] * 0.25)
        CFG['read_filter']['mismatch'] = 0
        CFG['read_filter']['mincount'] = 5

    ### settings for accepted cassette exons
    CFG['cassette_exon'] = []
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

