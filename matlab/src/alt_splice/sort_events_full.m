function event_list = sort_events_full(event_list),

    if length(event_list) == 0,
        return;
    end;

    if strcmp(event_list(1).event_type, 'exon_skip'), 
        sort_list = [vertcat(event_list.exon_pre_col) vertcat(event_list.exon_col) vertcat(event_list.exon_aft_col)];
    elseif strcmp(event_list(1).event_type, 'mult_exon_skip'), 
        tmp =[];
        for i = 1:length(event_list),
            tmp(end + 1, :) = [event_list(i).exons_col(1, 1) event_list(i).exons_col(end, end)];
        end;
        sort_list = [vertcat(event_list.exon_pre_col) tmp vertcat(event_list.exon_aft_col)];
    elseif strcmp(event_list(1).event_type, 'intron_retention'), 
        sort_list = [vertcat(event_list.exon1_col) vertcat(event_list.intron_col) vertcat(event_list.exon2_col)];
    elseif strcmp(event_list(1).event_type, 'alt_5prime') || strcmp(event_list(1).event_type, 'alt_3prime'), 
        sort_list = [vertcat(event_list.exon_const_col) vertcat(event_list.intron1_col) vertcat(event_list.exon_alt1_col) vertcat(event_list.intron2_col) vertcat(event_list.exon_alt2_col)];
    elseif strcmp(event_list(1).event_type, 'mutex_exons'),
        sort_list = [vertcat(event_list.exon_pre_col) vertcat(event_list.exon1_col) vertcat(event_list.exon2_col) vertcat(event_list.exon_aft_col)];
    end;
    
    chr_list = [event_list.chr_num] ;
    strand_list = [event_list.strand];
    sort_list = [chr_list' double(strand_list') sort_list] ;
    [tmp, idx] = sortrows(sort_list) ;
    event_list = event_list(idx) ;

