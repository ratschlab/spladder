function event_list = sort_events_full(event_list),

    if len(event_list) == 0m
        return;
    end;

    %%% TODO switch by event type not by presence of certain fields
   
    if isfield(event_list, 'exon_col') && isfield(event_list, 'exon_pre_col') && isfield(event_list, 'exon_aft_col'),
        sort_list = [vertcat(event_list.exon_pre_col) vertcat(event_list.exon_col) vertcat(event_list.exon_aft_col)];
    elseif isfield(event_list, 'exons_col') && isfield(event_list, 'exon_pre_col') && isfield(event_list, 'exon_aft_col'),
        tmp =[];
        for i = 1:length(event_list),
            tmp(end + 1, :) = [event_list(i).exons_col(1, 1) event_list(i).exons_col(end, end)];
        end;
        sort_list = [vertcat(event_list.exon_pre_col) tmp vertcat(event_list.exon_aft_col)];
    elseif isfield(event_list, 'intron_col') && isfield(event_list, 'exon1_col') && isfield(event_list, 'exon2_col'),
        sort_list = [vertcat(event_list.exon1_col) vertcat(event_list.intron_col) vertcat(event_list.exon2_col)];
    elseif isfield(event_list, 'intron1_col') && isfield(event_list, 'intron2_col') && isfield(event_list, 'exon_alt1_col') && isfield(event_list, 'exon_alt2_col') && isfield(event_list, 'exon_const_col'),
        sort_list = [vertcat(event_list.exon_const_col) vertcat(event_list.intron1_col) vertcat(event_list.exon_alt1_col) vertcat(event_list.intron2_col) vertcat(event_list.exon_alt2_col)];
    end;
    
    chr_list = [event_list.chr_num] ;
    strand_list = [event_list.strand];
    sort_list = [chr_list' double(strand_list') sort_list] ;
    [tmp, idx] = sortrows(sort_list) ;
    event_list = event_list(idx) ;

