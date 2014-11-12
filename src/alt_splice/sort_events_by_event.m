function event_list = sort_events_by_event(event_list),
   
    if isfield(event_list, 'exon_col') && isfield(event_list, 'exon_pre_col') && isfield(event_list, 'exon_aft_col'),
        pre_tmp = vertcat(event_list.exon_pre_col);
        aft_tmp = vertcat(event_list.exon_aft_col);
        sort_list = [pre_tmp(:, 2) vertcat(event_list.exon_col) aft_tmp(:, 1)];
        clear pre_tmp aft_tmp;
    elseif isfield(event_list, 'exons_col') && isfield(event_list, 'exon_pre_col') && isfield(event_list, 'exon_aft_col'),
        sort_list = [];
        for i = 1:length(event_list)
            sort_list = [sort_list; [event_list(i).exon_pre_col(2) event_list(i).exons_col(1:2) event_list(i).exons_col(end-1:end) event_list(i).exon_aft_col(1)]];
        %    sort_list(end + 1) = str2num([sprintf('%i', event_list(i).exon_pre_col(2)), sprintf('%i', event_list(i).exons), sprintf('%i', event_list(i).exon_aft_col(1))]); 
        end;
        %sort_list = sort_list';
    elseif isfield(event_list, 'intron_col') && isfield(event_list, 'exon1_col') && isfield(event_list, 'exon2_col'),
        sort_list = [vertcat(event_list.intron_col)];
    elseif isfield(event_list, 'intron1_col') && isfield(event_list, 'intron2_col') && isfield(event_list, 'exon_alt1_col') && isfield(event_list, 'exon_alt2_col') && isfield(event_list, 'exon_const_col'),
        int_tmp = [vertcat(event_list.intron1_col) vertcat(event_list.intron2_col)];
        %%% keep introns half open for sorting
        int_tmp(:, 2) = int_tmp(:, 2) + 1;
        int_tmp(:, 4) = int_tmp(:, 4) + 1;
        sort_list = zeros(size(int_tmp, 1), 3); 
        for i = 1:size(int_tmp, 1),
            sort_list(i, :) = sort(unique(int_tmp(i, :)));
        end;
        clear int_tmp;
    end;
    
    chr_list = [event_list.chr_num] ;
    strand_list = [event_list.strand];
    sort_list = [chr_list' double(strand_list') sort_list] ;
    [tmp, idx] = sortrows(sort_list) ;
    event_list = event_list(idx) ;

