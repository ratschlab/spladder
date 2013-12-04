function ev = verify_all_events(ev, strain_idx, list_bam, verify, event_type, conf_filter, is_half_open) ;
% ev = verify_all_events(ev, strain_idx, list_bam, verify, event_type, conf_filter, is_half_open) ;

if nargin < 7,
    is_half_open = 0;
end;

out_fn = '' ;

fieldnames = {'intron', 'exon1', 'exon2', 'exon_alt1', 'exon_alt2', 'exon_const', 'intron1', 'intron2', 'exon', 'exon_pre', 'exon_aft', 'exons'};

%%% set parameters if called by rproc
if nargin==1,
    PAR = ev ;
    ev = PAR.ev ;
    strain_idx = PAR.strain_idx ;
    list_bam = PAR.list_bam ;
    verify = PAR.verify ;
    if isfield(PAR, 'out_fn'),
        out_fn = PAR.out_fn ;
    end ;
    if_half_open = 0;
    if isfield(PAR, 'is_half_open')
        is_half_open = PAR.is_half_open;
    end ;
    event_type = PAR.event_type;
    conf_filter = PAR.conf_filter;
end ;

%%% verify the events if demanded
if verify,
    for j = 1:length(strain_idx),
        s_idx = strain_idx(j);
        fprintf('%i/%i\r', j, length(strain_idx));
        for i = 1:length(ev),
            ev_tmp = ev(i) ;
            for f_idx = 1:length(fieldnames),
                field = fieldnames{f_idx};
                if isfield(ev, field) && size(ev_tmp.(field), 1) > 1,
                    ev_tmp.(field) = ev_tmp.(field)(s_idx, :);
                end;
            end;
            if strcmp(event_type, 'exon_skip'),
                [ev(i).verified(j,:), ev(i).info(j)] = verify_exon_skip(ev_tmp, list_bam(:,s_idx), is_half_open, conf_filter) ;
            elseif strcmp(event_type, 'alt_3prime') || strcmp(event_type, 'alt_5prime') || strcmp(event_type, 'alt_prime'),
                [ev(i).verified(j,:), ev(i).info(j)] = verify_alt_prime(ev_tmp, list_bam(:,s_idx), is_half_open, conf_filter) ;
            elseif strcmp(event_type, 'intron_retention'),
                [ev(i).verified(j,:), ev(i).info(j)] = verify_intron_retention(ev_tmp, list_bam(:,s_idx), is_half_open, conf_filter) ;
            elseif strcmp(event_type, 'mult_exon_skip'),
                [ev(i).verified(j,:), ev(i).info(j)] = verify_mult_exon_skip(ev_tmp, list_bam(:,s_idx), is_half_open, conf_filter) ;
            end;
        end ;
    end ;
end ;

if ~isempty(out_fn),
    save(out_fn, 'ev');
end;
