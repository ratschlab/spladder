function ev = verify_all_events(ev, event_type, CFG) ;
% ev = verify_all_events(ev, event_type, CFG) ;

out_fn = '' ;

fieldnames = {'intron', 'exon1', 'exon2', 'exon_alt1', 'exon_alt2', 'exon_const', 'intron1', 'intron2', 'exon', 'exon_pre', 'exon_aft', 'exons'};

%%% set parameters if called by rproc
if nargin==1,
    PAR = ev ;
    ev = PAR.ev ;
    if isfield(PAR, 'out_fn'),
        out_fn = PAR.out_fn ;
    end ;
    event_type = PAR.event_type;
    CFG = PAR.CFG;
end ;

%%% load counts
prune_tag = '';
if CFG.do_prune,
    prune_tag = '_pruned';
end;
validate_tag = '';
if CFG.validate_splicegraphs,
    validate_tag = '.validated';
end;
fn_count = sprintf('%s/spladder/genes_graph_conf%i.%s%s%s.count.mat', CFG.out_dirname, CFG.confidence_level, CFG.merge_strategy, validate_tag, prune_tag);
load(fn_count);
fn_genes = sprintf('%s/spladder/genes_graph_conf%i.%s%s%s.mat', CFG.out_dirname, CFG.confidence_level, CFG.merge_strategy, validate_tag, prune_tag);
load(fn_genes);

%%% verify the events if demanded
if CFG.verify_alt_events,
    for s_idx = 1:length(CFG.strains),
        fprintf('%i/%i\r', s_idx, length(CFG.strains));
        for i = 1:length(ev),
            ev_tmp = ev(i) ;
            for f_idx = 1:length(fieldnames),
                field = fieldnames{f_idx};
                if isfield(ev, field) && size(ev_tmp.(field), 1) > 1,
                    ev_tmp.(field) = ev_tmp.(field)(s_idx, :);
                end;
            end;
            if strcmp(event_type, 'exon_skip'),
                [ev(i).verified(s_idx,:), ev(i).info(s_idx)] = verify_exon_skip(ev_tmp, genes, counts(s_idx, :), CFG) ;
            elseif strcmp(event_type, 'alt_3prime') || strcmp(event_type, 'alt_5prime') || strcmp(event_type, 'alt_prime'),
                [ev(i).verified(s_idx,:), ev(i).info(s_idx)] = verify_alt_prime(ev_tmp, genes, counts(s_idx, :), CFG) ;
            elseif strcmp(event_type, 'intron_retention'),
                [ev(i).verified(s_idx,:), ev(i).info(s_idx)] = verify_intron_retention(ev_tmp, genes, counts(s_idx, :), CFG) ;
            elseif strcmp(event_type, 'mult_exon_skip'),
                [ev(i).verified(s_idx,:), ev(i).info(s_idx)] = verify_mult_exon_skip(ev_tmp, genes, counts(s_idx, :), CFG) ;
            end;
        end ;
    end ;
end ;

if ~isempty(out_fn),
    save(out_fn, 'ev');
end;
