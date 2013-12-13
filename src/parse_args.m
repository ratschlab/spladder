function CFG = parse_args(ARGS, CFG)

%%% are we running Octave or Matlab?
if size(ver('Octave'), 1)
    IS_OCT = 1;
else
    IS_OCT = 0;
end;
CFG.IS_OCT = IS_OCT;

%%% turn off warnings
if IS_OCT,
    warning('off', 'Octave:precedence-change');
    warning('off', 'Octave:function-name-clash');
    warning('off', '');
    warning('off', 'Octave:num-to-str');
    warning('off', 'Octave:function-name-clash');
    warning('off', 'Octave:divide-by-zero');
    warning('off', 'Octave:future-time-stamp');
    warning('off', 'solve_qp:constraints');
    warning('off', 'Octave:assign-as-truth-value');
    warning('off', 'Octave:matlab-incompatible');
else
    warning('off', 'MATLAB:typeaheadBufferOverflow');
end;

%%% Split Parameter string
if IS_OCT,
    ARGS = strsplit(ARGS, ';');
else
    ARGS = regexp(ARGS, ';', 'split');
end;

%%% Split values in parameter fields
for i = 1:length(ARGS),
    
    if isempty(ARGS{i}),
        continue;
    end;

    if IS_OCT,
        VALS = strsplit(ARGS{i}, ':');
    else
        VALS = regexp(ARGS{i}, ':', 'split');
    end;

    if length(VALS) > 2
        error(['ERROR: more than one field for variable: ' VALS{1} ':' VALS{:} '\n\tInput arguments MUST NOT contain colons!'])
    end;

    if strcmp(VALS{1}(1:2), 'S_'),
        eval([VALS{1} ' = "VALS{2}";']);
    elseif strcmp(VALS{1}(1:2), 'F_'),
        if strcmp(VALS{2}, 'y'),
            eval([VALS{1} ' = 1;']);
        elseif strcmp(VALS{2}, 'n'),
            eval([VALS{1} ' = 0;']);
        else
            error(sprintf('Unknown option for booloean input flag: %s\nplease use n for no or y for yes\n', VALS{2}));
        end;
    elseif strcmp(VALS{1}(1:2), 'I_'),
        eval([VALS{1} ' = str2num(VAL{2});']);
    else
        error(['ERROR: Unknown variable:' VALS{1} ':' VALS{2} '\n']);
    end;
end;

%%% get switches
CFG.do_insert_intron_retentions = F_INSERT_IR;
CFG.do_insert_cassette_exons = F_INSERT_CE;
CFG.do_insert_intron_edges = F_INSERT_IE;
CFG.do_remove_short_exons = F_REMOVE_SE;
CFG.do_infer_splice_graph = F_INFER_SG;
CFG.verbose = F_VERBOSE;
CFG.debug = F_DEBUG;

CFG.insert_intron_iterations = I_INSERT_INTRON_ITER;
CFG.confidence_level = I_CONFIDENCE;

%%% settings for the alt splice part
CFG.merge_strategy = S_MERGE_STRATEGY;
CFG.validate_splicegraphs = F_VALIDATE_SG;
CFG.same_genestruct_for_all_samples = F_SHARE_GENESTRUCT;
if IS_OCT,
    CFG.replicate_idxs = strsplit(S_REPLICATE_IDX, ',');
else
    CFG.replicate_idxs = regexp(S_REPLICATE_IDX, ',', 'split');
end;
CFG.curate_alt_prime_events = F_CURATE_ALTPRIME;

%%% open log file, if specified
if isempty(S_LOG_FNAME),
    CFG.log_fname = '';
%    CFG.fd_log = 1;
else
    CFG.log_fname = S_LOG_FNAME;
%    CFG.fd_log = fopen(S_LOG_FNAME, 'w');
end;

CFG.anno_fname = S_ANNO_FNAME;
CFG.out_dirname = S_OUT_DIRNAME;
CFG.user_settings = S_USER_FNAME;

CFG.no_reset_conf = 0;

if IS_OCT,
    CFG.bam_fnames = strsplit(S_BAM_FNAME, ',');
else
    CFG.bam_fnames = regexp(S_BAM_FNAME, ',', 'split');
end;

CFG.reference_strain = S_REFERENCE_STRAIN;
if strcmp(S_REFERENCE_STRAIN, '-'),
    ref_tag = '';
else
    ref_tag = sprintf('%s:', S_REFERENCE_STRAIN);
end;

CFG.list_config = {};
CFG.samples = {};

for i = 1:length(CFG.bam_fnames),
    if ~strcmp(S_EXPERIMENT_LABEL, '-'),
        CFG.samples{end + 1} = [S_EXPERIMENT_LABEL '_' regexprep(regexprep(CFG.bam_fnames{i}, '.*/', ''), '.bam', '')]; 
    else
        CFG.samples{end + 1} = [regexprep(regexprep(CFG.bam_fnames{i}, '.*/', ''), '.bam', '')]; 
    end;
    CFG.strains{end + 1} = [ref_tag CFG.samples{end}] ; 
end;

%%% rproc options
CFG.rproc = F_RPROC;
CFG.options_rproc = struct()
CFG.options_rproc.mem_req_resubmit  = [30000 60000 80000];
CFG.options_rproc.time_req_resubmit = [60*60 80*60 90*60];
CFG.options_rproc.resubmit = 3;
CFG.options_rproc.priority = 100;
CFG.options_rproc.addpaths = CFG.paths;

%%% set defaults for now
CFG.do_prune = 0;
CFG.do_gen_isoforms = 0;
CFG.do_merge_all = 0;

CFG.sg_min_edge_count = 1;

