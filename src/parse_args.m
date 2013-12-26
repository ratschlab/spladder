function CFG = parse_args(ARGS, CFG)

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

%%% load all default settings
default_settings

%%% Split Parameter string
ARGS = split_string(ARGS, ';');

%%% Split values in parameter fields
for i = 1:length(ARGS),
    
    if isempty(ARGS{i}),
        continue;
    end;

    VALS = split_string(ARGS{i}, ':');

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
if exist(F_INSERT_IR, 'var') == 1, CFG.do_insert_intron_retentions = F_INSERT_IR; end;
if exist(F_INSERT_CE, 'var') == 1, CFG.do_insert_cassette_exons = F_INSERT_CE; end;
if exist(F_INSERT_IE, 'var') == 1, CFG.do_insert_intron_edges = F_INSERT_IE; end;
if exist(F_REMOVE_SE, 'var') == 1, CFG.do_remove_short_exons = F_REMOVE_SE; end;
if exist(F_INFER_SG, 'var') == 1, CFG.do_infer_splice_graph = F_INFER_SG; end;
if exist(F_VERBOSE, 'var') == 1, CFG.verbose = F_VERBOSE; end;
if exist(F_DEBUG, 'var') == 1, CFG.debug = F_DEBUG; end;

if exist(I_INSERT_INTRON_ITER, 'var') == 1, CFG.insert_intron_iterations = I_INSERT_INTRON_ITER; end;
if exist(I_CONFIDENCE, 'var') == 1, CFG.confidence_level = I_CONFIDENCE; end;

%%% settings for the alt splice part
if exist(S_MERGE_STRATEGY, 'var') == 1, CFG.merge_strategy = S_MERGE_STRATEGY; end;
if exist(F_VALIDATE_SG, 'var') == 1, CFG.validate_splicegraphs = F_VALIDATE_SG; end;
if exist(F_SHARE_GENESTRUCT, 'var') == 1, CFG.same_genestruct_for_all_samples = F_SHARE_GENESTRUCT; end;
if exist(S_REPLICATE_IDX, 'var') == 1, CFG.replicate_idxs = split_string(S_REPLICATE_IDX, ','); end;
if exist(F_CURATE_ALTPRIME, 'var') == 1, CFG.curate_alt_prime_events = F_CURATE_ALTPRIME; end;

%%% open log file, if specified
if exist(S_LOG_FNAME, 'var') == 1,
    CFG.log_fname = S_LOG_FNAME;
    CFG.fd_log = fopen(S_LOG_FNAME, 'w');
end;

if exist(S_USER_FNAME, 'var') == 1, CFG.user_settings = S_USER_FNAME; end;

%%% mandatory parameters
CFG.bam_fnames = split_string(S_BAM_FNAME, ',');
CFG.anno_fname = S_ANNO_FNAME;
CFG.out_dirname = S_OUT_DIRNAME;

if exist(S_REFERENCE_STRAIN, 'var') == 1,
    CFG.reference_strain = S_REFERENCE_STRAIN;
    ref_tag = sprintf('%s:', S_REFERENCE_STRAIN);
else
    ref_tag = '';
end;

CFG.list_config = {};
CFG.samples = {};
for i = 1:length(CFG.bam_fnames),
    if exist(S_EXPERIMENT_LABEL, 'var') == 1,
        CFG.samples{end + 1} = [S_EXPERIMENT_LABEL '_' regexprep(regexprep(CFG.bam_fnames{i}, '.*/', ''), '.bam', '')]; 
    else
        CFG.samples{end + 1} = [regexprep(regexprep(CFG.bam_fnames{i}, '.*/', ''), '.bam', '')]; 
    end;
    CFG.strains{end + 1} = [ref_tag CFG.samples{end}] ; 
end;

%%% rproc options
if exist(F_PROC, 'var') == 1,
    CFG.rproc = F_RPROC;
    CFG.options_rproc = struct()
    CFG.options_rproc.mem_req_resubmit  = [30000 60000 80000];
    CFG.options_rproc.time_req_resubmit = [60*60 80*60 90*60];
    CFG.options_rproc.resubmit = 3;
    CFG.options_rproc.priority = 100;
    CFG.options_rproc.addpaths = CFG.paths;
end;
