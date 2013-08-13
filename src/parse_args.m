function CFG = parse_args(ARGS, CFG)

%%% are we running Octave or Matlab?
if size(ver('Octave'), 1)
    IS_OCT = 1;
else
    IS_OCT = 0;
end;

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
    elseif strcmp(VALS{1}(1:2), 'I_'),
        eval([VALS{1} ' = str2num(VALS{2});']);
    else
        error(['ERROR: Unknown variable:' VALS{1} ':' VALS{2} '\n']);
    end;
end;

%%% get switches
CFG.do_insert_intron_retentions = I_INSERT_IR;
CFG.do_insert_cassette_exons = I_INSERT_CE;
CFG.do_insert_intron_edges = I_INSERT_IE;
CFG.do_remove_short_exons = I_REMOVE_SE;
CFG.do_infer_splice_graph = I_INFER_SG;
CFG.verbose = I_VERBOSE;
CFG.debug = I_DEBUG;

CFG.insert_intron_iterations = I_INSERT_INTRON_ITER;
CFG.confidence_level = I_CONFIDENCE;

%%% open log file, if specified
if isempty(S_LOG_FNAME),
    CFG.fd_log = 1;
else
    CFG.fd_log = fopen(S_LOG_FNAME, 'w');
end;
CFG.anno_fname = S_ANNO_FNAME;
CFG.out_dir = S_OUT_DIR;
CFG.user_settings = S_USER_FNAME;

if IS_OCT,
    CFG.bam_fnames = strsplit(S_BAM_FNAME, ',');
else
    CFG.bam_fnames = regexp(S_BAM_FNAME, ',', 'split');
end;

