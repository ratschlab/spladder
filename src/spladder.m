function genes = spladder(ARGS)
% function genes = spladder(ARGS)

%%% parse parameters from ARGS string
if isstruct(ARGS),
    CFG = ARGS;
    CFG = parse_args('', CFG);
else
    CFG = parse_args(ARGS);
end;

%%% add dependencies provided in config section
if isfield(CFG, 'paths'),
    addpath(CFG.paths{:});
end;

%%% load confidence level settings
if CFG.no_reset_conf == 0,
    CFG = set_confidence_level(CFG);
end;

%%% do not compute components of merged set, if result file already exists
fn_out_merge = '';
prune_tag = '';
if strcmp(CFG.merge_strategy, 'merge_graphs'),
    if CFG.do_prune,
        prune_tag = '_pruned';
    end;
    fn_out_merge = sprintf('%s/spladder/genes_graph_conf%i.%s%s.mat', CFG.out_dirname, CFG.confidence_level, CFG.merge_strategy, prune_tag);
end;

if ~isfield(CFG, 'spladder_infile') && ~exist(fn_out_merge, 'file'),
    %%% iterate over files, if merge strategy is single
    if strcmp(CFG.merge_strategy, 'single') || strcmp(CFG.merge_strategy, 'merge_graphs'),
        idxs = 1:length(CFG.samples);
    else
        idxs = 1;
    end;
        
    %%% set parallelization
    if CFG.rproc,
        jobinfo = rproc_empty() ;
        job_nr = 1;
    end;

    %%% create out-directory
    [tmp, tmp] = mkdir(CFG.out_dirname);

    %%% create spladder sub-directory
    [tmp, tmp] = mkdir(CFG.out_dirname,'spladder');

    %%% create truncation sub-directory, if necessary
    if CFG.detect_trunc,
        [tmp, tmp] = mkdir(CFG.out_dirname, 'truncations');
    end;

    for idx = idxs,
        CFG_ = CFG;
        if ~strcmp(CFG.merge_strategy, 'merge_bams'),
            CFG.bam_fnames = CFG.bam_fnames(idx);
            CFG.samples = CFG.samples(idx);
        end;
        %%% truncation detection mode?
        if CFG.detect_trunc,
            if ~strcmp(CFG.merge_strategy, 'merge_bams'),
                CFG.out_fname = sprintf('%s/truncations/genes_graph_conf%i.%s.tsv', CFG.out_dirname, CFG.confidence_level, CFG.samples{1});
            else
                CFG.out_fname = sprintf('%s/truncations/genes_graph_conf%i.%s.tsv', CFG.out_dirname, CFG.confidence_level, CFG.merge_strategy);
            end;
            if exist(CFG.out_fname, 'file'),
                fprintf(1, '%s exists\n', CFG.out_fname);
            else
                if CFG.rproc,
                    jobinfo(job_nr) = rproc('detect_truncations', CFG, 25000, CFG.options_rproc, 60) ;
                    job_nr = job_nr + 1;
                else
                    detect_truncations(CFG);
                end;
            end;
        else
            %%% spladder output files
            if ~strcmp(CFG.merge_strategy, 'merge_bams'),
                CFG.out_fname = sprintf('%s/spladder/genes_graph_conf%i.%s.mat', CFG.out_dirname, CFG.confidence_level, CFG.samples{1});
            else
                CFG.out_fname = sprintf('%s/spladder/genes_graph_conf%i.%s.mat', CFG.out_dirname, CFG.confidence_level, CFG.merge_strategy);
            end;

            %%% assemble out filename to check if we are already done
            fn_out = CFG.out_fname;
            if CFG.do_prune,
                fn_out = strrep(fn_out, '.mat', '_pruned.mat');
            end;
            if CFG.do_gen_isoforms,
                fn_out = strrep(fn_out, '.mat', '_with_isoforms.mat');
            end;

            if exist(fn_out, 'file'),
                fprintf('All result files already exist.\n');
            else
                if CFG.rproc,
                    jobinfo(job_nr) = rproc('spladder_core', CFG, 30000, CFG.options_rproc, 12*60) ;
                    job_nr = job_nr + 1;
                else
                    spladder_core(CFG);
                end;
            end;
        end;
        CFG = CFG_;
    end;

    %%% collect results after parallelization
    if CFG.rproc,
        jobinfo = rproc_wait(jobinfo, 30, 1, -1) ;
    end;

    if CFG.detect_trunc,
        return
    end;

    %%% merge parts if necessary
    if strcmp(CFG.merge_strategy, 'merge_graphs'),
        run_merge(CFG);
    end;
end;

%%% determine count output file
if ~isfield(CFG, 'spladder_infile'),
    if CFG.validate_splicegraphs,
        fn_in_count = sprintf('%s/spladder/genes_graph_conf%i.%s%s.validated.mat', CFG.out_dirname, CFG.confidence_level, CFG.merge_strategy, prune_tag);
    else
        fn_in_count = sprintf('%s/spladder/genes_graph_conf%i.%s%s.mat', CFG.out_dirname, CFG.confidence_level, CFG.merge_strategy, prune_tag);
    end;
else
    fn_in_count = CFG.spladder_infile;
end;
fn_out_count = [strrep(fn_in_count, '.mat', ''), '.count.mat'];

%%% count segment graph
if ~exist(fn_out_count, 'file'),
    count_graph_coverage_wrapper(fn_in_count, fn_out_count, CFG);
end;

%%% count intron coverage phenotype
if CFG.count_intron_cov,
    fn_out_intron_count = strrep(fn_out_count, 'mat', 'introns.mat');
    count_intron_coverage_wrapper(fn_in_count, fn_out_intron_count, CFG);
end;

%%% handle alternative splicing part
if CFG.run_as_analysis,
    alt_genes_collect(CFG);

    for idx = 1:length(CFG.event_types),
        alt_genes_analyze(CFG, CFG.event_types{idx});
    end;
end;
