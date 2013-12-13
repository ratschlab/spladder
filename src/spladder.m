function genes = spladder(ARGS)
% function genes = spladder(ARGS)

%%% parse parameters from ARGS string
if isstruct(ARGS),
    CFG = ARGS;
else
    CFG = parse_args(ARGS, CFG);
end;

%%% add dependencies provided in config section
if isfield(CFG, 'paths'),
    addpath(CFG.paths{:});
end;

%%% load confidence level settings
if CFG.no_reset_conf == 0,
    CFG = set_confidence_level(CFG);
end;

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

%%% create spladder sub-directory
[tmp, tmp] = mkdir(CFG.out_dirname,'spladder');

for idx = idxs,
    CFG_ = CFG;
    if ~strcmp(CFG.merge_strategy, 'merge_bams'),
        CFG.bam_fnames = CFG.bam_fnames(idx);
        CFG.samples = CFG.samples(idx);
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
            jobinfo(job_nr) = rproc('spladder_core', CFG, 10000, CFG.options_rproc, 40*60) ;
            job_nr = job_nr + 1;
        else
            spladder_core(CFG);
        end;
    end;

    CFG = CFG_;
end;

%%% collect results after parallelization
if CFG.rproc,
    jobinfo = rproc_wait(jobinfo, 30, 1, 1) ;
end;

%%% merge parts if necessary
if strcmp(CFG.merge_strategy, 'merge_graphs'),
    run_merge(CFG);
end;

%%% handle alternative splicing part
alt_genes_collect(CFG);

for idx = 1:length(CFG.event_types),
    alt_genes_analyze(CFG, CFG.event_types{idx});
end;
