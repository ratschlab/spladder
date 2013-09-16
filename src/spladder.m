function genes = spladder(ARGS)
% function genes = spladder(ARGS)

% TODO: switches do_prune, min_orflen, find_orfs

%%% parse parameters from ARGS string
if isstruct(ARGS),
    CFG = ARGS;
else
    CFG = parse_args(ARGS, CFG);
end;

%%% add dependencies provided in config section
if isfield(CFG, 'paths'),
    for i = 1:length(CFG.paths),
        addpath(CFG.paths{i});
    end;
end;

%%% load confidence level settings
if CFG.no_reset_conf == 0,
    CFG = set_confidence_level(CFG);
end;

%%% check if result file exists and start gen graph step if necessary
if ~fexist(CFG.out_fname),
    if ~isfield(CFG, 'genes'),
        load(CFG.anno_fname, 'genes');
    else
        genes = CFG.genes ;
    end;
	genes = gen_graphs(genes, CFG) ;

	fprintf('saving genes to %s\n', CFG.out_fname) ;
	save(out_fname, 'genes') ;
else
	fprintf('gen_graphs already done\n%s exists\nloading genes from %s\n', CFG.out_fname, CFG.out_fname);
	load(out_fname, 'genes') ;
end ;

if do_prune,
    CFG.out_fname = strrep(CFG.out_fname, '.mat', '_pruned.mat');
end;
if do_prune && ~fexist(CFG.out_fname),
    genes = uniquify_splicegraph(genes);

    %%% prune graphs
    num_paths_before = count_all_paths(genes);
    genes = prune_graph(genes, bam_fnames);
    num_paths_after = count_all_paths(genes);

    %%% plot path reduction
    %plot(num_paths_before, num_paths_after, '.');
    %xlim([0, max(num_paths_before)]);
    %ylim([0, max(num_paths_before)]);
    %figure
    %plot(log10(num_paths_before), log10(num_paths_after), '.');
    %ylabel('num paths after (log10)');
    %xlabel('num paths before (log10)');
    %fn_fig = sprintf('%s/num_paths_reduction_log10.eps', fileparts(fn_out));
    %fprintf('print path reduction stats to file: %s', fn_fig);
    %print('-depsc', fn_fig);

    %%% save pruned genes
    fprintf('saving genes to %s\n', CFG.out_fname) ;
    save(CFG.out_fname, 'genes') ;
else
    if do_prune,
        fprintf('pruning already done; loading genes\n');
        load(CFG.out_fname, 'genes') ;
    else
        fprintf('No pruning requested!\n');
    end;
end;

if do_gen_isoforms,
    fn_out = strrep(fn_out, '.mat', '_with_isoforms.mat');
end;
if do_gen_isoforms &&  ~fexist(fn_out)
    %%% generate isoforms
    fprintf('Generating all isoforms\n') ;
    genes = generate_isoforms(genes, PAR.conf);

    %%% re-constitute splicing graph
    fprintf('re-constituting simplified splice graph from generated transcripts\n') ;
    genes_unsimplified = genes ;
    genes = splice_graph(genes, conf.do_infer_splice_graph);

    %%% save splicing graph with isoforms
    fprintf('saving genes to %s\n', fn_out) ;
    save(fn_out, 'genes', 'genes_unsimplified') ;
else
    if do_gen_isoforms,
        fprintf('Generating all isoforms already done; loading genes\n');
    else
        fprintf(' Generating all isoforms not requested\n');
    end;
end;


%%%%% short
return

if find_orfs,
  genes = closed_to_half_open(genes);
  
  for i=1:length(genes),
    for j=1:length(genes(i).transcripts)
      genes(i).cds_exons{j}=[] ;
    end ;
    genes(i).cds_exons(j+1:end)=[] ;
  end ;

  % remove duplicates
  for i=1:length(genes),
    genes(i) = reduce_number_of_transcripts(genes(i), []) ;
  end ;
  
  % infer orf
  rel_cutoff=1 ;
  genes = correct_tis_stop(genes, genome_info, min_orflen, rel_cutoff, 0);

  genes = half_open_to_closed(genes);
end;


