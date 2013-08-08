function [genes, inserted] = gen_graphs(genes, fn_bam, fn_genome_config, conf)
% genes = gen_graphs(genes, fn_bam, fn_genome_config, conf, log_fname)

if nargin == 1,
    PAR = genes ;
    genes  = PAR.genes ;
    fn_bam = PAR.fn_bam ;
    fn_genome_config = PAR.fn_genome_config ;
    conf = PAR.conf ;
    if isfield(PAR, 'log_fname'),
        log_fname = PAR.log_fname ;
    else
        log_fname = '';
    end;
end ;

if (nargin < 4 && nargin > 1) || ~isfield(CFG, 'intron_filter')
	CFG.intron_filter = [] ;
end ;
if (nargin < 4 && nargin > 1) || ~isfield(CFG, 'intron_retention')
	CFG.intron_retention = [] ;
end ;
if (nargin < 4 && nargin > 1) || ~isfield(CFG, 'intron_edges')
	CFG.intron_edges = [] ;
end ;
if (nargin < 4 && nargin > 1) || ~isfield(CFG, 'remove_exons')
	CFG.remove_exons = [] ;
end ;

%%% init the stats for inserted elements
inserted.cassette_exon = 0 ;
inserted.intron_retention = 0 ;
inserted.intron_in_exon = 0 ;
inserted.alt_53_prime = 0 ;
inserted.exon_skip = 0 ;
inserted.gene_merge = 0 ;
inserted.new_terminal_exon = 0 ;

% build splice graph for all genes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(CFG.fd_log, 'generate splice graph\n');
genes = splice_graph(genes, CFG);

% append list of introns supported by RNA-seq data to 
% the genes structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(CFG.fd_log, 'load introns from file\n');
genes = append_intron_list(genes, fn_bam, CFG);
%%% check feasibility
genes = make_introns_feasible(genes, fn_bam, CFG);

if CFG.do_insert_cassette_exons,
	fprintf(CFG.fd_log, 'inserting cassette exons\n') ;
    CFG_ = CFG;
    if isfield(CFG, 'cassette_exon') && isfield(CFG.cassette_exon, 'read_filter')
        disp('using CFG.cassette_exon.read_filter') ;
        CFG.intron_filter = CFG.cassette_exon.read_filter;
    end ;
	[genes, inserted] = insert_cassette_exons(genes, fn_bam, CFG);
    CFG = CFG_;
	fprintf(CFG.fd_log, 'inserted %i casette exons\n', inserted.cassette_exon) ;
end; 

if conf.do_insert_intron_retentions,
	fprintf(CFG.fd_log, 'inserting intron retentions\n') ;
    read_filter = conf.intron_filter ;
    CFG_ = CFG;
    if isfield(CFG.intron_retention, 'read_filter')
        disp('using CFG.intron_retention.read_filter') ;
        CFG.intron_filter = CFG.intron_retention.read_filter ;
    end ;
	[genes, inserted] = insert_intron_retentions(genes, fn_bam, CFG);
    CFG = CFG_;
	fprintf(CFG.fd_log, 'inserted %i intron retentions\n', inserted.intron_retention) ;
end ;

if conf.do_remove_short_exons,
	fprintf(CFG.fd_log, 'removing short exons\n') ;

	genes = remove_short_exons(genes, CFG);

	for i = 1:length(genes)
		for j = 1:size(genes(i).splicegraph{1},2)
			if ~(genes(i).splicegraph{1}(2, j) - genes(i).splicegraph{1}(1, j) >= CFG.remove_exons.min_exon_len_remove),
				warning('could not remove all short exons') ;
			end ;
		end ; 
	end ; 
end ;


% test all exons if the reading frame is larger if exon is skipped
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fprintf('find exons to skip to elongate reading frame\n');
%genes = insert_cds_exon_skips(genes, genome_info)
%genes = splice_graph(genes);


% sanity checking
for i = 1:length(genes),
	assert(all(genes(i).splicegraph{1}(1, :) <= genes(i).splicegraph{1}(2, :)));
end ;

if CFG.do_insert_intron_edges,
    % append list of introns supported by RNA-seq data to 
    % the genes structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(CFG.fd_log, 'load introns from file\n');
    %%% re-append intron list
    genes = append_intron_list(genes, fn_bam, CFG);
    %%% check feasibility
    genes = make_introns_feasible(genes, fn_bam, CFG);

    for chr_idx = unique([genes.chr_num]),
        tmp_genes = genes([genes.chr_num] == chr_idx);
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for iter = 1:CFG.insert_intron_iterations,
            fprintf(fd_log, 'insert intron edges - chr %i - iteration %i/5\n', chr_idx, iter) ;
            [genes_mod, inserted_] = insert_intron_edges(tmp_genes, fn_bam, conf, fd_log);

            inserted.intron_in_exon = inserted.intron_in_exon + inserted_.intron_in_exon ;
            inserted.alt_53_prime = inserted.alt_53_prime + inserted_.alt_53_prime ;
            inserted.exon_skip = inserted.exon_skip + inserted_.exon_skip ;
            inserted.gene_merge = inserted.gene_merge + inserted_.gene_merge ;
            inserted.new_terminal_exon = inserted.new_terminal_exon + inserted_.new_terminal_exon ;

            % in case any exon was inserted that already existed, we merge them into one exon 
            genes_mod = merge_duplicate_exons(genes_mod, 1) ;

            % inserted
            if isequal(genes_mod, tmp_genes)
                break ;
            end ;
            tmp_genes = genes_mod ;
        end ;
        genes([genes.chr_num] == chr_idx) = genes_mod;
    end ;
end ;

return

