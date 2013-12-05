function [genes, inserted] = gen_graphs(genes, CFG)
% genes = gen_graphs(genes, CFG)

if nargin == 1,
    PAR = genes ;
    genes  = PAR.genes ;
    CFG = PAR.CFG ;
end ;

%%% init the stats for inserted elements
inserted.cassette_exon = 0 ;
inserted.intron_retention = 0 ;
inserted.intron_in_exon = 0 ;
inserted.alt_53_prime = 0 ;
inserted.exon_skip = 0 ;
inserted.gene_merge = 0 ;
inserted.new_terminal_exon = 0 ;

%%% init log stream
if isempty(CFG.log_fname),
    CFG.fd_log = 1;
else
    CFG.fd_log = fopen(CFG.log_fname, 'w');
end;

% build splice graph for all genes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(CFG.fd_log, 'Generate splice graph ...\n');
genes = splice_graph(genes, CFG);
fprintf(CFG.fd_log, '...done.\n\n');

% append list of introns supported by RNA-seq data to 
% the genes structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(CFG.fd_log, 'Load introns from file ...\n');
introns = get_intron_list(genes, CFG);
fprintf(CFG.fd_log, '...done.\n\n');

%%% check feasibility
introns = make_introns_feasible(introns, genes, CFG);

for i = 1:length(genes),
    genes(i).introns = introns(i, :);
end;

if CFG.do_insert_cassette_exons,
	fprintf(CFG.fd_log, '\t ...inserting cassette exons\n') ;
    CFG_ = CFG;
    if isfield(CFG, 'cassette_exon') && isfield(CFG.cassette_exon, 'read_filter')
        CFG.read_filter = CFG.cassette_exon.read_filter;
    end ;
	[genes, inserted_] = insert_cassette_exons(genes, CFG);
    inserted.cassette_exon = inserted_.cassette_exon;
    CFG = CFG_;
	fprintf(CFG.fd_log, '\tinserted %i casette exons\n', inserted.cassette_exon) ;
end; 

if CFG.do_insert_intron_retentions,
	fprintf(CFG.fd_log, '\t ... inserting intron retentions\n') ;
    CFG_ = CFG;
    if isfield(CFG.intron_retention, 'read_filter')
        CFG.read_filter = CFG.intron_retention.read_filter ;
    end ;
	[genes, inserted_] = insert_intron_retentions(genes, CFG);
    inserted.intron_retention = inserted_.intron_retention;

    CFG = CFG_;
	fprintf(CFG.fd_log, '\tinserted %i intron retentions\n', inserted.intron_retention) ;
end ;

if CFG.do_remove_short_exons,
	fprintf(CFG.fd_log, '... removing short exons\n') ;

	genes = remove_short_exons(genes, CFG);

	for i = 1:length(genes)
		for j = 1:size(genes(i).splicegraph{1},2)
			if ~(genes(i).splicegraph{1}(2, j) - genes(i).splicegraph{1}(1, j) >= CFG.remove_exons.min_exon_len_remove),
				warning('WARNING: could not remove all short exons') ;
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
    % re-set list of introns supported by RNA-seq data to 
    % the genes structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:length(genes),
        genes(i).introns = introns(i, :);
    end;

    for chr_idx = unique([genes.chr_num]),
        tmp_genes = genes([genes.chr_num] == chr_idx);
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isfield(CFG, 'insert_intron_iterations'),
            CFG.insert_intron_iterations = 5;
        end;
        for iter = 1:CFG.insert_intron_iterations,
            fprintf(CFG.fd_log, 'insert intron edges - chr %i - iteration %i/5\n', chr_idx, iter) ;
            [genes_mod, inserted_] = insert_intron_edges(tmp_genes, CFG);

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

%%% print summary to log file
fn = fieldnames(inserted);
fprintf(CFG.fd_log, 'Inserted:\n');
for i = 1:length(fn),
    fprintf(CFG.fd_log, '\t%s:\t%i\n', fn{i}, inserted.(fn{i}));
end;

if CFG.fd_log > 1,
    fclose(CFG.fd_log);
end;

return

