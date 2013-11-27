function alt_genes_collect(CFG)
% alt_genes_collect(CFG)

    %%% add paths
    %addpath ~/git/projects/2012/mGene_core/experiments/thaliana/
    %addpath ~/git/projects/2012/mGene_core/experiments/thaliana/tools/

    %%% init empty event fields
    idx_intron_reten = {} ;
    intron_intron_reten = {} ;
    idx_exon_skip = {} ;
    exon_exon_skip = {} ;
    intron_reten_pos = {} ;
    exon_skip_pos = {} ;

    idx_alt_end_5prime = {} ;
    exon_alt_end_5prime = {} ;
    alt_end_5prime_pos = {} ;
    idx_alt_end_3prime = {} ;
    exon_alt_end_3prime = {} ;
    alt_end_3prime_pos = {} ;

    idx_mult_exon_skip = {} ;
    exon_mult_exon_skip = {} ;
    id_mult_exon_skip = {} ;
    mult_exon_skip_pos = {};

    validate_tag = '';
    if isfield(CFG, 'validate_splicegraphs') && CFG.validate_splicegraphs,
        validate_tag = '.validated';
    end;

    for i = 1:length(CFG.samples),
        if CFG.same_genestruct_for_all_samples == 1 && i == 2,
            break;
        end;

        strain = strains{i} ;

        for ridx = CFG.replicate_idxs,
            if length(CFG.replicate_idxs) > 1,
                rep_tag = sprintf('_R%i', r_idx);
            else
                rep_tag = '';
            end;

            if strcmp(CFG.merge_strategy, 'single')
                genes_fnames{ridx, i} = sprintf('%s/spladder/genes_graph_conf%i%s.%s.mat', CFG.result_dir, CFG.confidence_level, rep_tag, CFG.samples{i});
            else
                genes_fnames{ridx, i} = sprintf('%s/spladder/genes_graph_conf%i%s.merged%s.mat', CFG.result_dir, CFG.confidence_level, rep_tag, validate_tag);
            end ;

            %%% define outfile names
            fn_out_ir = sprintf('%s/%s_intron_retention%s_C%i.mat', CFG.result_dir, CFG.merge_strategy, rep_tag, CFG.confidence_level) ;
            fn_out_es = sprintf('%s/%s_exon_skip%s_C%i.mat', CFG.result_dir, CFG.merge_strategy, rep_tag, CFG.confidence_level) ;
            fn_out_mes = sprintf('%s/%s_mult_exon_skip%s_C%i.mat', CFG.result_dir, CFG.merge_strategy, rep_tag, CFG.confidence_level) ;
            fn_out_a5 = sprintf('%s/%s_alt_5prime%s_C%i.mat', CFG.result_dir, CFG.merge_strategy, rep_tag, CFG.confidence_level) ;
            fn_out_a3 = sprintf('%s/%s_alt_3prime%s_C%i.mat', CFG.result_dir, CFG.merge_strategy, rep_tag, CFG.confidence_level) ;

            intron_reten_pos{j, i} = [] ;
            exon_skip_pos{j, i} = [] ;
            mult_exon_skip_pos{j, i} = [] ;
            alt_end_5prime_pos{j, i} = [] ;
            alt_end_3prime_pos{j, i} = [] ;

            fprintf('\nconfidence %i / sample %i / replicate %i\n', CFG.confidence_level, i, j);

            if fexist(genes_fnames{j, i}),
                fprintf(1, 'Loading gene structure from %s ...\n', genes_fnames{j, i});
                L = load(genes_fnames{j, i}) ;
                fprintf(1, '... done.\n');

                %%% detect intron retentions from splicegraph
                if ~fexist(fn_out_ir),
                    [idx_intron_reten{j, i}, intron_intron_reten{j, i}] = detect_intronreten(L.genes, find([L.genes.is_alt])) ;
                    for k = 1:length(idx_intron_reten{j, i}),
                        gene = L.genes(idx_intron_reten{j, i}(k)) ;

                        %%% perform liftover between strains if necessary
                        exons = gene.splicegraph{1} ;
                        exons_col = convert_strain_pos_intervals(gene.chr_num, gene.splicegraph{1}', strain, reference_strain)' ;
                        exons_col_pos = convert_strain_pos(gene.chr_num, gene.splicegraph{1}', strain, reference_strain)' ;
                        if ~isequal(size(exons_col), size(exons_col_pos)), 
                            fprintf('skipping non-mappable intron retention event\n') ;
                            continue ;
                        end 

                        %%% build intron retention data structure
                        intron_reten_pos{j,i}(end+1).chr = gene.chr ;
                        intron_reten_pos{j,i}(end).chr_num = gene.chr_num ;
                        intron_reten_pos{j,i}(end).strand = gene.strand ;
                        intron_reten_pos{j,i}(end).strain = {strain} ;
                        intron_reten_pos{j,i}(end).exon1 = [exons(1,intron_intron_reten{j,i}(1,k)) ...
                                                            exons(2,intron_intron_reten{j,i}(1,k))] ;
                        intron_reten_pos{j,i}(end).exon1_col = [exons_col(1,intron_intron_reten{j,i}(1,k)) ...
                                                                exons_col(2,intron_intron_reten{j,i}(1,k))] ;
                        intron_reten_pos{j,i}(end).exon2 = [exons(1,intron_intron_reten{j,i}(2,k)) ...
                                                            exons(2,intron_intron_reten{j,i}(2,k))] ;
                        intron_reten_pos{j,i}(end).exon2_col = [exons_col(1,intron_intron_reten{j,i}(2,k)) ...
                                                                exons_col(2,intron_intron_reten{j,i}(2,k))] ;
                        intron_reten_pos{j,i}(end).intron = [exons(2,intron_intron_reten{j,i}(1,k))+1 ...
                                                             exons(1,intron_intron_reten{j,i}(2,k))-1] ;
                        intron_reten_pos{j,i}(end).intron_col = [exons_col(2,intron_intron_reten{j,i}(1,k))+1 ...
                                                                 exons_col(1,intron_intron_reten{j,i}(2,k))-1] ;
                        intron_reten_pos{j,i}(end).intron_col_pos = [exons_col_pos(2,intron_intron_reten{j,i}(1,k))+1 ...
                                                                     exons_col_pos(1,intron_intron_reten{j,i}(2,k))-1] ;
                        intron_reten_pos{j,i}(end).p_values = [] ;
                        if isfield(gene, 'name')
                            intron_reten_pos{j,i}(end).gene_name = {gene.name};
                        end;
                        if isfield(gene, 'transcript_type')
                            intron_reten_pos{j,i}(end).gene_type = {gene.transcript_type};
                        end;
                    end ;
                else
                    fprintf(1, '%s already exists\n', fn_out_ir);
                end;

                %%% detect exon_skips from splicegraph
                if ~fexist(fn_out_es),
                    [idx_exon_skip{j,i}, exon_exon_skip{j,i}] = detect_exonskips(L.genes, find([L.genes.is_alt])) ;
                    for k = 1:length(idx_exon_skip{j,i}),
                        gene = L.genes(idx_exon_skip{j,i}(k)) ;

                        %%% perform liftover between strains if necessary
                        exons = gene.splicegraph{1} ;
                        exons_col = convert_strain_pos_intervals(gene.chr_num, gene.splicegraph{1}', strain, reference_strain)' ;
                        exons_col_pos = convert_strain_pos(gene.chr_num, gene.splicegraph{1}', strain, reference_strain)' ;
                        if ~isequal(size(exons_col), size(exons_col_pos)), 
                            fprintf('skipping non-mappable exon skip event\n') ;
                            continue ;
                        end 

                        %%% build exon skip data structure
                        exon_skip_pos{j,i}(end+1).chr = gene.chr ;
                        exon_skip_pos{j,i}(end).chr_num = gene.chr_num ;
                        exon_skip_pos{j,i}(end).strand = gene.strand ;
                        exon_skip_pos{j,i}(end).strain = {strain} ;
                        exon_skip_pos{j,i}(end).exon_pre = [exons(1,exon_exon_skip{j,i}(1,k)) ...
                                                            exons(2,exon_exon_skip{j,i}(1,k))] ;
                        exon_skip_pos{j,i}(end).exon_pre_col = [exons_col(1,exon_exon_skip{j,i}(1,k)) ...
                                                                exons_col(2,exon_exon_skip{j,i}(1,k))] ;
                        exon_skip_pos{j,i}(end).exon = [exons(1,exon_exon_skip{j,i}(2,k)) ...
                                                        exons(2,exon_exon_skip{j,i}(2,k))] ;
                        exon_skip_pos{j,i}(end).exon_aft = [exons(1,exon_exon_skip{j,i}(3,k)) ...
                                                            exons(2,exon_exon_skip{j,i}(3,k))] ;
                        exon_skip_pos{j,i}(end).exon_aft_col = [exons_col(1,exon_exon_skip{j,i}(3,k)) ...
                                                                exons_col(2,exon_exon_skip{j,i}(3,k))] ;
                        exon_skip_pos{j,i}(end).exon_col = [exons_col(1,exon_exon_skip{j,i}(2,k)) ...
                                                            exons_col(2,exon_exon_skip{j,i}(2,k))] ;
                        exon_skip_pos{j,i}(end).exon_col_pos = [exons_col_pos(1,exon_exon_skip{j,i}(2,k)) ...
                                                                exons_col_pos(2,exon_exon_skip{j,i}(2,k))] ;
                        exon_skip_pos{j,i}(end).p_values=[] ;
                        if isfield(gene, 'name')
                            exon_skip_pos{j,i}(end).gene_name = {gene.name};
                        end;
                        if isfield(gene, 'transcript_type')
                            exon_skip_pos{j,i}(end).gene_type = {gene.transcript_type};
                        end;
                    end ;
                else
                    fprintf(1, '%s already exists\n', fn_out_es);
                end;

                %%% detect alternative intron_ends from splicegraph
                if ~fexist(fn_out_a5) || ~fexist(fn_out_a3) ,
                    [idx_alt_end_5prime{j,i}, exon_alt_end_5prime{j, i}, idx_alt_end_3prime{j, i}, exon_alt_end_3prime{j,i}] = detect_altprime(L.genes, find([L.genes.is_alt])) ;
                    %%% handle 5 prime events
                    for k = 1:length(idx_alt_end_5prime{j,i}),
                        gene = L.genes(idx_alt_end_5prime{j,i}(k));

                        %%% perform liftover between strains if necessary
                        exons = gene.splicegraph{1} ;
                        exons_col = convert_strain_pos_intervals(gene.chr_num, gene.splicegraph{1}', strain, reference_strain)' ;
                        exons_col_pos = convert_strain_pos(gene.chr_num, gene.splicegraph{1}', strain, reference_strain)' ;
                        if ~isequal(size(exons_col), size(exons_col_pos)), 
                            fprintf('skipping non-mappable alt 5 prime event\n') ;
                            continue ;
                        end;
                        
                        for k1 = 1:length(exon_alt_end_5prime{j, i}(k).fiveprimesites) - 1,
                            for k2 = k1 + 1:length(exon_alt_end_5prime{j, i}(k).fiveprimesites),

                                exon_alt1_col = exons_col(:, exon_alt_end_5prime{j, i}(k).fiveprimesites(k1))';
                                exon_alt2_col = exons_col(:, exon_alt_end_5prime{j, i}(k).fiveprimesites(k2))';

                                %%% check if exons overlap
                                if (exon_alt1_col(1) > exon_alt2_col(2) || exon_alt1_col(2) < exon_alt2_col(1))
                                    continue;
                                end;

                                alt_end_5prime_pos{j, i}(end + 1).chr = gene.chr ;
                                alt_end_5prime_pos{j, i}(end).chr_num = gene.chr_num ;
                                alt_end_5prime_pos{j, i}(end).strand = gene.strand ;
                                alt_end_5prime_pos{j, i}(end).strain = {strain};
                                alt_end_5prime_pos{j, i}(end).exon_const = exons(:, exon_alt_end_5prime{j, i}(k).threeprimesite(1))'; 
                                alt_end_5prime_pos{j, i}(end).exon_const_col = exons_col(:, exon_alt_end_5prime{j, i}(k).threeprimesite(1))'; 
                                alt_end_5prime_pos{j, i}(end).exon_alt1 = exons(:, exon_alt_end_5prime{j, i}(k).fiveprimesites(k1))'; 
                                alt_end_5prime_pos{j, i}(end).exon_alt1_col = exon_alt1_col;
                                alt_end_5prime_pos{j, i}(end).exon_alt2 = exons(:, exon_alt_end_5prime{j, i}(k).fiveprimesites(k2))'; 
                                alt_end_5prime_pos{j, i}(end).exon_alt2_col = exon_alt2_col;
                                if gene.strand == '+',
                                    alt_end_5prime_pos{j, i}(end).intron1 = [exons(2, exon_alt_end_5prime{j, i}(k).fiveprimesites(k1)) + 1 ...
                                                                             exons(1, exon_alt_end_5prime{j, i}(k).threeprimesite) - 1];
                                    alt_end_5prime_pos{j, i}(end).intron1_col = [exons_col(2, exon_alt_end_5prime{j, i}(k).fiveprimesites(k1)) + 1 ...
                                                                                 exons_col(1, exon_alt_end_5prime{j, i}(k).threeprimesite) - 1];
                                    alt_end_5prime_pos{j, i}(end).intron2 = [exons(2, exon_alt_end_5prime{j, i}(k).fiveprimesites(k2)) + 1 ...
                                                                             exons(1, exon_alt_end_5prime{j, i}(k).threeprimesite) - 1];
                                    alt_end_5prime_pos{j, i}(end).intron2_col = [exons_col(2, exon_alt_end_5prime{j, i}(k).fiveprimesites(k2)) + 1 ...
                                                                                 exons_col(1, exon_alt_end_5prime{j, i}(k).threeprimesite) - 1];
                                else
                                    alt_end_5prime_pos{j, i}(end).intron1 = [exons(2, exon_alt_end_5prime{j, i}(k).threeprimesite) + 1 ... 
                                                                             exons(1, exon_alt_end_5prime{j, i}(k).fiveprimesites(k1)) - 1];
                                    alt_end_5prime_pos{j, i}(end).intron1_col = [exons_col(2, exon_alt_end_5prime{j, i}(k).threeprimesite) + 1 ... 
                                                                                 exons_col(1, exon_alt_end_5prime{j, i}(k).fiveprimesites(k1)) - 1];
                                    alt_end_5prime_pos{j, i}(end).intron2 = [exons(2, exon_alt_end_5prime{j, i}(k).threeprimesite) + 1 ... 
                                                                             exons(1, exon_alt_end_5prime{j, i}(k).fiveprimesites(k2)) - 1];
                                    alt_end_5prime_pos{j, i}(end).intron2_col = [exons_col(2, exon_alt_end_5prime{j, i}(k).threeprimesite) + 1 ... 
                                                                                 exons_col(1, exon_alt_end_5prime{j, i}(k).fiveprimesites(k2)) - 1];
                                end;
                                alt_end_5prime_pos{j, i}(end).p_values = [];
                                if isfield(gene, 'name')
                                    alt_end_5prime_pos{j, i}(end).gene_name = {gene.name};
                                end;
                                if isfield(gene, 'transcript_type')
                                    alt_end_5prime_pos{j,i}(end).gene_type = {gene.transcript_type};
                                end;
                            end;
                        end;
                    end;
                    %%% handle 3 prime events
                    for k = 1:length(idx_alt_end_3prime{j,i}),
                        gene = L.genes(idx_alt_end_3prime{j,i}(k));

                        %%% perform liftover between strains if necessary
                        exons = gene.splicegraph{1} ;
                        exons_col = convert_strain_pos_intervals(gene.chr_num, gene.splicegraph{1}', strain, reference_strain)' ;
                        exons_col_pos = convert_strain_pos(gene.chr_num, gene.splicegraph{1}', strain, reference_strain)' ;
                        if ~isequal(size(exons_col), size(exons_col_pos)), 
                            fprintf('skipping non-mappable alt 3 prime event\n') ;
                            continue ;
                        end;

                        for k1 = 1:length(exon_alt_end_3prime{j, i}(k).threeprimesites) - 1,
                            for k2 = k1 + 1:length(exon_alt_end_3prime{j, i}(k).threeprimesites),
                                exon_alt1_col = exons_col(:, exon_alt_end_3prime{j, i}(k).threeprimesites(k1))';
                                exon_alt2_col = exons_col(:, exon_alt_end_3prime{j, i}(k).threeprimesites(k2))';

                                %%% check if exons overlap
                                if (exon_alt1_col(1) > exon_alt2_col(2) || exon_alt1_col(2) < exon_alt2_col(1))
                                    continue;
                                end;

                                alt_end_3prime_pos{j, i}(end + 1).chr = gene.chr ;
                                alt_end_3prime_pos{j, i}(end).chr_num = gene.chr_num ;
                                alt_end_3prime_pos{j, i}(end).strand = gene.strand ;
                                alt_end_3prime_pos{j, i}(end).strain = {strain};
                                alt_end_3prime_pos{j, i}(end).exon_const = exons(:, exon_alt_end_3prime{j, i}(k).fiveprimesite(1))'; 
                                alt_end_3prime_pos{j, i}(end).exon_const_col = exons_col(:, exon_alt_end_3prime{j, i}(k).fiveprimesite(1))'; 
                                alt_end_3prime_pos{j, i}(end).exon_alt1 = exons(:, exon_alt_end_3prime{j, i}(k).threeprimesites(k1))'; 
                                alt_end_3prime_pos{j, i}(end).exon_alt1_col = exon_alt1_col;
                                alt_end_3prime_pos{j, i}(end).exon_alt2 = exons(:, exon_alt_end_3prime{j, i}(k).threeprimesites(k2))'; 
                                alt_end_3prime_pos{j, i}(end).exon_alt2_col = exon_alt2_col; 
                                if gene.strand == '+',
                                    alt_end_3prime_pos{j, i}(end).intron1 = [exons(2, exon_alt_end_3prime{j, i}(k).fiveprimesite) + 1 ... 
                                                                             exons(1, exon_alt_end_3prime{j, i}(k).threeprimesites(k1)) - 1];
                                    alt_end_3prime_pos{j, i}(end).intron1_col = [exons_col(2, exon_alt_end_3prime{j, i}(k).fiveprimesite) + 1 ... 
                                                                                 exons_col(1, exon_alt_end_3prime{j, i}(k).threeprimesites(k1)) - 1];
                                    alt_end_3prime_pos{j, i}(end).intron2 = [exons(2, exon_alt_end_3prime{j, i}(k).fiveprimesite) + 1 ... 
                                                                             exons(1, exon_alt_end_3prime{j, i}(k).threeprimesites(k2)) - 1];
                                    alt_end_3prime_pos{j, i}(end).intron2_col = [exons_col(2, exon_alt_end_3prime{j, i}(k).fiveprimesite) + 1 ... 
                                                                                 exons_col(1, exon_alt_end_3prime{j, i}(k).threeprimesites(k2)) - 1];
                                else
                                    alt_end_3prime_pos{j, i}(end).intron1 = [exons(2, exon_alt_end_3prime{j, i}(k).threeprimesites(k1)) + 1, ...
                                                                             exons(1, exon_alt_end_3prime{j, i}(k).fiveprimesite) - 1];
                                    alt_end_3prime_pos{j, i}(end).intron1_col = [exons_col(2, exon_alt_end_3prime{j, i}(k).threeprimesites(k1)) + 1, ...
                                                                                 exons_col(1, exon_alt_end_3prime{j, i}(k).fiveprimesite) - 1];
                                    alt_end_3prime_pos{j, i}(end).intron2 = [exons(2, exon_alt_end_3prime{j, i}(k).threeprimesites(k2)) + 1, ...
                                                                             exons(1, exon_alt_end_3prime{j, i}(k).fiveprimesite) - 1];
                                    alt_end_3prime_pos{j, i}(end).intron2_col = [exons_col(2, exon_alt_end_3prime{j, i}(k).threeprimesites(k2)) + 1, ...
                                                                                 exons_col(1, exon_alt_end_3prime{j, i}(k).fiveprimesite) - 1];
                                end;
                                alt_end_3prime_pos{j, i}(end).p_values = [];
                                if isfield(gene, 'name')
                                    alt_end_3prime_pos{j, i}(end).gene_name = {gene.name};
                                end;
                                if isfield(gene, 'transcript_type')
                                    alt_end_3prime_pos{j,i}(end).gene_type = {gene.transcript_type};
                                end;
                            end;
                        end;
                    end;
                else
                    fprintf(1, '%s and %s already exists\n', fn_out_a5, fn_out_a3);
                end;

                %%% detect multiple_exon_skips from splicegraph
                if ~fexist(fn_out_mes),
                    [idx_mult_exon_skip{j,i}, exon_mult_exon_skip{j,i}, id_mult_exon_skip{j, i}] = detect_multipleskips(L.genes, find([L.genes.is_alt])) ;
                    %id_mult_exon_skip{j, i} = [];
                    for k = unique(id_mult_exon_skip{j, i}),
                    %for k = 1:length(idx_mult_exon_skip{j,i}),
                        k_ = find(id_mult_exon_skip{j, i} == k);
                        gene = L.genes(idx_mult_exon_skip{j,i}(k_(1))) ;

                        %%% perform liftover between strains if necessary
                        exons = gene.splicegraph{1} ;
                        exons_col = convert_strain_pos_intervals(gene.chr_num, gene.splicegraph{1}', strain, reference_strain)' ;
                        exons_col_pos = convert_strain_pos(gene.chr_num, gene.splicegraph{1}', strain, reference_strain)' ;
                        if ~isequal(size(exons_col), size(exons_col_pos)), 
                            fprintf('skipping non-mappable exon skip event\n') ;
                            continue ;
                        end 

                        %%% build multiple exon skip data structure
                        mult_exon_skip_pos{j,i}(end+1).chr = gene.chr ;
                        mult_exon_skip_pos{j,i}(end).chr_num = gene.chr_num ;
                        mult_exon_skip_pos{j,i}(end).strand = gene.strand ;
                        mult_exon_skip_pos{j,i}(end).strain = {strain} ;
                        mult_exon_skip_pos{j,i}(end).exon_pre = [exons(1,exon_mult_exon_skip{j,i}(k_(1))) ...
                                                                 exons(2,exon_mult_exon_skip{j,i}(k_(1)))] ;
                        mult_exon_skip_pos{j,i}(end).exon_pre_col = [exons_col(1,exon_mult_exon_skip{j,i}(k_(1))) ...
                                                                     exons_col(2,exon_mult_exon_skip{j,i}(k_(1)))] ;
                        mult_exon_skip_pos{j,i}(end).exons = [];
                        mult_exon_skip_pos{j,i}(end).exons_col = [];
                        mult_exon_skip_pos{j,i}(end).exons_col_pos = [];
                        for l = 2:(length(k_) - 1),
                            mult_exon_skip_pos{j,i}(end).exons = [ mult_exon_skip_pos{j,i}(end).exons, exons(1:2,exon_mult_exon_skip{j,i}(k_(l)))'] ;
                            mult_exon_skip_pos{j,i}(end).exons_col = [ mult_exon_skip_pos{j,i}(end).exons_col, exons_col(1:2,exon_mult_exon_skip{j,i}(k_(l)))'] ;
                            mult_exon_skip_pos{j,i}(end).exons_col_pos = [ mult_exon_skip_pos{j,i}(end).exons_col_pos, exons_col_pos(1:2,exon_mult_exon_skip{j,i}(k_(l)))'] ;
                        end;
                        mult_exon_skip_pos{j,i}(end).exon_aft = [exons(1,exon_mult_exon_skip{j,i}(k_(end))) ...
                                                                 exons(2,exon_mult_exon_skip{j,i}(k_(end)))] ;
                        mult_exon_skip_pos{j,i}(end).exon_aft_col = [exons_col(1,exon_mult_exon_skip{j,i}(k_(end))) ...
                                                                     exons_col(2,exon_mult_exon_skip{j,i}(k_(end)))] ;
                        mult_exon_skip_pos{j,i}(end).p_values=[] ;
                        if isfield(gene, 'name')
                            mult_exon_skip_pos{j,i}(end).gene_name = {gene.name};
                        end;
                        if isfield(gene, 'transcript_type')
                            mult_exon_skip_pos{j,i}(end).gene_type = {gene.transcript_type};
                        end;
                    end ;
                else
                    fprintf(1, '%s already exists\n', fn_out_mes);
                end;

            %%% genes file does not exist
            else
                fprintf('result file not found: %s\n', genes_fnames{j, i}) ;
                idx_intron_reten{j,i}=[] ;
                intron_intron_reten{j,i}=[] ;
                idx_exon_skip{j,i}=[] ;
                exon_exon_skip{j,i}=[] ;
                idx_mult_exon_skip{j,i}=[] ;
                exon_mult_exon_skip{j,i}=[] ;
            end  % if gene file exists;
        end ; % replicate_idxs
    end ; % samples 

    %%% combine events for all samples
    for ridx = replicate_idxs,

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMBINE INTRON RETENTIONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~fexist(fn_out_ir),
            intron_reten_pos_all = intron_reten_pos(replicate_idx,:) ;
            intron_reten_pos_all = [intron_reten_pos_all{:}] ;
            [intron_reten_pos_all(:).event_type] = deal('intron_retention');

            idx_valid_col = zeros(1,length(intron_reten_pos_all)) ;
            for i = 1:length(intron_reten_pos_all),
                idx_valid_col(i) = all(intron_reten_pos_all(i).intron_col_pos>0) ;
            end ;
            intron_reten_pos_all = intron_reten_pos_all(idx_valid_col~=0) ;
           
            %%% sort events by all coordinates
            intron_reten_pos_all = sort_events_full(intron_reten_pos_all) ;
           
            %%% make intron retentions unique by strain
            fprintf('\nMake intron retention events unique by strain\n');
            intron_reten_pos_all = make_unique_by_strain(intron_reten_pos_all);
            
            %%% sort events by event coordinates
            intron_reten_pos_all = sort_events_by_event(intron_reten_pos_all) ;
           
            %%% make intron retentions unique by event
            fprintf('\nMake intron retention events unique by event\n');
            intron_reten_pos_all = make_unique_by_event(intron_reten_pos_all);

            %%% count detected strains
            for i = 1:length(intron_reten_pos_all), 
                intron_reten_pos_all(i).num_detected = length(intron_reten_pos_all(i).strain) ;
            end ;
            
            %%% store intron retentions
            fprintf('saving intron retentions to %s\n', fn_out_ir) ;
            save(fn_out_ir, 'intron_reten_pos_all') ;
        else
            fprintf(1, '%s already exists\n', fn_out_ir);
        end;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMBINE EXON SKIPS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~fexist(fn_out_es),
            exon_skip_pos_all = exon_skip_pos(replicate_idx,:) ;
            exon_skip_pos_all = [exon_skip_pos_all{:}] ;
            [exon_skip_pos_all(:).event_type] = deal('exon_skip');

            idx_valid_col = zeros(1,length(exon_skip_pos_all)) ;
            for i = 1:length(exon_skip_pos_all),
                idx_valid_col(i) = all(exon_skip_pos_all(i).exon_col_pos>0) ;
            end ;
            exon_skip_pos_all = exon_skip_pos_all(idx_valid_col~=0) ;
           
            %%% sort exon skip events by all coordinates
            exon_skip_pos_all = sort_events_full(exon_skip_pos_all) ;
            
            %%% make exon skip events unique by strain
            fprintf('\nMake exon skip events unique by strain\n');
            exon_skip_pos_all = make_unique_by_strain(exon_skip_pos_all);

            %%% sort exon skip events by event coordinates
            exon_skip_pos_all = sort_events_by_event(exon_skip_pos_all) ;
            
            %%% make exon skip events unique by strain
            fprintf('\nMake exon skip events unique by event\n');
            exon_skip_pos_all = make_unique_by_event(exon_skip_pos_all);

            %%% count detected strains
            for i=1:length(exon_skip_pos_all),
                exon_skip_pos_all(i).num_detected = length(exon_skip_pos_all(i).strain) ;
            end ;
            
            %%% store exons skip events
            fprintf('saving exon skips to %s\n', fn_out_es) ;
            save(fn_out_es, 'exon_skip_pos_all') ;
        else
            fprintf(1, '%s already exists\n', fn_out_es);
        end;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMBINE MULTIPLE EXON SKIPS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~fexist(fn_out_mes),
            mult_exon_skip_pos_all = mult_exon_skip_pos(replicate_idx,:) ;
            mult_exon_skip_pos_all = [mult_exon_skip_pos_all{:}] ;
            [mult_exon_skip_pos_all(:).event_type] = deal('mult_exon_skip');

            if ~isempty([mult_exon_skip_pos_all.event_type]),
                idx_valid_col = zeros(1,length(mult_exon_skip_pos_all)) ;
                for i = 1:length(mult_exon_skip_pos_all),
                    idx_valid_col(i) = all(mult_exon_skip_pos_all(i).exons_col_pos(:) > 0) ;
                end ;
                mult_exon_skip_pos_all = mult_exon_skip_pos_all(idx_valid_col~=0) ;
               
                %%% sort exon skip events by all coordinates
                mult_exon_skip_pos_all = sort_events_full(mult_exon_skip_pos_all) ;
                
                %%% make multiple exon skip events unique by strain
                fprintf('\nMake multiple exon skip events unique by strain\n');
                mult_exon_skip_pos_all = make_unique_by_strain(mult_exon_skip_pos_all);

                %%% sort multiple exon skip events by event coordinates
                mult_exon_skip_pos_all = sort_events_by_event(mult_exon_skip_pos_all) ;
                
                %%% make multiple exon skip events unique by strain
                fprintf('\nMake multiple exon skip events unique by event\n');
                mult_exon_skip_pos_all = make_unique_by_event(mult_exon_skip_pos_all);

                %%% count detected strains
                for i=1:length(mult_exon_skip_pos_all),
                    mult_exon_skip_pos_all(i).num_detected = length(mult_exon_skip_pos_all(i).strain) ;
                end ;
            end ;
            
            %%% store exons skip events
            fprintf('saving multiple exon skips to %s\n', fn_out_mes) ;
            save(fn_out_mes, 'mult_exon_skip_pos_all') ;
        else
            fprintf(1, '%s already exists\n', fn_out_mes);
        end;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMBINE ALT FIVE PRIME EVENTS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~fexist(fn_out_a5),
            alt_end_5prime_pos_all = alt_end_5prime_pos(replicate_idx,:) ;
            alt_end_5prime_pos_all = [alt_end_5prime_pos_all{:}] ;
            [alt_end_5prime_pos_all(:).event_type] = deal('alt_5prime');
          
            %%% sort alt 5 prime events by all coordinates
            alt_end_5prime_pos_all = sort_events_full(alt_end_5prime_pos_all) ;
            
            %%% make alt 5 prime events unique by strain
            fprintf('\nMake alt 5 prime events unique by strain\n');
            alt_end_5prime_pos_all = make_unique_by_strain(alt_end_5prime_pos_all);

            %%% curate alt prime events
            %%% cut to min len, if alt exon lengths differ
            %%% remove, if alt exons do not overlap
            if CFG.curate_alt_prime_events == 1,
                alt_end_5prime_pos_all = curate_alt_prime(alt_end_5prime_pos_all);
            end

            %%% sort alt 5 prime events by event coordinates
            alt_end_5prime_pos_all = sort_events_by_event(alt_end_5prime_pos_all) ;
            
            %%% make alt 5 prime events unique by events
            fprintf('\nMake alt 5 prime events unique by event\n');
            alt_end_5prime_pos_all = make_unique_by_event(alt_end_5prime_pos_all);

            %%% count detected strains
            for i=1:length(alt_end_5prime_pos_all),
                alt_end_5prime_pos_all(i).num_detected = length(alt_end_5prime_pos_all(i).strain) ;
            end ;
            
            %%% store alt 5 prime events
            fprintf('saving alt 5 prime events to %s\n', fn_out_a5) ;
            save(fn_out_a5, 'alt_end_5prime_pos_all') ;
        else
            fprintf(1, '%s already exists\n', fn_out_a5);
        end;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMBINE ALT THREE PRIME EVENTS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~fexist(fn_out_a3),
            alt_end_3prime_pos_all = alt_end_3prime_pos(replicate_idx,:) ;
            alt_end_3prime_pos_all = [alt_end_3prime_pos_all{:}] ;
            [alt_end_3prime_pos_all(:).event_type] = deal('alt_3prime');
           
            %%% sort alt 3 prime events by all coordindates
            alt_end_3prime_pos_all = sort_events_full(alt_end_3prime_pos_all) ;
            
            %%% make alt 3 prime events unique by strain
            fprintf('\nMake alt 3 prime events unique\n');
            alt_end_3prime_pos_all = make_unique_by_strain(alt_end_3prime_pos_all);

            %%% curate alt prime events
            %%% cut to min len, if alt exon lengths differ
            %%% remove, if alt exons do not overlap
            if CFG.curate_alt_prime_events == 1,
                alt_end_3prime_pos_all = curate_alt_prime(alt_end_3prime_pos_all);
            end

            %%% sort alt 3 prime events by event coordinates
            alt_end_3prime_pos_all = sort_events_by_event(alt_end_3prime_pos_all) ;
            
            %%% make alt 3 prime events unique by event 
            fprintf('\nMake alt 3 prime events unique\n');
            alt_end_3prime_pos_all = make_unique_by_event(alt_end_3prime_pos_all);

            %%% count detected strains
            for i=1:length(alt_end_3prime_pos_all),
                alt_end_3prime_pos_all(i).num_detected = length(alt_end_3prime_pos_all(i).strain) ;
            end ;
            
            %%% store alt 3 prime events
            fprintf('saving alt 3 prime events to %s\n', fn_out_a3) ;
            save(fn_out_a3, 'alt_end_3prime_pos_all') ;
        else
            fprintf(1, '%s already exists\n', fn_out_a3);
        end;
    end ;
