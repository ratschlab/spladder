function alt_genes_collect(CFG)
% alt_genes_collect(CFG)

    %%% which events do we call
    do_exon_skip = 0;
    do_intron_retention = 0;
    do_mult_exon_skip = 0;
    do_alt_3prime = 0;
    do_alt_5prime = 0;
    for e_idx = 1:length(CFG.event_types),
        if strcmp('exon_skip', CFG.event_types{e_idx}),
            do_exon_skip = 1;
        elseif strcmp('intron_retention', CFG.event_types{e_idx}),
            do_intron_retention = 1;
        elseif strcmp('alt_3prime', CFG.event_types{e_idx}),
            do_alt_3prime = 1;
        elseif strcmp('alt_5prime', CFG.event_types{e_idx}),
            do_alt_5prime = 1;
        elseif strcmp('mult_exon_skip', CFG.event_types{e_idx}),
            do_mult_exon_skip = 1;
        else
            error(sprint('Unknown event type: %s\n', CFG.event_types{e_idx}));
        end;
    end;

    %%% init empty event fields
    if do_intron_retention, 
        idx_intron_reten = {} ;
        intron_intron_reten = {} ;
        intron_reten_pos = {} ;
    end;
    if do_exon_skip,
        idx_exon_skip = {} ;
        exon_exon_skip = {} ;
        exon_skip_pos = {} ;
    end;
    if do_alt_3prime || do_alt_5prime,
        idx_alt_end_5prime = {} ;
        exon_alt_end_5prime = {} ;
        alt_end_5prime_pos = {} ;
        idx_alt_end_3prime = {} ;
        exon_alt_end_3prime = {} ;
        alt_end_3prime_pos = {} ;
    end;
    if do_mult_exon_skip,
        idx_mult_exon_skip = {} ;
        exon_mult_exon_skip = {} ;
        id_mult_exon_skip = {} ;
        mult_exon_skip_pos = {};
    end;

    validate_tag = '';
    if isfield(CFG, 'validate_splicegraphs') && CFG.validate_splicegraphs,
        validate_tag = '.validated';
    end;

    %%% set offset that can be used for half open intervals
    if CFG.is_half_open,
        ho_offset = 1;
    else
        ho_offset = 0;
    end;

    for i = 1:length(CFG.samples),
        if CFG.same_genestruct_for_all_samples == 1 && i == 2,
            break;
        end;

        strain = CFG.strains{i} ;

        for ridx = CFG.replicate_idxs,
            if length(CFG.replicate_idxs) > 1,
                rep_tag = sprintf('_R%i', ridx);
            else
                rep_tag = '';
            end;

            if isfield(CFG, 'spladder_infile'),
                genes_fnames{ridx, i} = CFG.CFG.spladder_infile;
            elseif strcmp(CFG.merge_strategy, 'single')
                genes_fnames{ridx, i} = sprintf('%s/spladder/genes_graph_conf%i%s.%s.mat', CFG.out_dirname, CFG.confidence_level, rep_tag, CFG.samples{i});
            else
                genes_fnames{ridx, i} = sprintf('%s/spladder/genes_graph_conf%i%s.%s%s.mat', CFG.out_dirname, CFG.confidence_level, rep_tag, CFG.merge_strategy, validate_tag);
            end ;

            %%% define outfile names
            fn_out_ir = sprintf('%s/%s_intron_retention%s_C%i.mat', CFG.out_dirname, CFG.merge_strategy, rep_tag, CFG.confidence_level) ;
            fn_out_es = sprintf('%s/%s_exon_skip%s_C%i.mat', CFG.out_dirname, CFG.merge_strategy, rep_tag, CFG.confidence_level) ;
            fn_out_mes = sprintf('%s/%s_mult_exon_skip%s_C%i.mat', CFG.out_dirname, CFG.merge_strategy, rep_tag, CFG.confidence_level) ;
            fn_out_a5 = sprintf('%s/%s_alt_5prime%s_C%i.mat', CFG.out_dirname, CFG.merge_strategy, rep_tag, CFG.confidence_level) ;
            fn_out_a3 = sprintf('%s/%s_alt_3prime%s_C%i.mat', CFG.out_dirname, CFG.merge_strategy, rep_tag, CFG.confidence_level) ;

            intron_reten_pos{ridx, i} = [] ;
            exon_skip_pos{ridx, i} = [] ;
            mult_exon_skip_pos{ridx, i} = [] ;
            alt_end_5prime_pos{ridx, i} = [] ;
            alt_end_3prime_pos{ridx, i} = [] ;

            fprintf('\nconfidence %i / sample %i / replicate %i\n', CFG.confidence_level, i, ridx);

            if exist(genes_fnames{ridx, i}, 'file'),
                fprintf(1, 'Loading gene structure from %s ...\n', genes_fnames{ridx, i});
                L = load(genes_fnames{ridx, i}) ;
                fprintf(1, '... done.\n');

                %%% detect intron retentions from splicegraph
                if do_intron_retention,
                    if ~exist(fn_out_ir, 'file'),
                        [idx_intron_reten{ridx, i}, intron_intron_reten{ridx, i}] = detect_intronreten(L.genes, find([L.genes.is_alt])) ;
                        for k = 1:length(idx_intron_reten{ridx, i}),
                            gene = L.genes(idx_intron_reten{ridx, i}(k)) ;

                            %%% perform liftover between strains if necessary
                            exons = gene.splicegraph{1} ;
                            if ~isfield(CFG, 'reference_strain'),
                                exons_col = exons;
                                exons_col_pos = exons;
                            else
                                exons_col = convert_strain_pos_intervals(gene.chr_num, gene.splicegraph{1}', strain, CFG.reference_strain)' ;
                                exons_col_pos = convert_strain_pos(gene.chr_num, gene.splicegraph{1}', strain, CFG.reference_strain)' ;
                            end;
                            if ~isequal(size(exons_col), size(exons_col_pos)), 
                                fprintf('skipping non-mappable intron retention event\n') ;
                                continue ;
                            end 

                            %%% build intron retention data structure
                            intron_reten_pos{ridx,i}(end+1).chr = gene.chr ;
                            intron_reten_pos{ridx,i}(end).chr_num = gene.chr_num ;
                            intron_reten_pos{ridx,i}(end).strand = gene.strand ;
                            intron_reten_pos{ridx,i}(end).strain = {strain} ;
                            intron_reten_pos{ridx,i}(end).exon1 = [exons(1,intron_intron_reten{ridx,i}(1,k)) ...
                                                                exons(2,intron_intron_reten{ridx,i}(1,k))] ;
                            intron_reten_pos{ridx,i}(end).exon1_col = [exons_col(1,intron_intron_reten{ridx,i}(1,k)) ...
                                                                    exons_col(2,intron_intron_reten{ridx,i}(1,k))] ;
                            intron_reten_pos{ridx,i}(end).exon2 = [exons(1,intron_intron_reten{ridx,i}(2,k)) ...
                                                                exons(2,intron_intron_reten{ridx,i}(2,k))] ;
                            intron_reten_pos{ridx,i}(end).exon2_col = [exons_col(1,intron_intron_reten{ridx,i}(2,k)) ...
                                                                    exons_col(2,intron_intron_reten{ridx,i}(2,k))] ;
                            intron_reten_pos{ridx,i}(end).intron = [exons(2,intron_intron_reten{ridx,i}(1,k)) + 1 - ho_offset ...
                                                                 exons(1,intron_intron_reten{ridx,i}(2,k)) - 1 + ho_offset] ;
                            intron_reten_pos{ridx,i}(end).intron_col = [exons_col(2,intron_intron_reten{ridx,i}(1,k)) + 1 - ho_offset ...
                                                                     exons_col(1,intron_intron_reten{ridx,i}(2,k)) - 1 + ho_offset] ;
                            intron_reten_pos{ridx,i}(end).intron_col_pos = [exons_col_pos(2,intron_intron_reten{ridx,i}(1,k)) + 1 - ho_offset ...
                                                                         exons_col_pos(1,intron_intron_reten{ridx,i}(2,k)) - 1 + ho_offset] ;
                            intron_reten_pos{ridx,i}(end).p_values = [] ;
                            if isfield(gene, 'name')
                                intron_reten_pos{ridx,i}(end).gene_name = {gene.name};
                            end;
                            if isfield(gene, 'transcript_type')
                                intron_reten_pos{ridx,i}(end).gene_type = {gene.transcript_type};
                            end;
                        end ;
                    else
                        fprintf(1, '%s already exists\n', fn_out_ir);
                    end;
                end;

                %%% detect exon_skips from splicegraph
                if do_exon_skip,
                    if ~exist(fn_out_es, 'file'),
                        [idx_exon_skip{ridx,i}, exon_exon_skip{ridx,i}] = detect_exonskips(L.genes, find([L.genes.is_alt])) ;
                        for k = 1:length(idx_exon_skip{ridx,i}),
                            gene = L.genes(idx_exon_skip{ridx,i}(k)) ;

                            %%% perform liftover between strains if necessary
                            exons = gene.splicegraph{1} ;
                            if ~isfield(CFG, 'reference_strain'),
                                exons_col = exons;
                                exons_col_pos = exons;
                            else
                                exons_col = convert_strain_pos_intervals(gene.chr_num, gene.splicegraph{1}', strain, CFG.reference_strain)' ;
                                exons_col_pos = convert_strain_pos(gene.chr_num, gene.splicegraph{1}', strain, CFG.reference_strain)' ;
                            end;
                            if ~isequal(size(exons_col), size(exons_col_pos)), 
                                fprintf('skipping non-mappable exon skip event\n') ;
                                continue ;
                            end 

                            %%% build exon skip data structure
                            exon_skip_pos{ridx,i}(end+1).chr = gene.chr ;
                            exon_skip_pos{ridx,i}(end).chr_num = gene.chr_num ;
                            exon_skip_pos{ridx,i}(end).strand = gene.strand ;
                            exon_skip_pos{ridx,i}(end).strain = {strain} ;
                            exon_skip_pos{ridx,i}(end).exon_pre = [exons(1,exon_exon_skip{ridx,i}(1,k)) ...
                                                                exons(2,exon_exon_skip{ridx,i}(1,k))] ;
                            exon_skip_pos{ridx,i}(end).exon_pre_col = [exons_col(1,exon_exon_skip{ridx,i}(1,k)) ...
                                                                    exons_col(2,exon_exon_skip{ridx,i}(1,k))] ;
                            exon_skip_pos{ridx,i}(end).exon = [exons(1,exon_exon_skip{ridx,i}(2,k)) ...
                                                            exons(2,exon_exon_skip{ridx,i}(2,k))] ;
                            exon_skip_pos{ridx,i}(end).exon_aft = [exons(1,exon_exon_skip{ridx,i}(3,k)) ...
                                                                exons(2,exon_exon_skip{ridx,i}(3,k))] ;
                            exon_skip_pos{ridx,i}(end).exon_aft_col = [exons_col(1,exon_exon_skip{ridx,i}(3,k)) ...
                                                                    exons_col(2,exon_exon_skip{ridx,i}(3,k))] ;
                            exon_skip_pos{ridx,i}(end).exon_col = [exons_col(1,exon_exon_skip{ridx,i}(2,k)) ...
                                                                exons_col(2,exon_exon_skip{ridx,i}(2,k))] ;
                            exon_skip_pos{ridx,i}(end).exon_col_pos = [exons_col_pos(1,exon_exon_skip{ridx,i}(2,k)) ...
                                                                    exons_col_pos(2,exon_exon_skip{ridx,i}(2,k))] ;
                            exon_skip_pos{ridx,i}(end).p_values=[] ;
                            if isfield(gene, 'name')
                                exon_skip_pos{ridx,i}(end).gene_name = {gene.name};
                            end;
                            if isfield(gene, 'transcript_type')
                                exon_skip_pos{ridx,i}(end).gene_type = {gene.transcript_type};
                            end;
                        end ;
                    else
                        fprintf(1, '%s already exists\n', fn_out_es);
                    end;
                end;

                %%% detect alternative intron_ends from splicegraph
                if do_alt_3prime || do_alt_5prime, 
                    if ~exist(fn_out_a5, 'file') || ~exist(fn_out_a3, 'file'),
                        [idx_alt_end_5prime{ridx,i}, exon_alt_end_5prime{ridx, i}, idx_alt_end_3prime{ridx, i}, exon_alt_end_3prime{ridx,i}] = detect_altprime(L.genes, find([L.genes.is_alt])) ;
                        %%% handle 5 prime events
                        for k = 1:length(idx_alt_end_5prime{ridx,i}),
                            gene = L.genes(idx_alt_end_5prime{ridx,i}(k));

                            %%% perform liftover between strains if necessary
                            exons = gene.splicegraph{1} ;
                            if ~isfield(CFG, 'reference_strain'),
                                exons_col = exons;
                                exons_col_pos = exons;
                            else
                                exons_col = convert_strain_pos_intervals(gene.chr_num, gene.splicegraph{1}', strain, CFG.reference_strain)' ;
                                exons_col_pos = convert_strain_pos(gene.chr_num, gene.splicegraph{1}', strain, CFG.reference_strain)' ;
                            end;
                            if ~isequal(size(exons_col), size(exons_col_pos)), 
                                fprintf('skipping non-mappable alt 5 prime event\n') ;
                                continue ;
                            end;
                            
                            for k1 = 1:length(exon_alt_end_5prime{ridx, i}(k).fiveprimesites) - 1,
                                for k2 = k1 + 1:length(exon_alt_end_5prime{ridx, i}(k).fiveprimesites),

                                    exon_alt1_col = exons_col(:, exon_alt_end_5prime{ridx, i}(k).fiveprimesites(k1))';
                                    exon_alt2_col = exons_col(:, exon_alt_end_5prime{ridx, i}(k).fiveprimesites(k2))';

                                    %%% check if exons overlap
                                    if (exon_alt1_col(1) > exon_alt2_col(2) || exon_alt1_col(2) < exon_alt2_col(1))
                                        continue;
                                    end;

                                    alt_end_5prime_pos{ridx, i}(end + 1).chr = gene.chr ;
                                    alt_end_5prime_pos{ridx, i}(end).chr_num = gene.chr_num ;
                                    alt_end_5prime_pos{ridx, i}(end).strand = gene.strand ;
                                    alt_end_5prime_pos{ridx, i}(end).strain = {strain};
                                    alt_end_5prime_pos{ridx, i}(end).exon_const = exons(:, exon_alt_end_5prime{ridx, i}(k).threeprimesite(1))'; 
                                    alt_end_5prime_pos{ridx, i}(end).exon_const_col = exons_col(:, exon_alt_end_5prime{ridx, i}(k).threeprimesite(1))'; 
                                    alt_end_5prime_pos{ridx, i}(end).exon_alt1 = exons(:, exon_alt_end_5prime{ridx, i}(k).fiveprimesites(k1))'; 
                                    alt_end_5prime_pos{ridx, i}(end).exon_alt1_col = exon_alt1_col;
                                    alt_end_5prime_pos{ridx, i}(end).exon_alt2 = exons(:, exon_alt_end_5prime{ridx, i}(k).fiveprimesites(k2))'; 
                                    alt_end_5prime_pos{ridx, i}(end).exon_alt2_col = exon_alt2_col;
                                    if gene.strand == '+',
                                        alt_end_5prime_pos{ridx, i}(end).intron1 = [exons(2, exon_alt_end_5prime{ridx, i}(k).fiveprimesites(k1)) + 1 - ho_offset, ...
                                                                                    exons(1, exon_alt_end_5prime{ridx, i}(k).threeprimesite) - 1 + ho_offset];
                                        alt_end_5prime_pos{ridx, i}(end).intron1_col = [exons_col(2, exon_alt_end_5prime{ridx, i}(k).fiveprimesites(k1)) + 1 - ho_offset, ...
                                                                                        exons_col(1, exon_alt_end_5prime{ridx, i}(k).threeprimesite) - 1 + ho_offset];
                                        alt_end_5prime_pos{ridx, i}(end).intron2 = [exons(2, exon_alt_end_5prime{ridx, i}(k).fiveprimesites(k2)) + 1 - ho_offset, ...
                                                                                    exons(1, exon_alt_end_5prime{ridx, i}(k).threeprimesite) - 1 + ho_offset];
                                        alt_end_5prime_pos{ridx, i}(end).intron2_col = [exons_col(2, exon_alt_end_5prime{ridx, i}(k).fiveprimesites(k2)) + 1 - ho_offset, ...
                                                                                        exons_col(1, exon_alt_end_5prime{ridx, i}(k).threeprimesite) - 1 + ho_offset];
                                    else
                                        alt_end_5prime_pos{ridx, i}(end).intron1 = [exons(2, exon_alt_end_5prime{ridx, i}(k).threeprimesite) + 1 - ho_offset, ... 
                                                                                    exons(1, exon_alt_end_5prime{ridx, i}(k).fiveprimesites(k1)) - 1 + ho_offset];
                                        alt_end_5prime_pos{ridx, i}(end).intron1_col = [exons_col(2, exon_alt_end_5prime{ridx, i}(k).threeprimesite) + 1 - ho_offset, ... 
                                                                                        exons_col(1, exon_alt_end_5prime{ridx, i}(k).fiveprimesites(k1)) - 1 + ho_offset];
                                        alt_end_5prime_pos{ridx, i}(end).intron2 = [exons(2, exon_alt_end_5prime{ridx, i}(k).threeprimesite) + 1 - ho_offset, ... 
                                                                                    exons(1, exon_alt_end_5prime{ridx, i}(k).fiveprimesites(k2)) - 1 + ho_offset];
                                        alt_end_5prime_pos{ridx, i}(end).intron2_col = [exons_col(2, exon_alt_end_5prime{ridx, i}(k).threeprimesite) + 1 - ho_offset, ... 
                                                                                        exons_col(1, exon_alt_end_5prime{ridx, i}(k).fiveprimesites(k2)) - 1 + ho_offset];
                                    end;
                                    alt_end_5prime_pos{ridx, i}(end).p_values = [];
                                    if isfield(gene, 'name')
                                        alt_end_5prime_pos{ridx, i}(end).gene_name = {gene.name};
                                    end;
                                    if isfield(gene, 'transcript_type')
                                        alt_end_5prime_pos{ridx,i}(end).gene_type = {gene.transcript_type};
                                    end;
                                end;
                            end;
                        end;
                        %%% handle 3 prime events
                        for k = 1:length(idx_alt_end_3prime{ridx,i}),
                            gene = L.genes(idx_alt_end_3prime{ridx,i}(k));

                            %%% perform liftover between strains if necessary
                            exons = gene.splicegraph{1} ;
                            if ~isfield(CFG, 'reference_strain'),
                                exons_col = exons;
                                exons_col_pos = exons;
                            else
                                exons_col = convert_strain_pos_intervals(gene.chr_num, gene.splicegraph{1}', strain, CFG.reference_strain)' ;
                                exons_col_pos = convert_strain_pos(gene.chr_num, gene.splicegraph{1}', strain, CFG.reference_strain)' ;
                            end;
                            if ~isequal(size(exons_col), size(exons_col_pos)), 
                                fprintf('skipping non-mappable alt 3 prime event\n') ;
                                continue ;
                            end;

                            for k1 = 1:length(exon_alt_end_3prime{ridx, i}(k).threeprimesites) - 1,
                                for k2 = k1 + 1:length(exon_alt_end_3prime{ridx, i}(k).threeprimesites),
                                    exon_alt1_col = exons_col(:, exon_alt_end_3prime{ridx, i}(k).threeprimesites(k1))';
                                    exon_alt2_col = exons_col(:, exon_alt_end_3prime{ridx, i}(k).threeprimesites(k2))';

                                    %%% check if exons overlap
                                    if (exon_alt1_col(1) > exon_alt2_col(2) || exon_alt1_col(2) < exon_alt2_col(1))
                                        continue;
                                    end;

                                    alt_end_3prime_pos{ridx, i}(end + 1).chr = gene.chr ;
                                    alt_end_3prime_pos{ridx, i}(end).chr_num = gene.chr_num ;
                                    alt_end_3prime_pos{ridx, i}(end).strand = gene.strand ;
                                    alt_end_3prime_pos{ridx, i}(end).strain = {strain};
                                    alt_end_3prime_pos{ridx, i}(end).exon_const = exons(:, exon_alt_end_3prime{ridx, i}(k).fiveprimesite(1))'; 
                                    alt_end_3prime_pos{ridx, i}(end).exon_const_col = exons_col(:, exon_alt_end_3prime{ridx, i}(k).fiveprimesite(1))'; 
                                    alt_end_3prime_pos{ridx, i}(end).exon_alt1 = exons(:, exon_alt_end_3prime{ridx, i}(k).threeprimesites(k1))'; 
                                    alt_end_3prime_pos{ridx, i}(end).exon_alt1_col = exon_alt1_col;
                                    alt_end_3prime_pos{ridx, i}(end).exon_alt2 = exons(:, exon_alt_end_3prime{ridx, i}(k).threeprimesites(k2))'; 
                                    alt_end_3prime_pos{ridx, i}(end).exon_alt2_col = exon_alt2_col; 
                                    if gene.strand == '+',
                                        alt_end_3prime_pos{ridx, i}(end).intron1 = [exons(2, exon_alt_end_3prime{ridx, i}(k).fiveprimesite) + 1 - ho_offset, ... 
                                                                                    exons(1, exon_alt_end_3prime{ridx, i}(k).threeprimesites(k1)) - 1 + ho_offset];
                                        alt_end_3prime_pos{ridx, i}(end).intron1_col = [exons_col(2, exon_alt_end_3prime{ridx, i}(k).fiveprimesite) + 1 - ho_offset, ... 
                                                                                        exons_col(1, exon_alt_end_3prime{ridx, i}(k).threeprimesites(k1)) - 1 + ho_offset];
                                        alt_end_3prime_pos{ridx, i}(end).intron2 = [exons(2, exon_alt_end_3prime{ridx, i}(k).fiveprimesite) + 1 - ho_offset, ... 
                                                                                    exons(1, exon_alt_end_3prime{ridx, i}(k).threeprimesites(k2)) - 1 + ho_offset];
                                        alt_end_3prime_pos{ridx, i}(end).intron2_col = [exons_col(2, exon_alt_end_3prime{ridx, i}(k).fiveprimesite) + 1 - ho_offset, ... 
                                                                                        exons_col(1, exon_alt_end_3prime{ridx, i}(k).threeprimesites(k2)) - 1 + ho_offset];
                                    else
                                        alt_end_3prime_pos{ridx, i}(end).intron1 = [exons(2, exon_alt_end_3prime{ridx, i}(k).threeprimesites(k1)) + 1 - ho_offset, ...
                                                                                    exons(1, exon_alt_end_3prime{ridx, i}(k).fiveprimesite) - 1 + ho_offset];
                                        alt_end_3prime_pos{ridx, i}(end).intron1_col = [exons_col(2, exon_alt_end_3prime{ridx, i}(k).threeprimesites(k1)) + 1 - ho_offset, ...
                                                                                        exons_col(1, exon_alt_end_3prime{ridx, i}(k).fiveprimesite) - 1 + ho_offset];
                                        alt_end_3prime_pos{ridx, i}(end).intron2 = [exons(2, exon_alt_end_3prime{ridx, i}(k).threeprimesites(k2)) + 1 - ho_offset, ...
                                                                                    exons(1, exon_alt_end_3prime{ridx, i}(k).fiveprimesite) - 1 + ho_offset];
                                        alt_end_3prime_pos{ridx, i}(end).intron2_col = [exons_col(2, exon_alt_end_3prime{ridx, i}(k).threeprimesites(k2)) + 1 - ho_offset, ...
                                                                                        exons_col(1, exon_alt_end_3prime{ridx, i}(k).fiveprimesite) - 1 + ho_offset];
                                    end;
                                    alt_end_3prime_pos{ridx, i}(end).p_values = [];
                                    if isfield(gene, 'name')
                                        alt_end_3prime_pos{ridx, i}(end).gene_name = {gene.name};
                                    end;
                                    if isfield(gene, 'transcript_type')
                                        alt_end_3prime_pos{ridx,i}(end).gene_type = {gene.transcript_type};
                                    end;
                                end;
                            end;
                        end;
                    else
                        fprintf(1, '%s and %s already exists\n', fn_out_a5, fn_out_a3);
                    end;
                end;

                %%% detect multiple_exon_skips from splicegraph
                if do_mult_exon_skip,
                    if ~exist(fn_out_mes, 'file'),
                        [idx_mult_exon_skip{ridx,i}, exon_mult_exon_skip{ridx,i}, id_mult_exon_skip{ridx, i}] = detect_multipleskips(L.genes, find([L.genes.is_alt])) ;
                        if ~isempty(id_mult_exon_skip{ridx, i}), % necessary for Octave
                            for k = unique(id_mult_exon_skip{ridx, i}),
                                k_ = find(id_mult_exon_skip{ridx, i} == k);
                                gene = L.genes(idx_mult_exon_skip{ridx,i}(k_(1))) ;

                                %%% perform liftover between strains if necessary
                                exons = gene.splicegraph{1} ;
                                if ~isfield(CFG, 'reference_strain'),
                                    exons_col = exons;
                                    exons_col_pos = exons;
                                else
                                    exons_col = convert_strain_pos_intervals(gene.chr_num, gene.splicegraph{1}', strain, CFG.reference_strain)' ;
                                    exons_col_pos = convert_strain_pos(gene.chr_num, gene.splicegraph{1}', strain, CFG.reference_strain)' ;
                                end;
                                if ~isequal(size(exons_col), size(exons_col_pos)), 
                                    fprintf('skipping non-mappable exon skip event\n') ;
                                    continue ;
                                end 

                                %%% build multiple exon skip data structure
                                mult_exon_skip_pos{ridx,i}(end+1).chr = gene.chr ;
                                mult_exon_skip_pos{ridx,i}(end).chr_num = gene.chr_num ;
                                mult_exon_skip_pos{ridx,i}(end).strand = gene.strand ;
                                mult_exon_skip_pos{ridx,i}(end).strain = {strain} ;
                                mult_exon_skip_pos{ridx,i}(end).exon_pre = [exons(1,exon_mult_exon_skip{ridx,i}(k_(1))) ...
                                                                         exons(2,exon_mult_exon_skip{ridx,i}(k_(1)))] ;
                                mult_exon_skip_pos{ridx,i}(end).exon_pre_col = [exons_col(1,exon_mult_exon_skip{ridx,i}(k_(1))) ...
                                                                             exons_col(2,exon_mult_exon_skip{ridx,i}(k_(1)))] ;
                                mult_exon_skip_pos{ridx,i}(end).exons = [];
                                mult_exon_skip_pos{ridx,i}(end).exons_col = [];
                                mult_exon_skip_pos{ridx,i}(end).exons_col_pos = [];
                                for l = 2:(length(k_) - 1),
                                    mult_exon_skip_pos{ridx,i}(end).exons = [ mult_exon_skip_pos{ridx,i}(end).exons, exons(1:2,exon_mult_exon_skip{ridx,i}(k_(l)))'] ;
                                    mult_exon_skip_pos{ridx,i}(end).exons_col = [ mult_exon_skip_pos{ridx,i}(end).exons_col, exons_col(1:2,exon_mult_exon_skip{ridx,i}(k_(l)))'] ;
                                    mult_exon_skip_pos{ridx,i}(end).exons_col_pos = [ mult_exon_skip_pos{ridx,i}(end).exons_col_pos, exons_col_pos(1:2,exon_mult_exon_skip{ridx,i}(k_(l)))'] ;
                                end;
                                mult_exon_skip_pos{ridx,i}(end).exon_aft = [exons(1,exon_mult_exon_skip{ridx,i}(k_(end))) ...
                                                                         exons(2,exon_mult_exon_skip{ridx,i}(k_(end)))] ;
                                mult_exon_skip_pos{ridx,i}(end).exon_aft_col = [exons_col(1,exon_mult_exon_skip{ridx,i}(k_(end))) ...
                                                                             exons_col(2,exon_mult_exon_skip{ridx,i}(k_(end)))] ;
                                mult_exon_skip_pos{ridx,i}(end).p_values=[] ;
                                if isfield(gene, 'name')
                                    mult_exon_skip_pos{ridx,i}(end).gene_name = {gene.name};
                                end;
                                if isfield(gene, 'transcript_type')
                                    mult_exon_skip_pos{ridx,i}(end).gene_type = {gene.transcript_type};
                                end;
                            end ;
                        end ;
                    else
                        fprintf(1, '%s already exists\n', fn_out_mes);
                    end;
                end;
            %%% genes file does not exist
            else
                fprintf('result file not found: %s\n', genes_fnames{ridx, i}) ;
                idx_intron_reten{ridx,i}=[] ;
                intron_intron_reten{ridx,i}=[] ;
                idx_exon_skip{ridx,i}=[] ;
                exon_exon_skip{ridx,i}=[] ;
                idx_mult_exon_skip{ridx,i}=[] ;
                exon_mult_exon_skip{ridx,i}=[] ;
            end  % if gene file exists;
        end ; % replicate_idxs
    end ; % samples 

    %%% combine events for all samples
    for ridx = CFG.replicate_idxs,

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMBINE INTRON RETENTIONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if do_intron_retention,
            if ~exist(fn_out_ir, 'file'),
                intron_reten_pos_all = intron_reten_pos(ridx,:) ;
                intron_reten_pos_all = [intron_reten_pos_all{:}] ;
                [intron_reten_pos_all(:).event_type] = deal('intron_retention');

                %%% post process event structure by sorting and making events unique
                events_all = post_process_event_struct(intron_reten_pos_all);

                %%% store intron retentions
                fprintf('saving intron retentions to %s\n', fn_out_ir) ;
                save(fn_out_ir, 'events_all') ;
            else
                fprintf(1, '%s already exists\n', fn_out_ir);
            end;
        end;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMBINE EXON SKIPS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if do_exon_skip,
            if ~exist(fn_out_es, 'file'),
                exon_skip_pos_all = exon_skip_pos(ridx,:) ;
                exon_skip_pos_all = [exon_skip_pos_all{:}] ;
                [exon_skip_pos_all(:).event_type] = deal('exon_skip');

                %%% post process event structure by sorting and making events unique
                events_all = post_process_event_struct(exon_skip_pos_all);

                %%% store exons skip events
                fprintf('saving exon skips to %s\n', fn_out_es) ;
                save(fn_out_es, 'events_all') ;
            else
                fprintf(1, '%s already exists\n', fn_out_es);
            end;
        end;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMBINE MULTIPLE EXON SKIPS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if do_mult_exon_skip,
            if ~exist(fn_out_mes, 'file'),
                mult_exon_skip_pos_all = mult_exon_skip_pos(ridx,:) ;
                mult_exon_skip_pos_all = [mult_exon_skip_pos_all{:}] ;
                [mult_exon_skip_pos_all(:).event_type] = deal('mult_exon_skip');

                if ~isempty([mult_exon_skip_pos_all.event_type]),
                    %%% post process event structure by sorting and making events unique
                    events_all = post_process_event_struct(mult_exon_skip_pos_all);
                else
                    events_all = mult_exon_skip_pos_all;
                end ;
                
                %%% store exons skip events
                fprintf('saving multiple exon skips to %s\n', fn_out_mes) ;
                save(fn_out_mes, 'events_all') ;
            else
                fprintf(1, '%s already exists\n', fn_out_mes);
            end;
        end;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMBINE ALT FIVE PRIME EVENTS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if do_alt_5prime, 
            if ~exist(fn_out_a5, 'file'),
                alt_end_5prime_pos_all = alt_end_5prime_pos(ridx,:) ;
                alt_end_5prime_pos_all = [alt_end_5prime_pos_all{:}] ;
                [alt_end_5prime_pos_all(:).event_type] = deal('alt_5prime');
              
                %%% post process event structure by sorting and making events unique
                events_all = post_process_event_struct(alt_end_5prime_pos_all);

                %%% curate alt prime events
                %%% cut to min len, if alt exon lengths differ
                %%% remove, if alt exons do not overlap
                if CFG.curate_alt_prime_events == 1,
                    events_all = curate_alt_prime(events_all);
                end

                %%% store alt 5 prime events
                fprintf('saving alt 5 prime events to %s\n', fn_out_a5) ;
                save(fn_out_a5, 'events_all') ;
            else
                fprintf(1, '%s already exists\n', fn_out_a5);
            end;
        end;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMBINE ALT THREE PRIME EVENTS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if do_alt_3prime, 
            if ~exist(fn_out_a3, 'file'),
                alt_end_3prime_pos_all = alt_end_3prime_pos(ridx,:) ;
                alt_end_3prime_pos_all = [alt_end_3prime_pos_all{:}] ;
                [alt_end_3prime_pos_all(:).event_type] = deal('alt_3prime');
               
                %%% post process event structure by sorting and making events unique
                events_all = post_process_event_struct(alt_end_3prime_pos_all);

                %%% curate alt prime events
                %%% cut to min len, if alt exon lengths differ
                %%% remove, if alt exons do not overlap
                if CFG.curate_alt_prime_events == 1,
                    events_all = curate_alt_prime(events_all);
                end

                %%% store alt 3 prime events
                fprintf('saving alt 3 prime events to %s\n', fn_out_a3) ;
                save(fn_out_a3, 'events_all') ;
            else
                fprintf(1, '%s already exists\n', fn_out_a3);
            end;
        end;
    end ;
