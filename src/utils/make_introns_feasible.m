function genes = make_introns_feasible(genes, fn_bam, CFG)
% genes = make_introns_feasible(genes, fn_bam, CFG)

    tmp1 = [];
    tmp2 = [];
    for i = 1:length(genes),
        tmp1 = [tmp1 length(genes(i).introns{1})];
        tmp2 = [tmp2 length(genes(i).introns{2})];
    end;
    
    unfeas = find(tmp1 > 200 | tmp2 > 200);
    fprintf(CFG.fd_log, 'found %i unfeasible genes\n', length(unfeas));

    while ~isempty(unfeas),
        %%% make filter more stringent
        CFG.intron_filter.exon_len = min(36, CFG.intron_filter.exon_len + 4);
        CFG.intron_filter.mincount = CFG.intron_filter.mincount + 10;
        CFG.intron_filter.mismatch = max(CFG.intron_filter.mismatch - 1, 0);

        %%% get new intron counts
        tmp_genes = append_intron_list(genes(unfeas), fn_bam, CFG);
        genes(unfeas) = tmp_genes;

        %%% stil unfeasible?
        tmp1 = [];
        tmp2 = [];
        for i = 1:length(genes),
            tmp1 = [tmp1 length(genes(i).introns{1})];
            tmp2 = [tmp2 length(genes(i).introns{2})];
        end;
        still_unfeas = find(tmp1 > 200 | tmp2 > 200);
        for i = setdiff(unfeas, still_unfeas),
            fprintf(CFG.fd_log, '[feasibility] set criteria for gene %s to: min_ex %i, min_conf %i, max_mism %i\n', genes(i).name, CFG.intron_filter.exon_len, CFG.intron_filter.mincount, CFG.intron_filter.mismatch);
        end;
        unfeas = still_unfeas;
    end;

