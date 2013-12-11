function introns = make_introns_feasible(introns, genes, CFG)
% introns = make_introns_feasible(introns, genes, CFG)

    tmp1 = [];
    tmp2 = [];
    for i = 1:length(introns),
        tmp1 = [tmp1 length(introns{i, 1})];
        tmp2 = [tmp2 length(introns{i, 2})];
    end;
    
    unfeas = find(tmp1 > 200 | tmp2 > 200);
    fprintf(CFG.fd_log, 'found %i unfeasible genes\n', length(unfeas));

    while ~isempty(unfeas),
        %%% make filter more stringent
        CFG.read_filter.exon_len = min(36, CFG.read_filter.exon_len + 4);
        CFG.read_filter.mincount = CFG.read_filter.mincount + 10;
        CFG.read_filter.mismatch = max(CFG.read_filter.mismatch - 1, 0);

        %%% get new intron counts
        tmp_introns = get_intron_list(genes(unfeas), CFG);
        introns(unfeas, :) = tmp_introns;

        %%% stil unfeasible?
        tmp1 = [];
        tmp2 = [];
        for i = 1:length(introns),
            tmp1 = [tmp1 length(introns{i, 1})];
            tmp2 = [tmp2 length(introns{i, 2})];
        end;
        still_unfeas = find(tmp1 > 200 | tmp2 > 200);
        for i = setdiff(unfeas, still_unfeas),
            fprintf(CFG.fd_log, '[feasibility] set criteria for gene %s to: min_ex %i, min_conf %i, max_mism %i\n', genes(i).name, CFG.read_filter.exon_len, CFG.read_filter.mincount, CFG.read_filter.mismatch);
        end;
        unfeas = still_unfeas;
    end;

