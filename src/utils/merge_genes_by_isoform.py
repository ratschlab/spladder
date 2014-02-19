def merge_genes_by_isoform(CFG):
    ### This script takes several gene structures and merges them into one. If the number of isoforms of an annotated transcript should 
    ### exceed max_num_isoforms, max_num_isoforms many are subsampled
    
    appended = False
    max_num_isoforms = 10

    ### generate merge list
    merge_list = []
    if CFG['do_prune']:
        prune_tag = '_pruned'
    else:
        prune_tag = ''

    ### add all single bam file isoforms
    for i in range(len(CFG['samples'])):
        merge_list.append('%s/spladder/genes_graph_conf%i.%s%s.mat' % (CFG['out_dirname'], CFG['confidence_level'], CFG['samples'][i], prune_tag)

    ### add also isoforms of all bam files combined
    fn_bam_merged = sprintf('%s/spladder/genes_graph_conf%i.merge_bams%s.mat', CFG.out_dirname, CFG.confidence_level, prune_tag);
    if CFG['merge_strategy'] == 'merge_all')  and os.path.exists(fn_bam_merged):
        merge_list.append(fn_bam_merged)

    for i in range(len(merge_list)):
        ### load gene structure from sample i
        print 'Loading %s ...' % merge_list[i]
        genes = cPickle.load(open(merge_list[i]), 'r')
        print '... done'

        ### sort
        name_list = sp.array([x.name for x in genes])
        s_idx = sp.argsort(name_list)
        name_list = name_list[s_idx]
        genes = genes[s_idx]

        ### jump over first sample - nothig to add yet 
        if i == 0:
            genes2 = genes.copy()
            del genes
            continue

        ### did we append genes in the last round? --> re-sort
        if appended:
            name_list = {genes2(:).name};
            [name_list, s_idx] = sort(name_list);
            genes2 = genes2(s_idx);
            appended = False

        ### iterate over current genes
        g_idx = 0
        print 'Processing ...'
        for j in range(genes.shape[0]):
            if j % 100 == 0:
                print '.',
                if j % 1000 == 0:
                    print '%i/%i' % (j, genes.shape[0])
            g_idx_ = g_idx

            while genes2[g_idx].name < genes[j].name:
                g_idx += 1
            # same gene
            if genes2[g_idx].name == genes[j].name:
                ### pairwise comparison of all exons
                for k in range(len(genes[j].exons)):
                    broken = False
                    for l in range(len(genes2[g_idx].exons)):
                        ### transcripts are not identical in size --> continue
                        if not sp.all(genes[j].exons[k].shape == genes2[g_idx].exons[l].shape):
                            continue 
                        ### found transcript in genes2 that is identical to current transcript --> break
                        if sp.all(genes[j].exons[k] == genes2[g_idx].exons[l]):
                            broken = True
                            break 
                    ### we did not find any identical transcript in genes2 --> append transcript
                    if not broken:
                        genes2[g_idx].exons.append(genes[j].exons[k])
                        genes2[g_idx].transcripts.append('%s_%i' % (genes[j].transcripts[k], i))
                        genes2[g_idx].start = min(genes2[g_idx].start, genes[j].exons[k][:, 0].min())
                        genes2[g_idx].stop = max(genes2[g_idx].stop, genes[j].exons[k][:, 1].max())
                        genes2[g_idx].splicegraph = Splicegraph(genes2[g_idx])
            ### we did non find the gene name --> append new gene to genes2
            elif genes2[g_idx].name > genes[j].name:
                g_idx = g_idx_
                genes2 = sp.c_[genes2, genes[j]]
                appended = True
        print '... done\n'
        del genes

    genes = genes2
    del genes2

    fn = '%s/spladder/genes_graph_conf%i.%s%s_merge_isoforms.mat' % (CFG['out_dirname'], CFG['confidence_level'], CFG['merge_strategy'], prune_tag)
    print 'Store genes at: %s' % fn
    cPickle.dump(genes, open(fn, 'w'), -1)

    ### subsample transcripts if neccessary 
    print 'Subsample genes ...'
    for i in range(genes.shape[0]):
        if len(genes[i].transcripts) > max_num_isoforms:
            print 'Subsample for gene %i from %i transcripts.' % (i, len(genes[i].transcripts))
            r_idx = randperm(length(genes(i).transcripts));
            for r_idx in sorted(ransom.sample(range(len(genes[i].transcripts)), len(genes[i].transcripts) - max_num_isoforms), reverse=True):
                del genes[i].transcripts[r_idx]
                del genes[i].exons[r_idx]
            genes[i].start = min(genes[i].start, genes[i].exons[0][:, 0].min())
            genes[i].stop = max(genes[i].stop, genes[i].exons[-1][:, 1].max())
            genes.splicegraph = Splicegraph(genes[i])
    print '... done\n'

    fn = '%s/spladder/genes_graph_conf%i.%s%s_merge_isoforms_subsampled.mat' % (CFG['out_dirname'], CFG['confidence_level'], CFG['merge_strategy'], prune_tag)
    print 'Store subsampled genes at: %s' % fn
    cPickle.dump(genes, open(fn, 'w'), -1)
