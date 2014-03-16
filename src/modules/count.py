import cPickle

def count_graph_coverage(genes, fn_bam=None, CFG=None, fn_out=None):
# [counts] = count_graph_coverage(genes, fn_bam, CFG, fn_out)

    if fn_bam is None and instanceof(genes, dict):
        PAR = genes
        genes = PAR['genes']
        fn_bam = PAR['fn_bam']
        if 'fn_out' in PAR:
            fn_out = PAR['fn_out'] 
        CFG = PAR['CFG']

    counts = sp.zeros((fn_bam.shape[0], genes.shape[0]), dtype='object')
    intron_tol = 0 

    for f in range(fn_bam.shape[0]):
        ### iterate over all genes and generate counts for
        ### the segments in the segment graph
        ### and the splice junctions in the splice graph
        for i in range(genes.shape[0]):
            print '.',
            gg = genes[i].copy()
            gg.start = gg.segmentgraph.segments.ravel().min()
            gg.stop = gg.segmentgraph.segments.ravel().max()

            ### add RNA-seq evidence to the gene structure
            (tracks, intron_list) = add_reads_from_bam(gg, fn_bam[f], ['exon_track','intron_list'], CFG['read_filter'], CFG['var_aware']);

            ### extract mean exon coverage for all segments
            counts[f, i] = Counts(gg.segmentgraph.segments.shape[1])

            for j in range(gg.segmentgraph.segments.shape[1]):
                idx = sp.arange(gg.segmentgraph.segments[0, j], gg.segmentgraph.segments[1, j) - gg.start
                counts[f, i].segments[j] = sp.mean(sp.sum(tracks[:, idx], axis=0))
                counts[f, i].seg_pos[j] = sp.sum(sp.sum(tracks[:, idx], axis=0) > 0)

            ### extract intron counts 
            k, l = sp.where(gg.segmentgraph.seg_edges == 1)
            for m in range(k.shape[0]):
                idx = sp.where((sp.absolute(intron_list[:, 0] - gg.segmentgraph.segments[1, k[m]]) <= intron_tol) & (sp.absolute(intron_list[:, 1] - gg.segmentgraph.segments[0, l[m]]) <= intron_tol))[0]
                if counts[f, i].edges.shape[0] == 0:
                    if idx.shape[0] > 0:
                        counts[f, i].edges = sp.array([sp.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), sp.sum(intron_list[idx, 2])])
                    else:
                        counts[f, i].edges = sp.array([sp.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), 0])
                else:
                    if idx.shape[0] > 0:
                        counts[f, i].edges = sp.r_[counts[f, i], [sp.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), sp.sum(intron_list[idx, 2])]]
                    else:
                        counts[f, i].edges = sp.r_[counts[f, i], [sp.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), 0]]

    if fn_out is not None:
        cPickle.dump(open(fn_out, 'w'), counts, -1)



def count_graph_coverage_wrapper(fname_in, fname_out, CFG):

    genes = cPickle.load(open(fname_in, 'r'))
    
    if genes.segmentgraph is None:
        genes = create_segment_graph(genes, CFG)
        cPickle.dump(open(fname_in, 'w'), genes, -1)

    if CFG['rproc']:
        for s_idx in range(CFG['strains'].shape[0]):
            print '%i/%i\r' % (j, CFG['strains'].shape[0])
            counts = count_graph_coverage(genes, CFG['list_bam'][s_idx], CFG)
    else:
        chunksize = 10
        jobinfo = rproc_empty(0)
        PAR = dict()
        PAR['CFG'] = CFG
        for c_idx in range(0, genes.shape[0], chunksize):
            cc_idx = min(genes.shape[0], c_idx + chunksize - 1)
            fn = fname_out.replace('.mat', '.chunk_%i_%i.mat' % (c_idx, cc_idx))
            if os.path.exists(fn):
                continue
            else:
                print 'submitting chunk %i to %i' % (c_idx, cc_idx)
                PAR['genes'] = genes[c_idx:cc_idx]
                PAR['fn_bam'] = CFG['bam_fnames']
                PAR['fn_out'] = fn
                PAR['CFG'] = CFG
                jobinfo.append(rproc('count_graph_coverage', PAR, 30000, CFG['options_rproc'], 48*60))

        jobinfo = rproc_wait(jobinfo, 30, 1, 1) ;

        ### merge results
        counts = sp.zeros((CFG['strains'].shape[0], genes.shape[0]), dtype='object')
        for c_idx in range(0, genes.shape[0], chunksize):
            cc_idx = min(genes.shape[0], c_idx + chunksize - 1)
            fn = fname_out.replace('.mat', '.chunk_%i_%i.mat' % (c_idx, cc_idx))
            if not os.path.exists(fn):
                print >> sys.stderr, 'ERROR: Not all chunks in counting graph coverage completed!'
                sys.exit(1)
            else:
                counts_ = cPickle.load(open(fn, 'r'))
                counts[:, c_idx:cc_idx] = counts_

    cPickle.dump(open(fname_out, 'w'), counts, -1)
