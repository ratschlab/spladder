from scipy.sparse import lil_matrix
from numpy.matlib import repmat
import scipy as sp
import sys
import operator

import multiprocessing as mp 
import signal as sig

from ..helpers import log_progress
from functools import reduce

def detect_multipleskips(genes, gidx, log=False, edge_limit=300):
    # [idx_multiple_skips, exon_multiple_skips] = detect_multipleskips(genes, idx_alt) ;

    idx_multiple_skips = []
    exon_multiple_skips = []
    for iix, ix in enumerate(gidx):
        if log:
            sys.stdout.write('.')
            if (iix + 1) % 50 == 0:
                sys.stdout.write(' - %i/%i, found %i\n' % (iix + 1, genes.shape[0] + 1, len(idx_multiple_skips)))
            sys.stdout.flush()

        genes[iix].from_sparse()
        num_exons = genes[iix].splicegraph.get_len()
        edges = genes[iix].splicegraph.edges.copy()
        labels = repmat(sp.arange(num_exons), num_exons, 1).T
        genes[iix].to_sparse()
        
        # adjecency matrix: upper half only
        A = sp.zeros((num_exons, num_exons))
        for i in range(num_exons - 1):
            for j in range(i + 1, num_exons):
                A[i, j] = edges[i, j]
        
        # possible starting and ending exons of a multiple exon skip
        Pairs = lil_matrix((num_exons, num_exons))
        Ai = sp.dot(sp.dot(A, A), A) #paths of length 3
        while sp.any(Ai.ravel() > 0):
            coords = sp.where((A > 0 ) & (Ai > 0)) # multiple skip
            Pairs[coords[0], coords[1]] = 1
            Ai = sp.dot(Ai, A)  # paths of length ..+1
        
        edge = sp.where(Pairs.toarray() == 1)
        

        if edge[0].shape[0] > edge_limit:
            print('\nWARNING: not processing gene %i (%s); has %i edges; current limit is %i; adjust edge_limit to include.' % (ix, genes[iix].name, edge[0].shape[0], edge_limit))
            continue
        
        for cnt in range(edge[0].shape[0]):
            exon_idx_first = edge[0][cnt]
            exon_idx_last = edge[1][cnt]
      
            if edges[exon_idx_first, exon_idx_last] == 1:
                
                # find all pairs shortest path
                exist_path = sp.triu(edges).astype('double')
                exist_path[exon_idx_first, exon_idx_last] = 0
                exist_path[exist_path == 0] = sp.inf
                # set diagonal to 0
                exist_path[sp.arange(exist_path.shape[0]), sp.arange(exist_path.shape[0])] = 0
                
                long_exist_path = -sp.triu(edges).astype('double')
                long_exist_path[exon_idx_first, exon_idx_last] = 0
                long_exist_path[long_exist_path == 0] = sp.inf
                # set diagonal to 0
                long_exist_path[sp.arange(long_exist_path.shape[0]), sp.arange(long_exist_path.shape[0])] = 0
                
                path = sp.isfinite(exist_path) * labels
                long_path = sp.isfinite(long_exist_path) * labels
                
                for k in range(num_exons):
                    for i in range(num_exons):
                        idx = sp.where((exist_path[i, k] + exist_path[k, :]) < exist_path[i, :])[0]
                        exist_path[i, idx] = exist_path[i, k] + exist_path[k, idx]
                        path[i, idx] = path[k, idx]
                        
                        idx = sp.where((long_exist_path[i,k] + long_exist_path[k, :]) < long_exist_path[i, :])[0]
                        long_exist_path[i, idx] = long_exist_path[i, k] + long_exist_path[k, idx]
                        long_path[i, idx] = long_path[k, idx]
                
                temp_ix = sp.isfinite(long_exist_path)
                long_exist_path[temp_ix] = -long_exist_path[temp_ix]
                
                if (exist_path[exon_idx_first, exon_idx_last] > 2) and sp.isfinite(exist_path[exon_idx_first, exon_idx_last]):
                    backtrace = sp.array([path[exon_idx_first, exon_idx_last]])
                    while backtrace[-1] > exon_idx_first:
                        backtrace = sp.r_[backtrace, path[exon_idx_first, backtrace[-1]]]
                    backtrace = backtrace[:-1]
                    backtrace = backtrace[::-1]
                    idx_multiple_skips.append(ix) #repmat(ix, 1, backtrace.shape[0] + 2))
                    exon_multiple_skips.append([exon_idx_first, backtrace, exon_idx_last])
                elif (long_exist_path[exon_idx_first, exon_idx_last] > 2) and sp.isfinite(long_exist_path[exon_idx_first, exon_idx_last]):
                    backtrace = sp.array([long_path[exon_idx_first, exon_idx_last]])
                    while backtrace[-1] > exon_idx_first:
                        backtrace = sp.r_[backtrace, long_path[exon_idx_first, backtrace[-1]]]
                    backtrace = backtrace[:-1]
                    backtrace = backtrace[::-1]
                    idx_multiple_skips.append(ix) #repmat(ix, 1, backtrace.shape[0] + 2))
                    exon_multiple_skips.append([exon_idx_first, backtrace, exon_idx_last])

    if log:
        print('Number of multiple exon skips:\t\t\t\t\t%d' % len(idx_multiple_skips))

    return (idx_multiple_skips, exon_multiple_skips)


def detect_intronreten(genes, gidx, log=False, edge_limit=1000):
    # [idx_intron_reten,intron_intron_reten] = detect_intronreten(genes) ;

    idx_intron_reten = []
    intron_intron_reten = []
    for iix, ix in enumerate(gidx):
        if log:
            sys.stdout.write('.')
            if (iix + 1) % 50 == 0:
                sys.stdout.write(' - %i/%i, found %i\n' % (iix + 1, genes.shape[0] + 1, len(idx_intron_reten)))
            sys.stdout.flush()

        genes[iix].from_sparse()
        num_exons = genes[iix].splicegraph.get_len()
        vertices = genes[iix].splicegraph.vertices
        edges = genes[iix].splicegraph.edges.copy()
        genes[iix].to_sparse()

        if edges.shape[0] > edge_limit:
            print('\nWARNING: not processing gene %i (%s); has %i edges; current limit is %i; adjust edge_limit to include.' % (ix, genes[iix].name, edges.shape[0], edge_limit))
            continue
        
        #introns  = []
        introns = sp.zeros((0, 2), dtype='int')
        for exon_idx in range(num_exons - 1):  # start of intron
            idx = sp.where(edges[exon_idx, exon_idx + 1 : num_exons] == 1)[0]
            if idx.shape[0] == 0:
                continue
            idx += (exon_idx + 1)
            for exon_idx2 in idx: # end of intron
                #is_intron_reten = False
                if sp.sum((introns[:, 0] == vertices[1, exon_idx]) & (introns[:, 1] == vertices[0, exon_idx2])) > 0:
                    continue

                ### find shortest fully overlapping exon
                iidx = sp.where((vertices[0, :] < vertices[1, exon_idx]) & (vertices[1, :] > vertices[0, exon_idx2]))[0]
                if len(iidx) > 0:
                    iidx = iidx[sp.argmin(vertices[1, iidx] - vertices[0, iidx])]
                    idx_intron_reten.append(ix)
                    intron_intron_reten.append([exon_idx, exon_idx2, iidx])
                    introns = sp.r_[introns, [[vertices[1, exon_idx], vertices[0, exon_idx2]]]]

                #for exon_idx1 in range(num_exons): # exon
                #    # check that the exon covers the intron
                #    if (vertices[1, exon_idx] > vertices[0, exon_idx1]) and (vertices[0, exon_idx2] < vertices[1, exon_idx1]):
                #        is_intron_reten = True
                #        long_exon = exon_idx1 
                #        for l in range(len(introns)):
                #            if (vertices[1, exon_idx] == introns[l][0]) and (vertices[0, exon_idx2] == introns[l][1]):
                #                is_intron_reten = False
                #if is_intron_reten:
                #    idx_intron_reten.append(ix)
                #    intron_intron_reten.append([exon_idx, exon_idx2, long_exon])
                #    introns.append([vertices[1, exon_idx], vertices[0, exon_idx2]])

    if log:
        print('\nNumber of intron retentions:\t\t\t\t\t%d' % len(idx_intron_reten))

    return (idx_intron_reten, intron_intron_reten)


def detect_exonskips(genes, gidx, log=False, edge_limit=1000):
    # [idx_exon_skips, exon_exon_skips] = detect_exonskips(genes) ;

    idx_exon_skips = []
    exon_exon_skips = []
    for iix, ix in enumerate(gidx):
        if log:
            sys.stdout.write('.')
            if (iix + 1) % 50 == 0:
                sys.stdout.write(' - %i/%i, found %i\n' % (iix + 1, genes.shape[0] + 1, len(idx_exon_skips)))
            sys.stdout.flush()

        genes[iix].from_sparse()
        num_exons = genes[iix].splicegraph.get_len()
        edges = genes[iix].splicegraph.edges.copy()
        genes[iix].to_sparse()

        if edges.shape[0] > edge_limit:
            print('\nWARNING: not processing gene %i (%s); has %i edges; current limit is %i; adjust edge_limit to include.' % (ix, genes[iix].name, edges.shape[0], edge_limit))
            continue
        
        for exon_idx in range(num_exons - 2): #first exon
            for exon_idx1 in range(exon_idx + 1, num_exons - 1): # middle exon
                for exon_idx2 in range(exon_idx1 + 1, num_exons): # last exon
                    if (edges[exon_idx, exon_idx1] == 1) and edges[exon_idx, exon_idx2] and edges[exon_idx1, exon_idx2]:
                        idx_exon_skips.append(ix)
                        exon_exon_skips.append([exon_idx, exon_idx1, exon_idx2])
    if log:
        print('\nNumber of single exon skips:\t\t\t\t\t%d' % len(idx_exon_skips))

    return (idx_exon_skips, exon_exon_skips)



def detect_altprime(genes, gidx, log=False, edge_limit=1000):
    # function [idx_alt_5prime,exon_alt_5prime, idx_alt_3prime,exon_alt_3prime] ...
    #    = detect_altprime(genes);
    #
    # detect the alternative 5 and 3 prime ends of the intron. Note that 5 prime refers to the left
    # and 3 prime to the right for a positive strand

    MIN_OVERLAP = 11 # two exons that are alternative should overlap by at least the length

    idx_alt_5prime = []
    idx_alt_3prime = []
    exon_alt_5prime = []
    exon_alt_3prime = []

    for iix, ix in enumerate(gidx):
        if log:
            sys.stdout.write('.')
            if (iix + 1) % 50 == 0:
                sys.stdout.write(' - %i/%i, found %i + %i\n' % (iix + 1, genes.shape[0] + 1, len(idx_alt_3prime), len(idx_alt_5prime)))
            sys.stdout.flush()

        genes[iix].from_sparse()
        num_exons = genes[iix].splicegraph.get_len()
        vertices = genes[iix].splicegraph.vertices
        edges = genes[iix].splicegraph.edges.copy()
        strand = genes[iix].strand
        genes[iix].to_sparse()

        if edges.shape[0] > edge_limit:
            print('\nWARNING: not processing gene %i (%s); has %i edges; current limit is %i; adjust edge_limit to include.' % (ix, genes[iix].name, edges.shape[0], edge_limit))
            continue
        
        # Find alternative sites on the right of the intron,
        # same site on the left of the intron.
        for exon_idx in range(num_exons - 2):
            rightsites = []
            rightidx = []
            nr_exons = sp.sum(edges[exon_idx, exon_idx + 1:])
            if nr_exons >= 2:
                which_exons = sp.where(edges[exon_idx, exon_idx + 1:])[0] + exon_idx + 1
                exons = vertices[:, which_exons]
                for i in range(nr_exons - 1):
                    for j in range(i + 1, nr_exons):
                        # check that the left splice site of the exons are different
                        # make sure exons overlap - either:
                        # - left splice site of exon(i) is in exon(j)
                        # - left splice site of exon(j) is in exon(i)
                        # note that the 'overlap' relationship is not transitive
                        if ((exons[0, i] != exons[0, j]) and 
                          (((exons[0, i] > exons[0, j]) and (exons[0, i] < exons[1, j])) or ((exons[0, j] > exons[0, i]) and (exons[0, j] < exons[1, i]))) and
                          (min(exons[1, i], exons[1, j]) - max(exons[0, i], exons[0, j]) >= MIN_OVERLAP)):

                            assert(not ((exons[0, i] == exons[0, j]) and (exons[1, i] == exons[1, j])));
                            assert(exons[0, i] != exons[0, j])

                            ### add new events to the list
                            if not exons[0, i] in rightsites:
                                rightsites.append(exons[0, i])
                                rightidx.append(which_exons[i])
                            if not exons[0, j] in rightsites:
                                rightsites.append(exons[0, j])
                                rightidx.append(which_exons[j])
           
            # construct the output
            if len(rightsites) >= 2:
                if strand == '+':
                    exon_alt_3prime.append({'fiveprimesite':exon_idx, 'threeprimesites':rightidx})
                    idx_alt_3prime.append(ix)
                if strand == '-':
                    exon_alt_5prime.append({'threeprimesite':exon_idx, 'fiveprimesites':rightidx})
                    idx_alt_5prime.append(ix)
        
        # Find alternative sites on the left of the intron,
        # same site on the right of the intron.
        for exon_idx in range(2, num_exons):
            nr_exons = sp.sum(edges[:exon_idx, exon_idx])
            leftsites = []
            leftidx = []
            if nr_exons >= 2:
                which_exons = sp.where(edges[:exon_idx, exon_idx])[0]
                exons = vertices[:, which_exons]
                for i in range(nr_exons - 1):
                    for j in range(i + 1, nr_exons):
                        # check that the 5prime sites are different
                        # make sure exons overlap - either:
                        # - right splice site of exon(i) is in exon(j)
                        # - right splice site of exon(j) is in exon(i)
                        # note that the 'overlap' relationship is not transitive
                        if ((exons[1, i] != exons[1, j]) and
                          (((exons[1, i] <= exons[1, j]) and (exons[1, i] >= exons[0,j])) or ((exons[1, j] <= exons[1, i]) and (exons[1, j] >= exons[0, i]))) and
                          (min(exons[1, i], exons[1, j]) - max(exons[0, i], exons[0, j]) >= MIN_OVERLAP)):

                            assert(not((exons[0, i] == exons[0, j]) and (exons[1, i] == exons[1, j])))
                            assert(exons[1, i] != exons[1, j])
                        
                            # add new events to the list
                            if not exons[1, i] in leftsites:
                                leftsites.append(exons[1, i])
                                leftidx.append(which_exons[i])
                            if not exons[1, j] in leftsites:
                                leftsites.append(exons[1, j])
                                leftidx.append(which_exons[j])

            # construct the output
            if len(leftsites) >= 2:
                if strand == '+':
                    exon_alt_5prime.append({'threeprimesite':exon_idx, 'fiveprimesites':leftidx})
                    idx_alt_5prime.append(ix) 
                if strand == '-':
                    exon_alt_3prime.append({'fiveprimesite':exon_idx, 'threeprimesites':leftidx})
                    idx_alt_3prime.append(ix)

    if log:
        print('\nNumber of alternative 5 prime sites:\t\t\t\t%d' % len(idx_alt_5prime))
        print('Number of alternative 3 prime sites:\t\t\t\t%d' % len(idx_alt_3prime))

    return (idx_alt_5prime, exon_alt_5prime, idx_alt_3prime, exon_alt_3prime)


def detect_xorexons(genes, gidx, log=False, edge_limit=1000):
    #[idx_xor_exons, exon_xor_exons] = detect_xorexons(genes);

    idx_xor_exons = []
    exon_xor_exons = [] ### 5primesite of first exon, the 2 skipped
                        ### exons, 3primesite of last exon %%%
    for iix, ix in enumerate(gidx):

        if log:
            sys.stdout.write('.')
            if (iix + 1) % 50 == 0:
                sys.stdout.write(' - %i/%i, found %i\n' % (iix + 1, gidx.shape[0], len(idx_xor_exons)))
            sys.stdout.flush()

        genes[iix].from_sparse()
        num_exons = genes[iix].splicegraph.get_len()
        edges = genes[iix].splicegraph.edges.copy()
        vertices = genes[iix].splicegraph.vertices
        genes[iix].to_sparse()

        if edges.shape[0] > edge_limit:
            print('\nWARNING: not processing gene %i (%s); has %i edges; current limit is %i; adjust edge_limit to include.' % (ix, genes[iix].name, edges.shape[0], edge_limit))
            continue
        
        for exon_idx1 in range(num_exons - 3):
            for exon_idx2 in range(exon_idx1 + 1, num_exons - 2):
                if edges[exon_idx1, exon_idx2] == 1:
                    for exon_idx3 in range(exon_idx2 + 1, num_exons - 1):
                        if (edges[exon_idx1, exon_idx3] == 1) and (edges[exon_idx2, exon_idx3] == 0) and (vertices[0, exon_idx3] >= vertices[1, exon_idx2]):
                            for exon_idx4 in range(exon_idx3 + 1, num_exons):
                                if (edges[exon_idx2, exon_idx4] == 1) and (edges[exon_idx3, exon_idx4] == 1):
                                    idx_xor_exons.append(ix)
                                    exon_xor_exons.append([exon_idx1, exon_idx2, exon_idx3, exon_idx4])

    if log:
        print('\n\nNumber of XOR exons:\t\t\t\t\t\t%i\n' % len(idx_xor_exons))

    return (idx_xor_exons, exon_xor_exons)


def detect_wrapper(genes, event_type, gidx, idx, edge_limit, log=False):

    if event_type == 'mutex_exons':
        return (detect_xorexons(genes, gidx, log, edge_limit), idx)
    elif event_type == 'exon_skip':
        return (detect_exonskips(genes, gidx, log, edge_limit), idx)
    elif event_type == 'alt_prime':
        return (detect_altprime(genes, gidx, log, edge_limit), idx) 
    elif event_type == 'intron_retention':
        return (detect_intronreten(genes, gidx, log, edge_limit), idx)
    elif event_type == 'mult_exon_skip':
        return (detect_multipleskips(genes, gidx, log, edge_limit), idx)
    

def detect_events(genes, event_type, idx, options):

    if options.parallel > 1:
        pool = mp.Pool(processes=options.parallel, initializer=lambda: sig.signal(sig.SIGINT, sig.SIG_IGN))
        binsize = 3
        maxsize = idx.shape[0]
        idx_chunks = [sp.arange(x, min(x + binsize, maxsize)) for x in range(0, maxsize, binsize)]
        result_list = sp.empty((len(idx_chunks), ), dtype='object')

        try:
            result = [pool.apply_async(detect_wrapper, args=(genes[idx[cidx]], event_type, idx[cidx], c, options.detect_edge_limit)) for c,cidx in enumerate(idx_chunks)]
            res_cnt = 0
            while result:
                tmp = result.pop(0).get()
                result_list[tmp[1]] = tmp[0]
                if options.verbose:
                    log_progress(res_cnt, len(idx_chunks))
                    res_cnt += 1
            if options.verbose:
                log_progress(len(idx_chunks), len(idx_chunks))
                print('')
            pool.terminate()
            pool.join()
        except KeyboardInterrupt:
            print('Keyboard Interrupt - exiting', file=sys.stderr)
            pool.terminate()
            pool.join()
            sys.exit(1)

        ### integrate results in list into a coherent results list
        if len(result_list) > 0:
            result_list = [reduce(operator.add, [x[i] for x in result_list]) for i in range(len(result_list[0]))]
        else:
            if event_type == 'alt_prime':
                result_list = [[], [], [], []]
            else:
                result_list = [[], []]
    else:
        result_list = detect_wrapper(genes[idx], event_type, idx, None, edge_limit=options.detect_edge_limit, log=options.verbose)[0]        

    return result_list
 

