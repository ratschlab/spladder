from scipy.sparse import lil_matrix
import numpy as np
from numba import jit
from numba.core import types
from numba.typed import List
from numba.typed import Dict

import sys
import operator

import multiprocessing as mp 
import signal as sig

from ..helpers import log_progress
from functools import reduce

@jit(nopython=True)
def _detect_multipleskips_fastcore(ix, gene_name, num_exons, edges, labels, edge_limit, idx_multiple_skips, exon_multiple_skips):

    # adjecency matrix: upper half only
    A = np.zeros((num_exons, num_exons))
    for i in range(num_exons - 1):
        for j in range(i + 1, num_exons):
            A[i, j] = edges[i, j]
    
    # possible starting and ending exons of a multiple exon skip
    #Pairs = lil_matrix((num_exons, num_exons))
    edge = set()
    Ai = np.dot(np.dot(A, A), A) #paths of length 3
    while np.any(Ai.ravel() > 0):
        coords = np.where((A > 0 ) & (Ai > 0)) # multiple skip
        for i in range(len(coords[0])):
            #Pairs[coords[0], coords[1]] = 1
            edge.add((coords[0][i], coords[1][i]))
        Ai = np.dot(Ai, A)  # paths of length ..+1
    
    if len(edge) > edge_limit:
        print('\nWARNING: not processing gene ' + str(ix) + ' (' + gene_name + '); has ' + str(len(edge)) + ' edges; current limit is ' + str(edge_limit))
        return
    
    bt_start = len(exon_multiple_skips[3]) - 1 ### subtract 1 because the first element is just the type placeholder
    for _edge in sorted(edge):
        exon_idx_first = _edge[0] #[cnt]
        exon_idx_last = _edge[1] #Y[cnt]
  
        if edges[exon_idx_first, exon_idx_last] == 1:
            
            # find all pairs shortest path
            exist_path = np.triu(edges).astype(np.double)
            exist_path[exon_idx_first, exon_idx_last] = 0
            for i1 in range(exist_path.shape[0]):
                for i2 in range(exist_path.shape[1]):
                    if exist_path[i1, i2] == 0:
                        exist_path[i1, i2] = np.inf
            # set diagonal to 0
            for i1 in range(exist_path.shape[0]):
                exist_path[i1, i1] = 0
            
            long_exist_path = -np.triu(edges).astype(np.double)
            long_exist_path[exon_idx_first, exon_idx_last] = 0
            long_exist_path[exon_idx_first, exon_idx_last] = 0
            for i1 in range(long_exist_path.shape[0]):
                for i2 in range(long_exist_path.shape[1]):
                    if long_exist_path[i1, i2] == 0:
                        long_exist_path[i1, i2] = np.inf
            # set diagonal to 0
            for i1 in range(long_exist_path.shape[0]):
                long_exist_path[i1, i1] = 0
            
            path = np.isfinite(exist_path) * labels
            long_path = np.isfinite(long_exist_path) * labels
            
            for k in range(num_exons):
                for i in range(num_exons):
                    idx = np.where((exist_path[i, k] + exist_path[k, :]) < exist_path[i, :])[0]
                    for _idx in idx:
                        exist_path[i, _idx] = exist_path[i, k] + exist_path[k, _idx]
                        path[i, _idx] = path[k, _idx]
                    
                    idx = np.where((long_exist_path[i,k] + long_exist_path[k, :]) < long_exist_path[i, :])[0]
                    for _idx in idx:
                        long_exist_path[i, _idx] = long_exist_path[i, k] + long_exist_path[k, _idx]
                        long_path[i, _idx] = long_path[k, _idx]
            
            _idx1, _idx2 = np.where(np.isfinite(long_exist_path))
            for i1, i2 in zip(_idx1, _idx2):
                long_exist_path[i1, i2] = -long_exist_path[i1, i2]
            
            if (exist_path[exon_idx_first, exon_idx_last] > 2) and np.isfinite(exist_path[exon_idx_first, exon_idx_last]):
                backtrace = [path[exon_idx_first, exon_idx_last]]
                while backtrace[-1] > exon_idx_first:
                    backtrace.append(path[exon_idx_first, backtrace[-1]])
                backtrace = backtrace[:-1]
                backtrace = backtrace[::-1]
                idx_multiple_skips.append(ix) 
                #exon_multiple_skips.append([exon_idx_first, backtrace, exon_idx_last])
                exon_multiple_skips[0].append(exon_idx_first)
                exon_multiple_skips[1].append(bt_start)
                exon_multiple_skips[2].append(bt_start + len(backtrace))
                bt_start = bt_start + len(backtrace)
                exon_multiple_skips[3].extend(backtrace)
                exon_multiple_skips[4].append(exon_idx_last)
            elif (long_exist_path[exon_idx_first, exon_idx_last] > 2) and np.isfinite(long_exist_path[exon_idx_first, exon_idx_last]):
                backtrace = [long_path[exon_idx_first, exon_idx_last]]
                while backtrace[-1] > exon_idx_first:
                    backtrace.append(long_path[exon_idx_first, backtrace[-1]])
                backtrace = backtrace[:-1]
                backtrace = backtrace[::-1]
                idx_multiple_skips.append(ix) 
                #exon_multiple_skips.append([exon_idx_first, backtrace, exon_idx_last])
                exon_multiple_skips[0].append(exon_idx_first)
                exon_multiple_skips[1].append(bt_start)
                exon_multiple_skips[2].append(bt_start + len(backtrace))
                bt_start = bt_start + len(backtrace)
                exon_multiple_skips[3].extend(backtrace)
                exon_multiple_skips[4].append(exon_idx_last)

def detect_multipleskips(genes, gidx, log=False, edge_limit=300):
    # [idx_multiple_skips, exon_multiple_skips] = detect_multipleskips(genes, idx_alt) ;

    idx_multiple_skips = List()
    idx_multiple_skips.append(0)
    exon_multiple_skips = List()
    for i in range(5): # first bt_start bt_end bt last
        tmp = List()
        tmp.append(0)
        exon_multiple_skips.append(tmp)

    for iix, ix in enumerate(gidx):
        if log:
            sys.stdout.write('.')
            if (iix + 1) % 50 == 0:
                sys.stdout.write(' - %i/%i, found %i\n' % (iix + 1, genes.shape[0] + 1, len(idx_multiple_skips) - 1))
            sys.stdout.flush()

        genes[iix].from_sparse()
        num_exons = genes[iix].splicegraph.get_len()
        edges = genes[iix].splicegraph.edges.copy()
        labels = np.hstack([np.arange(num_exons)[:, np.newaxis]] * num_exons)
        genes[iix].to_sparse()
     
        _detect_multipleskips_fastcore(ix, genes[iix].name, num_exons, edges, labels, edge_limit, idx_multiple_skips, exon_multiple_skips)

    #### assemble outputs
    _exon_multiple_skips = []
    for i in range(5):
        _exon_multiple_skips.append(np.array(exon_multiple_skips[i][1:]))
    exon_multiple_skips = _exon_multiple_skips

    _exon_multiple_skips = []
    for i in range(len(exon_multiple_skips[0])):
        _bt = exon_multiple_skips[3][exon_multiple_skips[1][i]:exon_multiple_skips[2][i]]
        _exon_multiple_skips.append([exon_multiple_skips[0][i], _bt, exon_multiple_skips[4][i]])
    exon_multiple_skips = _exon_multiple_skips

    if exon_multiple_skips == []:
        idx_multiple_skips = []
    else:
        idx_multiple_skips = [_ for _ in idx_multiple_skips[1:]]
        
    if log:
        print('Number of multiple exon skips:\t\t\t\t\t%d' % len(idx_multiple_skips))

    return (idx_multiple_skips, exon_multiple_skips)

@jit(nopython=True)
def _detect_intronreten_fastcore(ix, num_exons, vertices, edges, idx_intron_reten, intron_intron_reten):

    introns = np.zeros((0, 2), dtype=np.int64)
    for exon_idx in range(num_exons - 1):  # start of intron
        idx = np.where(edges[exon_idx, exon_idx + 1 : num_exons] == 1)[0]
        if idx.shape[0] == 0:
            continue
        idx += (exon_idx + 1)
        for exon_idx2 in idx: # end of intron
            #is_intron_reten = False
            if np.sum((introns[:, 0] == vertices[1, exon_idx]) & (introns[:, 1] == vertices[0, exon_idx2])) > 0:
                continue

            ### find shortest fully overlapping exon
            iidx = np.where((vertices[0, :] < vertices[1, exon_idx]) & (vertices[1, :] > vertices[0, exon_idx2]))[0]
            if len(iidx) > 0:
                iidx = iidx[np.argmin(vertices[1, :][iidx] - vertices[0, :][iidx])]
                idx_intron_reten.append(ix)
                #intron_intron_reten.append([exon_idx, exon_idx2, iidx])
                intron_intron_reten[0].append(exon_idx)
                intron_intron_reten[1].append(exon_idx2)
                intron_intron_reten[2].append(iidx)
                tmp = np.array([[vertices[1, exon_idx], vertices[0, exon_idx2]]])
                introns = np.vstack((introns, tmp)) #[[vertices[1, exon_idx], vertices[0, exon_idx2]]]])

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


def detect_intronreten(genes, gidx, log=False, edge_limit=1000):
    # [idx_intron_reten,intron_intron_reten] = detect_intronreten(genes) ;

    idx_intron_reten = List()
    idx_intron_reten.append(0)
    intron_intron_reten = List()
    for i in range(3):
        tmp = List()
        tmp.append(0)
        intron_intron_reten.append(tmp)

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
        
        _detect_intronreten_fastcore(ix, num_exons, vertices, edges, idx_intron_reten, intron_intron_reten)

    ### assemble outputs
    intron_intron_reten = [[intron_intron_reten[0][i], intron_intron_reten[1][i], intron_intron_reten[2][i]] for i in range(1, len(intron_intron_reten[0]))]
    if intron_intron_reten == []:
        idx_intron_reten = []
    else:
        idx_intron_reten = [_ for _ in idx_intron_reten[1:]]

    if log:
        print('\nNumber of intron retentions:\t\t\t\t\t%d' % len(idx_intron_reten))

    return (idx_intron_reten, intron_intron_reten)

@jit(nopython=True)
def _detect_exonskips_fast(ix, num_exons, edges, idx_exon_skips, exon_exon_skips):

    for exon_idx in range(num_exons - 2): #first exon
        for exon_idx1 in range(exon_idx + 1, num_exons - 1): # middle exon
            for exon_idx2 in range(exon_idx1 + 1, num_exons): # last exon
                if (edges[exon_idx, exon_idx1] == 1) and edges[exon_idx, exon_idx2] and edges[exon_idx1, exon_idx2]:
                    idx_exon_skips.append(ix)
                    exon_exon_skips[0].append(exon_idx)
                    exon_exon_skips[1].append(exon_idx1)
                    exon_exon_skips[2].append(exon_idx2)

def detect_exonskips(genes, gidx, log=False, edge_limit=1000):
    # [idx_exon_skips, exon_exon_skips] = detect_exonskips(genes) ;

    idx_exon_skips = List() 
    idx_exon_skips.append(0)
    exon_exon_skips = List()
    for i in range(3):
        tmp = List()
        tmp.append(0)
        exon_exon_skips.append(tmp)

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

        _detect_exonskips_fast(ix, num_exons, edges, idx_exon_skips, exon_exon_skips)

    ### assemble outputs
    exon_exon_skips = [[exon_exon_skips[0][i], exon_exon_skips[1][i], exon_exon_skips[2][i]] for i in range(1, len(exon_exon_skips[0]))]
    if exon_exon_skips == []:
        idx_exon_skips = []
    else:
        idx_exon_skips = [_ for _ in idx_exon_skips[1:]]

    if log:
        print('\nNumber of single exon skips:\t\t\t\t\t%d' % len(idx_exon_skips))

    return (idx_exon_skips, exon_exon_skips)

@jit(nopython=True)
def _detect_altprime_fastcore(ix, num_exons, vertices, edges, strand, idx_alt_5prime, idx_alt_3prime, exon_alt_5prime, exon_alt_3prime):

    MIN_OVERLAP = 11 # two exons that are alternative should overlap by at least the length

    # Find alternative sites on the right of the intron,
    # same site on the left of the intron.
    for exon_idx in range(num_exons - 2):
        rightsites = []
        rightidx = []
        nr_exons = np.sum(edges[exon_idx, exon_idx + 1:])
        if nr_exons >= 2:
            which_exons = np.where(edges[exon_idx, exon_idx + 1:])[0] + exon_idx + 1
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
                exon_alt_3prime.append({'fiveprimesite':np.array([exon_idx]), 'threeprimesites':np.array(rightidx)})
                idx_alt_3prime.append(ix)
            if strand == '-':
                exon_alt_5prime.append({'threeprimesite':np.array([exon_idx]), 'fiveprimesites':np.array(rightidx)})
                idx_alt_5prime.append(ix)
    
    # Find alternative sites on the left of the intron,
    # same site on the right of the intron.
    for exon_idx in range(2, num_exons):
        nr_exons = np.sum(edges[:exon_idx, exon_idx])
        leftsites = []
        leftidx = []
        if nr_exons >= 2:
            which_exons = np.where(edges[:exon_idx, exon_idx])[0]
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
                exon_alt_5prime.append({'threeprimesite':np.array([exon_idx]), 'fiveprimesites':np.array(leftidx)})
                idx_alt_5prime.append(ix) 
            if strand == '-':
                exon_alt_3prime.append({'fiveprimesite':np.array([exon_idx]), 'threeprimesites':np.array(leftidx)})
                idx_alt_3prime.append(ix)


def detect_altprime(genes, gidx, log=False, edge_limit=1000):
    # function [idx_alt_5prime,exon_alt_5prime, idx_alt_3prime,exon_alt_3prime] ...
    #    = detect_altprime(genes);
    #
    # detect the alternative 5 and 3 prime ends of the intron. Note that 5 prime refers to the left
    # and 3 prime to the right for a positive strand

    idx_alt_5prime = List()
    idx_alt_5prime.append(0)
    idx_alt_3prime = List()
    idx_alt_3prime.append(0)
    exon_alt_5prime = List()
    d = Dict.empty(
        key_type=types.unicode_type,
        value_type=types.int64[:],
    )
    exon_alt_5prime.append(d) #{'fiveprimesite':0, 'threeprimesites':np.arange(2)})
    exon_alt_3prime = List()
    exon_alt_3prime.append(d) #{'threeprimesite':0, 'fiveprimesites':np.arange(2)})

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

        _detect_altprime_fastcore(ix, num_exons, vertices, edges, strand, idx_alt_5prime, idx_alt_3prime, exon_alt_5prime, exon_alt_3prime)

    ### assemble outputs
    exon_alt_5prime = [{'threeprimesite':_['threeprimesite'][0], 'fiveprimesites':_['fiveprimesites']} for _ in exon_alt_5prime[1:]]
    exon_alt_3prime = [{'fiveprimesite':_['fiveprimesite'][0], 'threeprimesites':_['threeprimesites']} for _ in exon_alt_3prime[1:]]
    if exon_alt_5prime == []:
        idx_alt_5prime = []
    else:
        idx_alt_5prime = [_ for _ in idx_alt_5prime[1:]]
    if exon_alt_3prime == []:
        idx_alt_3prime = []
    else:
        idx_alt_3prime = [_ for _ in idx_alt_3prime[1:]]
       
    if log:
        print('\nNumber of alternative 5 prime sites:\t\t\t\t%d' % len(idx_alt_5prime))
        print('Number of alternative 3 prime sites:\t\t\t\t%d' % len(idx_alt_3prime))

    return (idx_alt_5prime, exon_alt_5prime, idx_alt_3prime, exon_alt_3prime)


@jit(nopython=True)
def _detect_xorexons_fastcore(ix, gene_name, num_exons, edges, vertices, edge_limit, idx_xor_exons, exon_xor_exons):
    for exon_idx1 in range(num_exons - 3):
        for exon_idx2 in range(exon_idx1 + 1, num_exons - 2):
            if edges[exon_idx1, exon_idx2] == 1:
                for exon_idx3 in range(exon_idx2 + 1, num_exons - 1):
                    if (edges[exon_idx1, exon_idx3] == 1) and (edges[exon_idx2, exon_idx3] == 0) and (vertices[0, exon_idx3] >= vertices[1, exon_idx2]):
                        for exon_idx4 in range(exon_idx3 + 1, num_exons):
                            if (edges[exon_idx2, exon_idx4] == 1) and (edges[exon_idx3, exon_idx4] == 1):
                                idx_xor_exons.append(ix)
                                #exon_xor_exons.append([exon_idx1, exon_idx2, exon_idx3, exon_idx4])
                                exon_xor_exons[0].append(exon_idx1)
                                exon_xor_exons[1].append(exon_idx2)
                                exon_xor_exons[2].append(exon_idx3)
                                exon_xor_exons[3].append(exon_idx4)


def detect_xorexons(genes, gidx, log=False, edge_limit=1000):
    #[idx_xor_exons, exon_xor_exons] = detect_xorexons(genes);

    idx_xor_exons = List() 
    idx_xor_exons.append(0)
    exon_xor_exons = List() ### 5primesite of first exon, the 2 skipped
                            ### exons, 3primesite of last exon %%%
    for i in range(4):
        tmp = List()
        tmp.append(0)
        exon_xor_exons.append(tmp)

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
        
        _detect_xorexons_fastcore(ix, genes[iix].name, num_exons, edges, vertices, edge_limit, idx_xor_exons, exon_xor_exons)

    ### assemble outputs
    exon_xor_exons = [[exon_xor_exons[0][i], exon_xor_exons[1][i], exon_xor_exons[2][i], exon_xor_exons[3][i]] for i in range(1, len(exon_xor_exons[0]))]
    if exon_xor_exons == []:
        idx_xor_exons = []
    else:
        idx_xor_exons = [_ for _ in idx_xor_exons[1:]]

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
        idx_chunks = [np.arange(x, min(x + binsize, maxsize)) for x in range(0, maxsize, binsize)]
        result_list = np.empty((len(idx_chunks), ), dtype='object')

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
 

