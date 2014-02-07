import pysam
import re
import scipy as sp
import scipy.sparse

def get_reads(fname, chr_name, start, stop, strand = None, filter = None, mapped = True, spliced = True, var_aware = None, collapse = False):
    
    infile = pysam.Samfile(fname, 'rb')

    ### vectors to build sparse matrix
    i = []
    j = []

    read_cnt = 0
    introns = []

    ### pysam query is zero based in position (results are as well), all intervals are pythonic half open
    for read in infile.fetch(chr_name, start, stop):
        if read.is_unmapped:
            continue
        if 3 in [x[0] for x in read.cigar]:
            if not spliced:
                continue
        else:
            if not mapped:
                continue

        if filter is not None:
            ### handle mismatches
            if var_aware:
                mm = sum([x[1] for x in read.tags if x[0] in ['XM', 'XG']])
            else:
                mm = sum([x[1] for x in read.tags if x[0] == 'NM'])
            if filter['mismatch'] < mm:
                continue
            ### handle min segment length
            ### remove all elements from CIGAR sting that do not 
            ### contribute to the segments (hard- and softclips and insertions)
            cig = re.sub(r'[0-9]*[HSI]', '', read.cigarstring)
            ### split the string at the introns and sum the remaining segment elements, compare to filter
            if min([sum([int(y) for y in re.split('[MD]', x)[:-1]]) for x in re.split('[0-9]*N', cig)]) < filter['min_exon_len']:
                continue
        
        ### get introns
        if 3 in [x[0] for x in read.cigar]:
            p = read.pos 
            for o in read.cigar:
                if o[0] == 3:
                    introns.append([p, p + o[1]])
                if o[0] in [0, 2, 3]:
                    p += o[1]

        ### get coverage
        for p in read.positions:
            if p - start >= 0:
                if p < start or p >= stop:
                    break
                else:
                    i.append(read_cnt)
                    j.append(p - start)

        read_cnt += 1

    i = sp.array(i)
    j = sp.array(j)

    ### construct sparse matrix
    read_matrix = scipy.sparse.coo_matrix((sp.ones(i.shape), (i, j)), shape = (read_cnt, stop - start), dtype = 'int')

    ### construct introns
    introns = sp.array(introns, dtype = 'int')

    if collapse:
        return (read_matrix.sum(axis = 0), introns)
    else:
        return (read_matrix, introns)
