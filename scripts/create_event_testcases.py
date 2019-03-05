import sys
import os
import numpy.random as npr
import scipy as sp


def rev_complement(seq):
    rev_dict =  {"A": "T", "T": "A", "C": "G", "G": "C"}  
    return ''.join([rev_dict[_] for _ in seq][::-1])


def create_test_genome(L, seed=23):

    ### generate a random DNA string of length L
    npr.seed(seed)
    map_dict = {0:'A', 1:'G', 2:'C', 3:'T'}
    return ''.join([map_dict[i] for i in npr.randint(0, 4, L)])


def _adapt_splice_sites(exons, genome, strand='+'):
    """Makes sure we always use canonical splice consensus in the genome"""

    for e in range(1, len(exons) - 1):
        if strand == '-':
            ### acceptor
            genome = genome[:exons[e-1][1]] + 'CT' + genome[exons[e-1][1] + 2:]
            ### donor
            genome = genome[:exons[e][0]-2] + 'AC' + genome[exons[e][0]:]
        else:
            ### donor
            genome = genome[:exons[e-1][1]] + 'GT' + genome[exons[e-1][1] + 2:]
            ### acceptor
            genome = genome[:exons[e][0]-2] + 'AG' + genome[exons[e][0]:]
        

def _get_sequence(exons, genome, strand='+'):
    """Returns the DNA sequence of a given transcript"""

    seq = []
    pos = []
    for exon in exons:
        seq.append(genome[exon[0]:exon[1]])
        pos.extend(range(exon[0], exon[1]))
    if strand == '-':
        return (rev_complement(''.join(seq)), pos)
    else:
        return (''.join(seq), pos)


def _write_fasta(fname, seq_id, seq, width=80):
    
    with open(fname, 'w') as fh:
        fh.write('>%s\n' % seq_id)
        for i in range(0, len(seq), width):
            fh.write(seq[i:i+width] + '\n')


def _write_fastq(fh, read_id, seq):

    fh.write('@%s/1\n' % (read_id))
    fh.write(seq + '\n')
    fh.write('+\n')
    fh.write(len(seq)*'G' + '\n')

def parse_abundances(fname):
    
    abundances = dict()
    for line in open(fname, 'r'):
        sl = line.strip().split('\t')
        abundances[(sl[0], sl[1])] = [int(x) for x in sl[2:]]
    return abundances


def parse_blocks(fname):
    blocks = []
    for line in open(fname, 'r'):
        blocks.append(line.strip().split('\t'))
    return blocks
        

def create_reads(outfile, genome, genes, abundances, idx, testname, readlen=50, seed=23):

    npr.seed(seed)

    ### genes is dictionary with 
    # key1 --> (gene ID, strand)
    # key2 --> trans ID
    # value --> list of ([start, stop)] of exon)
    ### abundances is dictionary with
    # key --> (gene ID, trans ID)
    # value --> list of abundance as transcript count (choose correct entry with idx)

    cov = sp.zeros((len(genome), ), dtype='int')
    with open(outfile, 'w') as fh:
        for gene in genes:
            for trans in genes[gene]:

                ### generate transcript sequence
                _adapt_splice_sites(genes[gene][trans], genome, gene[1])
                seq, pos = _get_sequence(genes[gene][trans], genome, gene[1])

                ### get number of reads
                readnum = int(len(seq) * abundances[(gene[0], trans)][idx] / readlen)

                ### sample from sequence
                for i, start in enumerate(npr.randint(0, len(seq) - readlen, readnum)):
                    _write_fastq(fh, '{}:{}:read_{}'.format(testname, seed, i+1), seq[start:start + readlen])
                    cov[pos[start:start + readlen]] += 1
   
def gen_annotation(blocks, outfile):

    ### collect gene boundaries
    gene_bounds = dict()
    for block in blocks:
        exons = [x.split(':') for x in block[4].split(',')]
        if not block[0] in gene_bounds:
            gene_bounds[block[0]] = (int(exons[0][0]), int(exons[-1][1]))
        else:
            gene_bounds[block[0]] = (min(int(exons[0][0]), gene_bounds[block[0]][0]), max(int(exons[-1][1]), gene_bounds[block[0]][1]))

    genes = dict()
    with open(outfile, 'w') as fh:
        for block in blocks:
            gene_id = block[0]
            trans_id = block[1]
            chrm = block[2]
            strand = block[3]
            exons = [x.split(':') for x in block[4].split(',')]

            if not (gene_id, strand) in genes:
                genes[(gene_id, strand)] = dict()
            genes[(gene_id, strand)][trans_id] = [[int(x[0]) - 1, int(x[1])] for x in exons]

            fh.write('\t'.join([chrm, 'sim', 'gene', str(gene_bounds[block[0]][0]), str(gene_bounds[block[0]][1]), '.', strand, '.', 'gene_id "{}"; gene_type "protein_coding"; gene_name "{}";\n'.format(gene_id, gene_id)]))
            fh.write('\t'.join([chrm, 'sim', 'transcript', exons[0][0], exons[-1][1], '.', strand, '.', 'gene_id "{}"; transcript_id "{}.{}"; gene_type "protein_coding"; gene_name "{}";\n'.format(gene_id, gene_id, trans_id, gene_id)]))
            for e,exon in enumerate(exons):
                fh.write('\t'.join([chrm, 'sim', 'exon', exon[0], exon[1], '.', strand, '.', 'gene_id "{}"; transcript_id "{}.{}"; gene_type "protein_coding"; gene_name "{}";\n'.format(gene_id, gene_id, trans_id, gene_id)])) 

    return genes
            

def main(argv):
    if len(argv) < 4:
        sys.stderr.write('Usage: %s <out_dir> <test_name> <blocks> <abundances>>\n' % argv[0])
        sys.exit(1)
    outdir = argv[1]
    testname = argv[2]
    blocks = parse_blocks(argv[3])
    abundances = parse_abundances(argv[4])

    ### generate genome
    genome = create_test_genome(10000)
    _write_fasta(os.path.join(outdir, 'genome.fa'), '1', genome)

    ### create annotation from blocks 
    genes = gen_annotation(blocks, os.path.join(outdir, '{}.gtf'.format(testname)))

    ### create fastq files
    myseed=1
    k = next(iter(abundances.keys()))
    for i in range(len(abundances[k])):
        outfile = os.path.join(outdir, '{}_{}_sample{}.fq'.format(testname, myseed, i+1))
        create_reads(outfile, genome, genes, abundances, i, testname, seed=myseed)

if __name__ == "__main__":
    main(sys.argv)

