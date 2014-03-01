import sys
import cPickle
import os
import pysam
import scipy as sp
from .classes.gene import Gene
from .classes.region import Region
from .classes.splicegraph import Splicegraph

def get_tags_gff3(tagline):
    """Extract tags from given tagline"""

    tags = dict()
    for t in tagline.split(';'):
        tt = t.split('=')
        tags[tt[0]] = tt[1]
    return tags


def get_tags_gtf(tagline):
    """Extract tags from given tagline"""

    tags = dict()
    for t in tagline.strip(';').split(';'):
        tt = t.strip(' ').split(' ')
        tags[tt[0]] = tt[1].strip('"')
    return tags


def init_genes_gtf(infile, CFG=None, outfile=None):
    # genes = init_genes_gtf(infile, CFG=None, outfile=None)

    """This function reads the gtf input file and returns the information in an
       internal data structure"""

    if CFG is not None and CFG['verbose']:
        print >> sys.stderr, "Parsing annotation from %s ..." % infile
    
    ### initial run to get the transcript to gene mapping
    if CFG is not None and CFG['verbose']:
        print >> sys.stderr, "... init structure"

    genes = dict()
    chrms = []

    for line in open(infile, 'r'):
        if line[0] == '#':
            continue
        sl = line.strip().split('\t')
        tags = get_tags_gtf(sl[8])
        if sl[2] in ['gene', 'Gene']:
            try:
                start = int(sl[3]) - 1
            except ValueError:
                start =  -1
            try:
                stop = int(sl[4])
            except ValueError:
                stop = -1
            genes[tags['gene_id']] = Gene(name=tags['gene_id'], start=start, stop=stop, chr=sl[0], strand=sl[6], source=sl[1], gene_type=tags['gene_type'])
            chrms.append(sl[0])

    chrms = sp.sort(sp.unique(chrms))
    for gid in genes:
        genes[gid].chr_num = sp.where(chrms == genes[gid].chr)[0]

    counter = 1
    for line in open(infile, 'r'):
        if CFG is not None and CFG['verbose'] and counter % 10000 == 0:
            print >> sys.stderr, '.',
        counter += 1        

        if line[0] == '#':
            continue
        sl = line.strip().split('\t')
        
        ### get start and end
        try:
            start = int(sl[3]) - 1
        except ValueError:
            start =  -1
        try:
            stop = int(sl[4])
        except ValueError:
            stop = -1

        ### get tags
        tags = get_tags_gtf(sl[8])

        ### add exons
        if sl[2] in ['exon', 'Exon']:
            trans_id = tags['transcript_id']
            gene_id = tags['gene_id']
            try:
                t_idx = genes[gene_id].transcripts.index(trans_id)
            except ValueError:
                t_idx = len(genes[gene_id].transcripts)
                genes[gene_id].transcripts.append(trans_id)
            genes[gene_id].add_exon(sp.array([int(sl[3]) - 1, int(sl[4]) - 1], dtype='int'), idx=t_idx)

    ### add splicegraphs
    for gene in genes:
        genes[gene].splicegraph = Splicegraph(genes[gene])

    ### convert to scipy array
    genes = sp.array([genes[gene] for gene in genes], dtype='object')

    if CFG is not None and CFG['verbose']:
        print >> sys.stderr, "... done"

    if outfile is not None:
        if CFG is not None and CFG['verbose']:
            print >> sys.stderr, "Storing gene structure in %s ..." % outfile

        cPickle.dump(genes, open(outfile, 'w'), -1)

        if CFG is not None and CFG['verbose']:
            print >> sys.stderr, "... done"

    return genes



def init_genes_gff3(infile, CFG=None, outfile=None):
    # genes = init_genes_gff3(infile, CFG=None, outfile=None)

    """This function reads the gff3 input file and returns the information in an
       internal data structure"""

    if CFG is not None and CFG['verbose']:
        print >> sys.stderr, "Parsing annotation from %s ..." % infile
    
    ### initial run to get the transcript to gene mapping
    if CFG is not None and CFG['verbose']:
        print >> sys.stderr, "... init structure"

    trans2gene = dict() ### dict with: keys = transcript IDs, values = gene IDs
    exon2trans = dict() ### dict with: keys = exon IDs, values = transcript IDs
    genes = dict()
    chrms = []

    for line in open(infile, 'r'):
        if line[0] == '#':
            continue
        sl = line.strip().split('\t')
        tags = get_tags_gff3(sl[8])
        if sl[2] in ['mRNA', 'transcript', 'mrna', 'miRNA', 'tRNA', 'snRNA', 'snoRNA', 'ncRNA', 'mRNA_TE_gene', 'rRNA', 'pseudogenic_transcript', 'transposon_fragment']:
            trans2gene[tags['ID']] = tags['Parent']
        elif sl[2] in ['exon', 'Exon']:
            exon2trans[tags['ID']] = tags['Parent']
        elif not 'Parent' in tags:
            try:
                start = int(sl[3]) - 1
            except ValueError:
                start =  -1
            try:
                stop = int(sl[4])
            except ValueError:
                stop = -1
            genes[tags['ID']] = Gene(name=tags['ID'], start=start, stop=stop, chr=sl[0], strand=sl[6], source=sl[1], gene_type=sl[2])
            chrms.append(sl[0])

    chrms = sp.sort(sp.unique(chrms))
    for gid in genes:
        genes[gid].chr_num = sp.where(chrms == genes[gid].chr)[0]

    counter = 1
    for line in open(infile, 'r'):
        if CFG is not None and CFG['verbose'] and counter % 10000 == 0:
            print >> sys.stderr, '.',
        counter += 1        

        if line[0] == '#':
            continue
        sl = line.strip().split('\t')
        
        ### get start and end
        try:
            start = int(sl[3]) - 1
        except ValueError:
            start =  -1
        try:
            stop = int(sl[4])
        except ValueError:
            stop = -1

        ### get tags
        tags = get_tags_gff3(sl[8])

        ### add exons
        if sl[2] == 'exon':
            trans_id = tags['Parent']
            gene_id = trans2gene[trans_id]
            try:
                t_idx = genes[gene_id].transcripts.index(trans_id)
            except ValueError:
                t_idx = len(genes[gene_id].transcripts)
                genes[gene_id].transcripts.append(trans_id)
            genes[gene_id].add_exon(sp.array([int(sl[3]) - 1, int(sl[4]) - 1], dtype='int'), idx=t_idx)

    ### add splicegraphs
    for gene in genes:
        genes[gene].splicegraph = Splicegraph(genes[gene])

    ### convert to scipy array
    genes = sp.array([genes[gene] for gene in genes], dtype='object')

    if CFG is not None and CFG['verbose']:
        print >> sys.stderr, "... done"

    if outfile is not None:
        if CFG is not None and CFG['verbose']:
            print >> sys.stderr, "Storing gene structure in %s ..." % outfile

        cPickle.dump(genes, open(outfile, 'w'), -1)

        if CFG is not None and CFG['verbose']:
            print >> sys.stderr, "... done"

    return genes

def init_regions(fn_bams):
    # regions=init_regions(fn_bams)

    regions = []
    processed = []
    cnt = 0

    if not isinstance(fn_bams, list):
        fn_bams = [fn_bams]

    for i in range(len(fn_bams)):
        if not os.path.exists(fn_bams[i]):
            continue
        else:
            ### load bamfile
            IN = pysam.Samfile(fn_bams[i], 'rb')
            header_info = IN.header['SQ']
            strands = ['+', '-']
            for c in range(len(header_info)):
                if not (header_info[c]['SN'], header_info[c]['LN']) in processed:
                    for s in strands:
                        region = Region()
                        region.chr = header_info[c]['SN']
                        region.chr_num = c
                        region.strand = s
                        region.start = 1
                        region.stop = header_info[c]['LN']
                        region.id = cnt
                        regions.append(region)
                        cnt += 1
                    processed.append((header_info[c]['SN'], header_info[c]['LN']))
            IN.close()
        break

    return sp.array(regions, dtype = 'object')
