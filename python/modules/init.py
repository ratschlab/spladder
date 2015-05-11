import sys
import cPickle
import os
import pysam
import scipy as sp

if __package__ is None:
    __package__ = 'modules'

from .classes.gene import Gene
from .classes.region import Region
from .classes.splicegraph import Splicegraph

def get_tags_gff3(tagline):
    """Extract tags from given tagline"""

    tags = dict()
    for t in tagline.strip(';').split(';'):
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
    # (genes, CFG) = init_genes_gtf(infile, CFG=None, outfile=None)

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
        if sl[2].lower() in ['gene', 'pseudogene', 'unprocessed_pseudogene', 'transposable_element_gene']:
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

    CFG = append_chrms(sp.sort(sp.unique(chrms)), CFG)

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
            except KeyError:
                if 'gene_type' in tags:
                    gene_type = tags['gene_type']
                elif 'gene_biotype' in tags:
                    gene_type = tags['gene_biotype']
                else:
                    gene_type = None
                print >> sys.stderr, 'WARNING: %s does not have gene level information for transcript %s - information has been inferred from tags'  % (infile, trans_id)
                genes[gene_id] = Gene(name=gene_id, start=start, stop=stop, chr=sl[0], strand=sl[6], source=sl[1], gene_type=gene_type)
                t_idx = len(genes[gene_id].transcripts)
                genes[gene_id].transcripts.append(trans_id)
            genes[gene_id].add_exon(sp.array([int(sl[3]) - 1, int(sl[4])], dtype='int'), idx=t_idx)

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

    return (genes, CFG)


def append_chrms(chrms, CFG):
    """Checks if chrm is in lookup table adds it if not"""

    if CFG is None:
        CFG = dict()

    if not 'chrm_lookup' in CFG:
        CFG['chrm_lookup'] = dict()

    for chrm in chrms:
        if not chrm in CFG['chrm_lookup']:
            CFG['chrm_lookup'][chrm] = len(CFG['chrm_lookup'])

    return CFG


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
    genes = dict()
    chrms = []

    for line in open(infile, 'r'):
        if line[0] == '#':
            continue
        sl = line.strip().split('\t')
        tags = get_tags_gff3(sl[8])
        if sl[2] in ['mRNA', 'transcript', 'mrna', 'miRNA', 'tRNA', 'snRNA', 'snoRNA', 'ncRNA', 'rRNA', 'pseudogenic_transcript', 'transposon_fragment', 'mRNA_TE_gene']:
            trans2gene[tags['ID']] = tags['Parent']
        elif sl[2] in ['gene', 'transposable_element_gene', 'pseudogene']:
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


    CFG = append_chrms(sp.sort(sp.unique(chrms)), CFG)

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
        if sl[2] in ['exon', 'pseudogenic_exon']:
            trans_id = tags['Parent']
            gene_id = trans2gene[trans_id]
            try:
                t_idx = genes[gene_id].transcripts.index(trans_id)
            except ValueError:
                t_idx = len(genes[gene_id].transcripts)
                genes[gene_id].transcripts.append(trans_id)
            genes[gene_id].add_exon(sp.array([int(sl[3]) - 1, int(sl[4])], dtype='int'), idx=t_idx)

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

    return (genes, CFG)

def parse_header(header_string):

    hd = dict()

    for line in header_string.strip('\n').split('\n'):
        sl = line.strip().split('\t')
        td = dict([x.split(':', 1) for x in sl[1:]])    
        try:
            hd[sl[0].strip('@')].append(td)
        except KeyError:
            hd[sl[0].strip('@')] = [td]
    return hd

def init_regions(fn_bams, CFG=None):
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
            #header_info = IN.header['SQ']
            header_info = parse_header(IN.text)['SQ']
            
            CFG = append_chrms([x['SN'] for x in header_info], CFG)

            strands = ['+', '-']
            for c in range(len(header_info)):
                if not (header_info[c]['SN'], header_info[c]['LN']) in processed:
                    for s in strands:
                        region = Region()
                        region.chr = header_info[c]['SN']
                        region.chr_num = CFG['chrm_lookup'][region.chr]
                        region.strand = s
                        region.start = 1
                        region.stop = header_info[c]['LN']
                        region.id = cnt
                        regions.append(region)
                        cnt += 1
                    processed.append((header_info[c]['SN'], header_info[c]['LN']))
            IN.close()
        break

    return (sp.array(regions, dtype = 'object'), CFG)
