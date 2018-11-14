import numpy as np
import os
import sys
from collections import namedtuple
from collections import Counter

GeneTable = namedtuple('GeneTable', ['gene_to_first_exon', 'ts_to_exons', 'gene_to_ts'])


def rev_complement(seq):
    rev_dict =  {"A": "T", "T": "A", "C": "G", "G": "C"}  
    return ''.join([rev_dict[_] for _ in seq][::-1])

def leq_strand(coord1, coord2, strand):
    if strand == "+":
        return coord1 <= coord2
    else:
        return coord1 >= coord2

def attribute_item_to_dict(a_item, file_type, feature_type):
    """  From attribute item in annotation file to get corresponding dictionary

    Parameters
    ----------
    a_item: str. attribute item
    file_type: str. Choose from {'gtf', 'gff', 'gff3'}
    feature_type: str. Extract other fields. We only
        consider 'exon', 'mRNA' and 'transcript'

    Returns
    -------
    gtf_dict: dict. store all the necessary data

    """
    gtf_dict = {}
    if file_type == 'gtf':
        attribute_list = a_item.split('; ')
        for attribute_pair in attribute_list:
            pair = attribute_pair.split(' ')
            gtf_dict[pair[0]] = pair[1][1:-1]
    elif file_type == 'gff3':
        attribute_list = a_item.split(';')
        for attribute_pair in attribute_list:
            pair = attribute_pair.split('=')
            gtf_dict[pair[0]] = pair[1]
    elif file_type == 'gff':
        gff_dict = {}
        attribute_list = a_item.split(';')
        for attribute_pair in attribute_list:
            pair = attribute_pair.split('=')
            gff_dict[pair[0]] = pair[1]  # delete "", currently now work on level 2
        if feature_type == 'exon':
            gtf_dict['transcript_id'] = gff_dict['Parent']
        elif feature_type in {'mRNA', 'transcript'}:  # mRNA or transcript
            gtf_dict['gene_id'] = gff_dict['geneID']
            gtf_dict['transcript_id'] = gff_dict['ID']
            gtf_dict['gene_type'] = gff_dict['gene_type']
            gtf_dict['transcript_type'] = gff_dict['transcript_type']

    return gtf_dict

def preprocess_ann(ann_path):
    """ Extract information from annotation file (.gtf, .gff and .gff3)

    Parameters
    ----------
    ann_path: str. Annotation file path

    Returns
    -------
    gene_table: NamedTuple.store the gene-transcript-cds mapping tables derived
        from .gtf file. has attribute ['gene_to_first_exon', 'ts_to_exons', 'gene_to_cds']
    """
    transcript_to_gene_dict = {}    # transcript -> gene id
    gene_to_transcript_dict = {}    # gene_id -> list of transcripts
    transcript_to_exon_dict = {}     # transcript -> list of exons
    transcript_first_exon_dict = {}  # transcript -> first exon of thetranscript
    gene_first_exon_dict = {}        # gene -> list of first exons

    file_type = ann_path.split('.')[-1]

    # collect information from annotation file
    for line in open(ann_path):
        if line[0] == '#':
            continue
        item = line.strip().split('\t')
        feature_type = item[2]
        attribute_item = item[-1]
        attribute_dict = attribute_item_to_dict(attribute_item, file_type, feature_type)
        # store relationship between gene ID and its transcript IDs
        if feature_type in ['transcript', 'mRNA']:
            gene_id = attribute_dict['gene_id']
            transcript_id = attribute_dict['transcript_id']
            # only use protein coding genes
            if attribute_dict['gene_type'] != 'protein_coding' or attribute_dict['transcript_type']  != 'protein_coding':
                continue
            assert (transcript_id not in transcript_to_gene_dict)
            transcript_to_gene_dict[transcript_id] = gene_id
            if gene_id in gene_to_transcript_dict:
                gene_to_transcript_dict[gene_id].add(transcript_id)
            else:
                gene_to_transcript_dict[gene_id] = {transcript_id}

        # Todo python is 0-based while gene annotation file(.gtf, .vcf, .maf) is one based
        elif feature_type == "exon":
            parent_ts = attribute_dict['transcript_id']
            strand_mode = item[6]
            exon_left = int(item[3])-1
            exon_right = int(item[4])-1
            if parent_ts in transcript_to_exon_dict:
                transcript_to_exon_dict[parent_ts].append((exon_left, exon_right))
            else:
                transcript_to_exon_dict[parent_ts] = [(exon_left, exon_right)]
            if strand_mode == "+" :
                trans_start, trans_stop = exon_left, exon_right
            else:
                trans_start, trans_stop = exon_right, exon_left

            # we only consider the start of the whole CoDing Segment
            if parent_ts not in transcript_first_exon_dict or \
               leq_strand(trans_start, transcript_first_exon_dict[parent_ts][0], strand_mode):
                transcript_first_exon_dict[parent_ts] = (trans_start, trans_stop, item)

    # collect first exons for all transcripts of a gene
    for ts_key in transcript_to_gene_dict:
        target_gene = transcript_to_gene_dict[ts_key]
        if target_gene not in gene_first_exon_dict:
            gene_first_exon_dict[target_gene] = []
        if ts_key in transcript_first_exon_dict:
            gene_first_exon_dict[target_gene].append(transcript_first_exon_dict[ts_key])

    # sort list of CDS exons per transcript
    for ts_key in transcript_to_exon_dict:
        transcript_to_exon_dict[ts_key] = sorted(transcript_to_exon_dict[ts_key], key=lambda coordpair: coordpair[0])

    genetable = GeneTable(gene_first_exon_dict, transcript_to_exon_dict, gene_to_transcript_dict)
    return genetable


# steps to add new transcript
# 1. modify annotation_pos.gtf file
# 2. modify create_test_genome file to make tailor-made gene
# 3. delete quick_test_data/spladder
# 4. create new splicegraph with command tools
# 4. run main_immuno.py
def create_test_genome(L, stop_position, data_dir, seed=1):

    ### generate a random DNA string of length L
    np.random.seed(seed)
    map_dict = {0:'A', 1:'G', 2:'C', 3:'T'}
    dna_list = [map_dict[i] for i in np.random.randint(0, 4, L)]
    for pos in stop_position:
        dna_list[pos:pos+3] = 'TAG'

    # modify
    dna_list[69] = 'C'  # remove stop codon
    dna_list[98] = 'G'  # remove stop codon
    dna_list[39] = 'T'  # create stop codon in mutation mode

    # create GT/AG for splicing site
    dna_list[25] = 'G'  # in case the dna_list[24:27] = 'TAG'
    dna_list[26:28] = 'GT'
    dna_list[29:31] = 'GT'
    dna_list[36:38] = 'AG'
    dna_list[50:52] = 'GT'
    dna_list[59:61] = 'AG'
    dna_list[64:66] = 'AG'
    dna_list[75:77] = 'GT'
    dna_list[85:87] = 'AG'

    pos_dna = ''.join(dna_list)
    neg_dna = rev_complement(pos_dna)

    pos_seq_file = os.path.join(data_dir, 'genome_pos.fa')
    neg_seq_file = os.path.join(data_dir, 'genome_neg.fa')
    f_pos = open(pos_seq_file, 'w')
    f_neg = open(neg_seq_file, 'w')
    f_pos.write('>chr1'+'\n')
    f_pos.write(pos_dna)
    f_pos.close()
    f_neg.write('>chr1'+'\n')
    f_neg.write(neg_dna)
    f_neg.close()
    return pos_dna


def create_neg_from_pos(data_dir, L):
    posgtf_file = os.path.join(data_dir, 'annotation_pos.gtf')
    neggtf_file = os.path.join(data_dir, 'annotation_neg.gtf')

    f = open(posgtf_file, 'r')
    lines = f.readlines()
    new_line_list = []
    for line in lines:
        item = line.split('\t')
        start_pos = L+1-int(item[3])
        end_pos = L+1-int(item[4])
        item[3] = str(end_pos)
        item[4] = str(start_pos)
        item[6] = '-'
        new_line = '\t'.join(item)
        new_line_list.append(new_line)
    f = open(neggtf_file,'w')
    f.writelines(new_line_list)
    f.close()


def create_bam_file(data_dir, seed=1):

    np.random.seed(seed)
    posref_path = os.path.join(data_dir, 'genome_pos.fa')
    with open(posref_path, 'r') as fpos:
        posref_seq = fpos.readlines()[1]
    posann_path = os.path.join(data_dir, 'annotation_pos.gtf')
    pos_genetable = preprocess_ann(posann_path)

    ### define the parameters for the simulation
    #expr = [10, 0, 25, 15, 10]
    expr = {'TRANS1.1':10,
            'TRANS1.2':0,
            'TRANS1.3':25,
            'TRANS1.4':15,
            'TRANS1.5':10}
    reads_L = 15
    cover = 100

    expr_file = os.path.join(data_dir, 'expr_ts.txt')
    f_expr = open(expr_file, 'w')
    total_reads_list = []
    ts_dict = {}
    expr_list = []
   
    for i, key_value_pair in enumerate(pos_genetable.ts_to_exons.items()):
        ts_name = key_value_pair[0]
        pos_exon_coord = key_value_pair[1]
        ts = ''
        for coord_pair in pos_exon_coord:
            start_pos = coord_pair[0]
            stop_pos = coord_pair[1]
            ts += posref_seq[start_pos:stop_pos+1]
        N = expr[ts_name] * cover
        pos = np.random.randint(0, len(ts) - reads_L, N)
        pos_counter = Counter(pos)
        expr_arr = np.zeros(len(ts))
        expr_list.append(expr_arr)
        ts_dict[ts_name] = ts
        for ipos in pos_counter.keys():
            expr_arr[ipos:ipos+reads_L] += pos_counter[ipos]
        expr_arr = [str(int(i_expr)) for i_expr in expr_arr]
        new_line = ts_name + '\t' + '\t'.join(expr_arr) + '\n'
        f_expr.write(new_line)
        reads_list = [ts[ipos:ipos+reads_L] for ipos in pos]
        total_reads_list.extend(reads_list)
    f_expr.close()
    output_dir = os.path.join(data_dir, '%s_%i.fq' % (testname, seed))
    f = open(output_dir, 'w')
    total_num = len(total_reads_list)
    for i in range(total_num-1, -1, -1):
        line1 = '@%s_%i_' % (testname, seed) + str((i+1)*2) + '/1'+'\n'
        f.write(line1)
        line2 = total_reads_list[i]+'\n'
        f.write(line2)
        line3 = '+\n'
        f.write(line3)
        line4 = reads_L*'G'+'\n'
        f.write(line4)


if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s <out_dir> [<test_name> (test1)] [<n_bams> (1)]\n' % sys.argv[0])
    sys.exit(1)
data_dir = sys.argv[1]
testname = 'test1'
if len(sys.argv) > 2:
    testname = sys.argv[2]
nbams = 1
if len(sys.argv) > 3:
    nbams = int(sys.argv[3])
create_test_genome(150, [73, 102], data_dir)
print("Successfully created test genome")
create_neg_from_pos(data_dir, 150)
print("Successfully created gtf file")
for i in range(nbams):
    create_bam_file(data_dir, seed=i+1)
print("Successfully created bam segments")
