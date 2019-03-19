import sys
import pickle
import os
import pysam
import scipy as sp
import h5py
import re

if __package__ is None:
    __package__ = 'modules'

from .classes.region import Region

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


def init_genes_gtf(options):

    """This function reads the gtf input file and returns the information in an
       internal data structure"""

    from .classes.gene import Gene
    from .classes.splicegraph import Splicegraph

    if options.verbose:
        print("Parsing annotation from %s ..." % options.annotation, file=sys.stderr)
    
    ### initial run to get the transcript to gene mapping
    if options.verbose:
        print("... init structure", file=sys.stderr)

    genes = dict()
    chrms = []

    for line in open(options.annotation, 'r'):
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
            if 'gene_type' in tags:
                gene_type = tags['gene_type']
            elif 'gene_biotype' in tags:
                gene_type = tags['gene_biotype']
            else:
                gene_type = None
            genes[tags['gene_id']] = Gene(name=tags['gene_id'], start=start, stop=stop, chr=sl[0], strand=sl[6], source=sl[1], gene_type=gene_type)
            chrms.append(sl[0])

    options = append_chrms(sp.sort(sp.unique(chrms)), options)

    counter = 1
    warn_infer_count = 0
    inferred_genes = False
    for line in open(options.annotation, 'r'):
        if options.verbose and counter % 10000 == 0:
            sys.stdout.write('.')
            if counter % 100000 == 0:
                sys.stdout.write(' %i lines processed\n' % counter)
            sys.stdout.flush()
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

                warn_infer_count += 1
                if warn_infer_count < 5:
                    print('WARNING: %s does not have gene level information for transcript %s - information has been inferred from tags'  % (options.annotation, trans_id), file=sys.stderr)
                elif warn_infer_count == 5:
                    print('WARNING: too many warnings for inferred tags', file=sys.stderr)
                    
                genes[gene_id] = Gene(name=gene_id, start=start, stop=stop, chr=sl[0], strand=sl[6], source=sl[1], gene_type=gene_type)
                t_idx = len(genes[gene_id].transcripts)
                genes[gene_id].transcripts.append(trans_id)
                inferred_genes = True
            genes[gene_id].add_exon(sp.array([int(sl[3]) - 1, int(sl[4])], dtype='int'), idx=t_idx)

    ### post-process in case we have inferred genes
    if warn_infer_count >= 5:
        print('\nWARNING: a total of %i cases had no gene level information annotated - information has been inferred from tags' % warn_infer_count, file=sys.stderr)
    if inferred_genes:
        for gene in genes:
            if len(genes[gene].exons) == 0:
                continue
            genes[gene].start = min([x.min() for x in genes[gene].exons])
            genes[gene].stop = max([x.max() for x in genes[gene].exons])

    ### add splicegraphs
    for gene in genes:
        genes[gene].splicegraph = Splicegraph(genes[gene])

    ### convert to scipy array
    genes = sp.array([genes[gene] for gene in genes], dtype='object')

    ### check for consistency and validity of genes
    genes = check_annotation(options, genes)

    if options.verbose:
        print("... done", file=sys.stderr)

    if options.verbose:
        print("Storing gene structure in %s.pickle ..." % options.annotation, file=sys.stderr)

    pickle.dump(genes, open(options.annotation + '.pickle', 'wb'), -1)

    if options.verbose:
        print("... done", file=sys.stderr)

    return (genes, options)


def append_chrms(chrms, options):
    """Checks if chrm is in lookup table adds it if not"""

    if not hasattr(options, 'chrm_lookup'):
        options.chrm_lookup = dict()

    for chrm in chrms:
        if not chrm in options.chrm_lookup:
            options.chrm_lookup[chrm] = len(options.chrm_lookup)

    return options


def init_genes_gff3(options):

    """This function reads the gff3 input file and returns the information in an
       internal data structure"""

    from .classes.gene import Gene
    from .classes.splicegraph import Splicegraph

    ### initial run to get the transcript to gene mapping
    if options.verbose:
        print("Parsing annotation from %s ..." % options.annotation, file=sys.stderr)
        print("... init structure", file=sys.stderr)

    trans2gene = dict() ### dict with: keys = transcript IDs, values = gene IDs
    genes = dict()
    chrms = []

    for line in open(options.annotation, 'r'):
        if line[0] == '#':
            continue
        sl = line.strip().split('\t')
        tags = get_tags_gff3(sl[8])
        if sl[2].lower() in ['chromosome', 'contig', 'supercontig']:
            continue
        ### middle layer
        if 'ID' in tags and 'Parent' in tags:
            trans2gene[tags['ID']] = tags['Parent']
        ### root layer
        elif 'ID' in tags and not 'Parent' in tags:
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

    options = append_chrms(sp.sort(sp.unique(chrms)), options)

    counter = 1
    for line in open(options.annotation, 'r'):

        if options.verbose and counter % 10000 == 0:
            sys.stdout.write('.')
            if counter % 100000 == 0:
                sys.stdout.write(' %i lines processed\n' % counter)
            sys.stdout.flush()
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

    ### check for consistency and validity of genes
    genes = check_annotation(options, genes)

    if options.verbose:
        print("... done", file=sys.stderr)

    if options.verbose:
        print("Storing gene structure in %s.pickle ..." % options.annotation, file=sys.stderr)

    pickle.dump(genes, open(options.annotation + '.pickle', 'wb'), -1)

    if options.verbose:
        print("... done", file=sys.stderr)

    return (genes, options)


def parse_header(header_string):

    hd = dict()

    for line in header_string.strip('\n').split('\n'):
        sl = line.strip().split('\t')
        ### ignore comment lines
        if sl[0] == '@CO':
            continue
        td = dict([x.split(':', 1) for x in sl[1:] if ':' in x])    
        try:
            hd[sl[0].strip('@')].append(td)
        except KeyError:
            hd[sl[0].strip('@')] = [td]
    return hd

def init_regions(fn_bams, conf, options=None, sparse_bam=False):

    regions = []
    processed = []
    cnt = 0

    if not isinstance(fn_bams, list):
        fn_bams = [fn_bams]

    for i in range(len(fn_bams)):
        if not os.path.exists(fn_bams[i]):
            continue
        else:
            if sparse_bam:
                ### load bamfile
                IN = h5py.File(re.sub(r'.bam$', '', fn_bams[i]) + '.conf_%i' % conf + '.filt.hdf5', 'r')
                
                strands = ['+', '-']
                for k in IN:
                    if not k.endswith('_reads_shp'):
                        continue
                    chrm = re.sub(r'_reads_shp$', '', k)
                    lens = IN[k][:][1]
                    if not (chrm, lens) in processed:
                        for s in strands:
                            region = Region()
                            region.chr = chrm
                            options = append_chrms([chrm], options)
                            region.chr_num = options.chrm_lookup[region.chr]
                            region.strand = s
                            region.start = 1
                            region.stop = lens
                            region.id = cnt
                            regions.append(region)
                            cnt += 1
                        processed.append((chrm, lens))
            else:
                ### load bamfile
                IN = pysam.Samfile(fn_bams[i], 'rb')
                #header_info = IN.header['SQ']
                header_info = parse_header(IN.text)['SQ']
                
                options = append_chrms([x['SN'] for x in header_info], options)

                strands = ['+', '-']
                for c in range(len(header_info)):
                    if not (header_info[c]['SN'], header_info[c]['LN']) in processed:
                        for s in strands:
                            region = Region()
                            region.chr = header_info[c]['SN']
                            region.chr_num = options.chrm_lookup[region.chr]
                            region.strand = s
                            region.start = 1
                            region.stop = header_info[c]['LN']
                            region.id = cnt
                            regions.append(region)
                            cnt += 1
                        processed.append((header_info[c]['SN'], header_info[c]['LN']))
                IN.close()
        break

    return (sp.array(regions, dtype = 'object'), options)


def check_annotation(options, genes):
    
    if options.verbose:
        print('\n... checking annotation')

    ### check whether genes have no exons annotated
    rm_ids = []
    for gene in genes:
        if len(gene.exons) == 0:
            rm_ids.append(gene.name)
    if len(rm_ids) > 0:
        print('WARNING: removing %i genes from given annotation that had no exons annotated:' % len(rm_ids), file=sys.stderr)
        print('list of excluded genes written to: %s' % (options.annotation + '.genes_excluded_no_exons'), file=sys.stderr)
        sp.savetxt(options.annotation + '.genes_excluded_no_exons', rm_ids, fmt='%s', delimiter='\t')
        gene_names = sp.array([x.name for x in genes], dtype='str')
        k_idx = sp.where(~sp.in1d(gene_names, rm_ids))[0]
        genes = genes[k_idx]

    ### check whether we run unstranded analysis and have to exclude overlapping gene annotations
    if options.filter_overlap_genes:
        rm_ids = []
        chrms = sp.array([x.chr for x in genes])
        strands = sp.array([x.strand for x in genes])
        starts = sp.array([x.start for x in genes], dtype='int')
        stops = sp.array([x.stop for x in genes], dtype='int')
        for c in sp.unique(chrms):
            for s in sp.unique(strands):
                c_idx = sp.where((chrms == c) & (strands == s))[0]
                for i in c_idx:
                    if sp.sum((starts[i] <= stops[c_idx]) & (stops[i] >= starts[c_idx])) > 1:
                        rm_ids.append(genes[i].name)
        if len(rm_ids) > 0:
            rm_ids = sp.unique(rm_ids)
            print('WARNING: removing %i genes from given annotation that overlap to each other:' % rm_ids.shape[0], file=sys.stderr)
            print('list of excluded genes written to: %s' % (options.annotation + '.genes_excluded_gene_overlap'), file=sys.stderr)
            sp.savetxt(options.annotation + '.genes_excluded_gene_overlap', rm_ids, fmt='%s', delimiter='\t')
            gene_names = sp.array([x.name for x in genes], dtype='str')
            k_idx = sp.where(~sp.in1d(gene_names, rm_ids))[0]
            genes = genes[k_idx]

    ### check whether exons are part of multiple genes
    if options.filter_overlap_exons:
        exon_map = dict()
        for i, g in enumerate(genes):
            for t in range(len(g.exons)):
                for e in range(g.exons[t].shape[0]):
                    k = ':'.join([g.chr, str(g.exons[t][e, 0]), str(g.exons[t][e, 1])])
                    if k in exon_map:
                        exon_map[k].append(g.name)
                    else:
                        exon_map[k] = [g.name]
        rm_ids = []
        for exon in exon_map:
            if sp.unique(exon_map[exon]).shape[0] > 1:
                rm_ids.extend(exon_map[exon])
        if len(rm_ids) > 0:
            rm_ids = sp.unique(rm_ids)
            print('WARNING: removing %i genes from given annotation that share exact exon coordinates:' % rm_ids.shape[0], file=sys.stderr)
            print('list of excluded exons written to: %s' % (options.annotation + '.genes_excluded_exon_shared'), file=sys.stderr)
            sp.savetxt(options.annotation + '.genes_excluded_exon_shared', rm_ids, fmt='%s', delimiter='\t')
            gene_names = sp.array([x.name for x in genes], dtype='str')
            k_idx = sp.where(~sp.in1d(gene_names, rm_ids))[0]
            genes = genes[k_idx]

    ### check whether exons within the same transcript overlap
    if options.filter_overlap_transcripts:
        rm_ids = []
        for i, g in enumerate(genes):
            for t in range(len(g.exons)):
                for e in range(g.exons[t].shape[0] - 1):
                    if sp.any(g.exons[t][e+1:, 0] < g.exons[t][e, 1]):
                        rm_ids.append(g.name)
        if len(rm_ids) > 0:
            rm_ids = sp.unique(rm_ids)
            print('WARNING: removing %i genes from given annotation that have at least one transcript with overlapping exons.' % rm_ids.shape[0], file=sys.stderr)
            print('list of excluded genes written to: %s' % (options.annotation + '.genes_excluded_exon_overlap'), file=sys.stderr)
            sp.savetxt(options.annotation + '.genes_excluded_exon_overlap', rm_ids, fmt='%s', delimiter='\t')
            gene_names = sp.array([x.name for x in genes], dtype='str')
            k_idx = sp.where(~sp.in1d(gene_names, rm_ids))[0]
            genes = genes[k_idx]

    ### do we have any genes left?
    if genes.shape[0] == 0:
        print('\nERROR: there are no valid genes left in the input. Please verify correctnes of input annotation.\n', file=sys.stderr)
        sys.exit(1)

    return genes


