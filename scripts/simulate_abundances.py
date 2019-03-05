import sys
import scipy as sp
import numpy.random as npr
npr.seed(23)

if len(sys.argv) < 2:
    print('Usage: %s <blocks>' % sys.argv[0])
    sys.exit(1)
blockfile = sys.argv[1]

N = 20
split = 0.5

### pick gene to be differentially expressed
### here, the transcript ratio is 10% vs 90%, whereas in others
### expression is uniformlly distributed
diff_genes = ['gene1']

### parse blocks
genes = dict()
for line in open(blockfile, 'r'):
    sl = line.strip().split('\t')
    if not sl[0] in genes:
        genes[sl[0]] = [sl[1]]
        continue
    genes[sl[0]].append(sl[1])

### simulate libary sizes
libsize = sp.array([0.5 * npr.randint(1,10) for _ in range(N)])

### simulate gene expressioned abundances (total transcript counts per gene)
gene_exp = dict()
for gene in genes:
    gene_exp[gene] = npr.randint(20,50)

### simulate individual transcript counts
for gene in genes:
    for tidx, tr in enumerate(genes[gene]):
        if gene in diff_genes:
            if tidx == 0:
                exp1 = gene_exp[gene] * 0.8 * libsize[:int(N*split)]
                exp2 = gene_exp[gene] * 0.2 * libsize[int(N*split):]
            else:
                exp1 = 0.2 * gene_exp[gene] / (len(genes[gene]) - 1) * libsize[:int(N*split)]
                exp2 = 0.8 * gene_exp[gene] / (len(genes[gene]) - 1) * libsize[int(N*split):]
            exp = sp.r_[exp1, exp2] 
        else:
            exp = gene_exp[gene] / len(genes[gene]) * libsize
        exp *= (npr.rand(exp.shape[0]) + 0.5)
        print('\t'.join([gene, tr, '\t'.join(exp.astype('int').astype('str'))]))
