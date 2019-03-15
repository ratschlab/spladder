import scipy as sp
import warnings
import intervaltree as it

if __package__ is None:
    __package__ = 'modules'

from .reads import *

def make_introns_feasible(introns, genes, options):

    tmp1 = sp.array([x.shape[0] for x in introns[:, 0]])
    tmp2 = sp.array([x.shape[0] for x in introns[:, 1]])

    if options.logfile == '-':
        fd_log = sys.stdout
    else:
        fd_log = open(options.logfile, 'w')
    
    unfeas = sp.where((tmp1 > 200) | (tmp2 > 200))[0]
    print('found %i unfeasible genes' % unfeas.shape[0], file=fd_log)

    while unfeas.shape[0] > 0:
        ### make filter more stringent
        options.read_filter['exon_len'] = min(36, options.read_filter['exon_len'] + 4)
        options.read_filter['mincount'] = 2 * options.read_filter['mincount']
        options.read_filter['mismatch'] = max(options.read_filter['mismatch'] - 1, 0)

        ### get new intron counts
        tmp_introns = get_intron_list(genes[unfeas], options)
        introns[unfeas, :] = tmp_introns

        ### still unfeasible?
        tmp1 = sp.array([x.shape[0] for x in introns[:, 0]])
        tmp2 = sp.array([x.shape[0] for x in introns[:, 1]])

        still_unfeas = sp.where((tmp1 > 200) | (tmp2 > 200))[0]
        idx = sp.where(~sp.in1d(unfeas, still_unfeas))[0]

        for i in unfeas[idx]:
            print('[feasibility] set criteria for gene %s to: min_ex %i, min_conf %i, max_mism %i' % (genes[i].name, options.read_filter['exon_len'], options.read_filter['mincount'], options.read_filter['mismatch']), file=fd_log)
        unfeas = still_unfeas;

    if fd_log != sys.stdout:
        fd_log.close()

    return introns


### remove introns overlapping to more than one gene
def filter_introns(introns, genes, options):
    
    ### build interval trees of all genes starts and ends
    chrms = sp.array([_.strand for _ in genes])
    strands = sp.array([_.chr for _ in genes])
    gene_trees = dict()
    for c in sp.unique(chrms):
        for s in sp.unique(strands):
            gene_trees[(c, s)] = it.IntervalTree()
            c_idx = sp.where((chrms == c) & (strands == s))[0]
            for i in c_idx:
                gene_trees[(c, s)][genes[i].start:genes[i].stop] = i

    ### match all introns agains trees and remove elements overlapping
    ### more than one gene on the same chr/strand
    cnt_tot = 0
    cnt_rem = 0
    strand_list = ['+', '-']
    offset = options.intron_edges['append_new_terminal_exons_len']
    for si, s in enumerate(strand_list):
        for i in range(introns.shape[0]):
            if introns[i, si].shape[0] == 0:
                continue
            k_idx = []
            cnt_tot += introns[i, si].shape[0]
            for j in range(introns[i, si].shape[0]):
                if len(gene_trees[(s, genes[i].chr)].overlap(introns[i, si][j, 0] - offset, introns[i, si][j, 1] + offset)) == 1:
                    k_idx.append(j)
            if len(k_idx) < introns[i, si].shape[0]:
                cnt_rem += (introns[i, si].shape[0] - len(k_idx))
                introns[i, si] = introns[i, si][k_idx, :]
    print('removed %i of %i (%.2f percent) introns overlapping to no or multiple genes' % (cnt_rem, cnt_tot, cnt_rem / float(max(cnt_tot, 1)) * 100))

    return introns


### determine count output file
def get_filename(which, options, sample_idx=None):
    """This function returns a filename generated from the current configuration"""

    ### init any tags
    prune_tag = ''
    if options.do_prune:
        prune_tag = '_pruned'
    validate_tag = ''
    if options.validate_sg:
        validate_tag = '.validated'

    ### iterate over return file types    
    if which in ['fn_count_in', 'fn_count_out']:
        if options.spladderfile == '-':
            if options.merge == 'single':
                fname = os.path.join(options.outdir, 'spladder', 'genes_graph_conf%i.%s%s.pickle' % (options.confidence, options.samples[sample_idx], prune_tag))
            else:
                if (options.qmode == 'single' and which != 'fn_count_in') or (options.qmode == 'collect' and which == 'fn_count_in'):
                    fname = os.path.join(options.outdir, 'spladder', 'genes_graph_conf%i.%s.%s%s%s.pickle' % (options.confidence, options.merge, options.samples[sample_idx], prune_tag, validate_tag))
                else:
                    fname = os.path.join(options.outdir, 'spladder', 'genes_graph_conf%i.%s%s%s.pickle' % (options.confidence, options.merge, prune_tag, validate_tag))
        else:
            fname = options.spladderfile
        
        if which == 'fn_count_in':
            if options.qmode == 'collect':
                return fname.replace('.pickle', '') + '.count.hdf5'
            else:
                return fname
        elif which == 'fn_count_out':
            return fname.replace('.pickle', '') + '.count.hdf5'
    elif which == 'fn_out_merge':
        if options.merge == 'merge_graphs':
            return os.path.join(options.outdir, 'spladder', 'genes_graph_conf%i.%s%s.pickle' % (options.confidence, options.merge, prune_tag))
        else:
            return ''
    elif which == 'fn_out_merge_val':
        return os.path.join(options.outdir, 'spladder', 'genes_graph_conf%i.%s%s%s.pickle' % (options.confidence, options.merge, validate_tag, prune_tag))

def compute_psi(counts, event_type, options):
    
    ### collect count data based on event type
    if event_type == 'exon_skip':
        #a = counts[:, 4] + counts[:, 5]
        #b = 2 * counts[:, 6]
        a = sp.c_[counts[:, 4], counts[:, 5]].min(axis=1)
        b = counts[:, 6]
    elif event_type == 'intron_retention':
        a = counts[:, 1] # intron cov
        b = counts[:, 4] # intron conf
    elif event_type in ['alt_3prime', 'alt_5prime']:
        a = counts[:, 3] # intron1 conf
        b = counts[:, 4] # intron2 conf
    elif event_type == 'mutex_exons':
        a = counts[:, 5] + counts[:, 7] # exon_pre_exon1_conf + exon1_exon_aft_conf
        b = counts[:, 6] + counts[:, 8] # exon_pre_exon2_conf + exon2_exon_aft_conf
    elif event_type == 'mult_exon_skip':
        a = counts[:, 4] + counts[:, 5] + counts[:, 7] # exon_pre_exon_conf + exon_exon_aft_conf + sum_inner_exon_conf
        b = (counts[:, 8] + 1) * counts[:, 6] # (num_inner_exon + 1) * exon_pre_exon_aft_conf
    else:
        raise Exception('Unknown event type: %s' % event_type)

    ### compute psi - catch div by 0 warning
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        psi = a / (a + b)  

    ### filter for sufficient read support
    n_idx = sp.where((a + b) < options.psi_min_reads)
    psi[n_idx] = sp.nan

    return (psi, a, b)


def log_progress(idx, total, bins=50):
    
    global TIME0

    #binsize = max(total / bins, 1)
    binsize = bins / total
    time1 = time.time()
    if idx == 0:
        TIME0 = time1
    progress = int((idx + 1) * binsize)
    sys.stdout.write('\r[' + ('#' * progress) + (' ' * (bins - progress)) + ']' + ' %i / %i (%.0f%%)' % (idx, total, float(idx) / max(total, 1) * 100) + ' - took %i sec (ETA: %i sec)' % (time1 - TIME0, int((bins - progress) * float(time1 - TIME0) / max(progress, 1))))
    sys.stdout.flush()

def codeUTF8(s):
    return s.view(sp.chararray).encode('utf-8')

def decodeUTF8(s):
    return s.view(sp.chararray).decode('utf-8')

