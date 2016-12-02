import scipy as sp
import scipy.io as scio
import statsmodels.api as sm
import statsmodels.sandbox as sms
import h5py
import sys
import os
import pdb
import cPickle
import warnings
import time
import datetime
from scipy.optimize import minimize_scalar
from scipy.special import polygamma 
from scipy.stats import chi2,nbinom
import numpy.random as npr

from modules.classes.gene import Gene
import modules.alt_splice.quantify as quantify
import modules.testing.likelihood as likelihood
import modules.settings as settings 
import modules.viz.diagnose as plot

import multiprocessing as mp 
import signal as sig

from modules.helpers import log_progress 

TIME0 = time.time()

class Dummy():
    """Dummy class to mimic matlab structs"""
    pass


def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'MANDATORY')
    required.add_option('-o', '--outdir', dest='outdir', metavar='DIR', help='spladder output directory', default='-')
    required.add_option('-a', '--conditionA', dest='conditionA', metavar='idA1,idA2,idA3,...', help='comma separated list of samples files for condition A', default='-')
    required.add_option('-b', '--conditionB', dest='conditionB', metavar='idB1,idB2,idB3,...', help='comma separated list of samples files for condition B', default='-')
    input = OptionGroup(parser, 'INPUT OPTIONS')
    input.add_option('-n', '--readlen', dest='readlen', metavar='INT', type='int', help='read length [50]', default=50)
    input.add_option('-c', '--confidence', dest='confidence', metavar='INT', type='int', help='confidence level (0 lowest to 3 highest) [3]', default=3)
    input.add_option('-t', '--event_types', dest='event_types', metavar='y|n', help='list of alternative splicing events to be tested [exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip]', default='exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip')
    input.add_option('-M', '--merge_strat', dest='merge', metavar='<STRAT>', help='merge strategy, where <STRAT> is one of: merge_bams, merge_graphs, merge_all [merge_graphs]', default='merge_graphs')
    input.add_option('-V', '--validate_sg', dest='validate_sg', metavar='y|n', help='validate splice graph [n]', default='n')
    input.add_option('-m', '--matlab', dest='matlab', metavar='y|n', help='input data was generated with matlab version [n]', default='n')
    input.add_option('-s', '--subset_samples', dest='subset_samples', metavar='y|n', help='gene expression counting will be only done on the tested subset of samples [n]', default='n')
    testing = OptionGroup(parser, 'TESTING OPTIONS')
    testing.add_option('-C', '--correction', dest='correction', metavar='STR', help='method for multiple testing correction (BH, Bonferroni, Holm, Hochberg, BY, TSBH) [BH]', default='BH')
    testing.add_option('-0', '--max_zero_frac', dest='max_0_frac', metavar='FLOAT', type='float', help='max fraction of 0 values per event isoform quantification over all tested samples [0.5]', default=0.5)
    testing.add_option('-i', '--min_count', dest='min_count', metavar='INT', help='min read count sum over all samples for an event isoform to be tested [10]', default=10)
    output = OptionGroup(parser, 'OUTPUT OPTIONS')
    output.add_option('-v', '--verbose', dest='verbose', metavar='y|n', help='verbosity', default='n')
    output.add_option('-d', '--debug', dest='debug', metavar='y|n', help='use debug mode [n]', default='n')
    output.add_option('-D', '--diagnose_plots', dest='diagnose_plots', metavar='y|n', help='generate diagnose plots [n]', default='n')
    output.add_option('--timestamp', dest='timestamp', metavar='y|n', help='add timestamp to output directory [n]', default='n')
    output.add_option('--labelA', dest='labelA', metavar='STRING', help='label for condition A (used for output naming)', default='condA')
    output.add_option('--labelB', dest='labelB', metavar='STRING', help='label for condition B (used for output naming)', default='condB')
    experimental = OptionGroup(parser, 'EXPERIMENTAL - BETA STATE')
    experimental.add_option('', '--parallel', dest='parallel', metavar='<INT>', type='int', help='use multiple processors [1]', default=1)
    experimental.add_option('', '--non_alt_norm', dest='non_alt_norm', metavar='y|n', help='only use non alternative exon segments for gene expression counting [n]', default='n')
    parser.add_option_group(required)
    parser.add_option_group(input)
    parser.add_option_group(testing)
    parser.add_option_group(output)
    parser.add_option_group(experimental)

    (options, args) = parser.parse_args()

    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)

    options.parser = parser
    return options


def get_non_alt_seg_ids_matlab(gene):

    tmp = sp.ones((gene.segmentgraph[0, 2].shape[0],), dtype='bool')
    for i in range(gene.segmentgraph[0, 2].shape[0] - 1):
        ### get index of last acceptor
        idx = sp.where(gene.segmentgraph[0, 2][i, i + 1:])[0]
        ### mask all segments between current segment and acceptor
        if idx.shape[0] > 0:
            tmp[i + 1:idx[-1] + i + 1] = 0

    return sp.where(tmp)[0]

def get_gene_expression(CFG, fn_out=None, strain_subset=None):

    if CFG['verbose']:
        sys.stdout.write('Quantifying gene expression ...\n')
    
    ### load gene information
    if CFG['is_matlab']:
        genes = scio.loadmat(CFG['fname_genes'], struct_as_record=False)['genes'][0, :]
        numgenes = len(genes)
    else:
        genes = cPickle.load(open(CFG['fname_genes'], 'r'))[0]
        numgenes = genes.shape[0]

    ### open hdf5 file containing graph count information
    IN = h5py.File(CFG['fname_count_in'], 'r')
    strains = IN['strains'][:].astype('str')
    ### sort by strain ID
    s_idx = sp.argsort(strains)
    if strain_subset is None:
        strain_idx = s_idx.copy()
    else:
        strain_idx = s_idx[sp.in1d(strains[s_idx], strain_subset)]
    gene_counts = sp.zeros((numgenes, strain_idx.shape[0]), dtype='float')
    gene_names = sp.array([x.name for x in genes], dtype='str')

    if CFG['is_matlab']:
        seg_lens = IN['seg_len'][:, 0]
        gene_ids_segs = IN['gene_ids_segs'][0, :].astype('int') - 1
    else:
        seg_lens = IN['seg_len'][:]
        gene_ids_segs = IN['gene_ids_segs'][:].astype('int')

    ### no longer assume that the gene_ids_segs are sorted by gene ID
    s_idx = sp.argsort(gene_ids_segs[:, 0], kind='mergesort')
    _, u_idx = sp.unique(gene_ids_segs[s_idx, 0], return_index=True)
    s_idx = s_idx[u_idx]

    ### iterate over genes
    for gidx, iidx in enumerate(s_idx):

        if CFG['verbose']:  
            log_progress(gidx, numgenes, 100)
            
        ### get idx of non alternative segments
        if CFG['is_matlab']:
            non_alt_idx = get_non_alt_seg_ids_matlab(genes[gidx])
            seg_idx = sp.arange(iidx, iidx + genes[gidx].segmentgraph[0, 2].shape[0])
            if len(seg_idx) == 0:
                continue
        else:
            non_alt_idx = genes[gidx].get_non_alt_seg_ids()
            seg_idx = sp.arange(iidx, iidx + genes[gidx].segmentgraph.seg_edges.shape[0])

        gene_idx = gene_ids_segs[seg_idx]
        if len(gene_idx.shape) > 0:
            gene_idx = gene_idx[0]

        if CFG['is_matlab']:
            assert(IN['gene_names'][gene_idx] == genes[gidx].name)
        else:
            assert(IN['gene_names'][:][gene_idx] == genes[gidx].name)
        assert(genes[gidx].name == gene_names[gidx])

        if CFG['non_alt_norm']:
            seg_idx = seg_idx[non_alt_idx]

        ### compute gene expression as the read count over all non alternative segments
        if CFG['is_matlab']:
            #gene_counts[gidx, :] = sp.dot(IN['segments'][:, seg_idx], IN['seg_len'][seg_idx, 0]) / sp.sum(IN['seg_len'][seg_idx, 0])
            gene_counts[gidx, :] = sp.dot(IN['segments'][:, seg_idx][strain_idx], seg_lens[seg_idx]) / CFG['read_length']
            #seg_offset += genes[gidx].segmentgraph[0, 2].shape[0]
        else:
            #gene_counts[gidx, :] = sp.dot(IN['segments'][seg_idx, :].T, IN['seg_len'][:][seg_idx]) / sp.sum(IN['seg_len'][:][seg_idx])
            if seg_idx.shape[0] > 1:
                gene_counts[gidx, :] = sp.squeeze(sp.dot(IN['segments'][seg_idx, :][:, strain_idx].T, seg_lens[seg_idx])) / CFG['read_length']
            else:
                gene_counts[gidx, :] = IN['segments'][seg_idx, :][strain_idx] * seg_lens[seg_idx] / CFG['read_length']
            #seg_offset += genes[gidx].segmentgraph.seg_edges.shape[0]

    IN.close()

    if CFG['verbose']:
        sys.stdout.write('\n... done.\n')

    ### write results to hdf5
    if fn_out is not None:
        OUT = h5py.File(fn_out, 'w')
        OUT.create_dataset(name='strains', data=strains[strain_idx])
        OUT.create_dataset(name='genes', data=gene_names)
        OUT.create_dataset(name='raw_count', data=gene_counts, compression="gzip")
        OUT.close()

    return (gene_counts, strains[strain_idx], gene_names)


def get_size_factors(gene_counts, CFG):

    if CFG['verbose']:
        print 'Estimating size factors'

    ### take geometric mean of counts
    gmean = sp.exp(sp.mean(sp.log(gene_counts + 1), axis=1))

    size_factors = []
    for i in xrange(gene_counts.shape[1]):
        idx = gene_counts[:, i] > 0
        size_factors.append(sp.median(gene_counts[idx, i] / gmean[idx]))

    size_factors = sp.array(size_factors, dtype='float')

    return size_factors


def re_quantify_events(CFG):
    """This is more a legacy function for testing that requantifies events on a given graph"""

    ### load events
    if CFG['is_matlab']:
        if CFG['fname_events'].endswith('mat'):
            try:
                ev = scio.loadmat(CFG['fname_events'], struct_as_record=False)['events_all'][0, :]
            except NotImplementedError:
                print >> sys.stderr, 'The event file in matlab format is too big to be loaded with python correctly. Please use the script events_to_hdf5.m in the matlab/src part of SplAdder to convert your event file to HDF5 and use it here instead.'
                sys.exit(1)
        else:
            ev = []
            IN = h5py.File(CFG['fname_events'], 'r', driver='core')
            shp = IN['chr'].shape[0]
            for i in range(shp):
                if CFG['verbose']:
                    sys.stderr.write('.')
                    sys.stderr.flush()
                    if (i + 1) % 100 == 0:
                        sys.stderr.write('%i/%i\n' % (i + 1, shp + 1))
                tmp = Dummy()
                for k in IN.keys():
                    if IN[k].shape[0] == shp:
                        exec('tmp.%s = IN[\'%s\'][%i]' % (k, k, i))
                    elif IN[k].shape[1] == shp:
                        exec('tmp.%s = IN[\'%s\'][:, %i]' % (k, k, i))
                    elif IN[k].shape[2] == shp:
                        exec('tmp.%s = IN[\'%s\'][:, :, %i]' % (k, k, i))
                    if k == 'gene_idx':
                        tmp.gene_idx = int(tmp.gene_idx[0])
                ev.append(tmp)
            IN.close()
            ev = sp.array(ev, dtype='object')
    else:
        ev = cPickle.load(open(CFG['fname_events'], 'r'))[0]

    cov = quantify.quantify_from_graph(ev, sp.arange(1000), 'exon_skip', CFG, fn_merge=sys.argv[1])

    return cov

def estimate_dispersion_chunk(gene_counts, matrix, sf, CFG, idx, log=False):

    disp_raw = sp.empty((idx.shape[0], 1), dtype='float')
    disp_raw.fill(sp.nan)
    disp_raw_conv = sp.zeros((idx.shape[0], 1), dtype='bool')

    for i in range(idx.shape[0]):

        if log:
            log_progress(i, idx.shape[0])

        disp = 0.1
        resp = gene_counts[i, :].astype('int')

        if sum(resp / sf) < CFG['min_count'] or sp.mean(resp == 0) > 0.6:
            continue

        for j in range(10):
            modNB  = sm.GLM(resp, matrix, family=sm.families.NegativeBinomial(alpha=disp), offset=sp.log(sf))
            result = modNB.fit()

            last_disp = disp
            yhat = result.mu
            sign = -1.0
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = minimize_scalar(likelihood.adj_loglikelihood_scalar, args=(matrix, resp, yhat, sign), method='Bounded', bounds=(0, 10.0), tol=1e-5)
            disp = res.x

            if abs(sp.log(disp) - sp.log(last_disp)) < 1e-4:
                disp_raw[i] = disp
                disp_raw_conv[i] = True
                break
        else:
            disp_raw[i] = disp
            disp_raw_conv[i] = False
    if log:
        log_progress(idx.shape[0], idx.shape[0])

    return (disp_raw, disp_raw_conv, idx)


def estimate_dispersion(gene_counts, matrix, sf, CFG, event_type):
    
    if CFG['verbose']:
        print 'Estimating raw dispersions'

    if CFG['parallel'] > 1:
        disp_raw = sp.empty((gene_counts.shape[0], 1), dtype='float')
        disp_raw.fill(sp.nan)
        disp_raw_conv = sp.zeros((gene_counts.shape[0], 1), dtype='bool')

        pool = mp.Pool(processes=CFG['parallel'], initializer=lambda: sig.signal(sig.SIGINT, sig.SIG_IGN))
        binsize = 30
        idx_chunks = [sp.arange(x, min(x + binsize, gene_counts.shape[0])) for x in range(0, gene_counts.shape[0], binsize)]

        try:
            result = [pool.apply_async(estimate_dispersion_chunk, args=(gene_counts[idx, :], matrix, sf, CFG, idx,)) for idx in idx_chunks]
            res_cnt = 0
            while result:
                tmp = result.pop(0).get()
                for i, j in enumerate(tmp[2]):
                    if CFG['verbose']:
                        log_progress(res_cnt, gene_counts.shape[0])
                        res_cnt += 1
                    disp_raw[j] = tmp[0][i]
                    disp_raw_conv[j] = tmp[1][i]
            if CFG['verbose']:
                log_progress(gene_counts.shape[0], gene_counts.shape[0])
                print ''
            pool.terminate()
            pool.join()
        except KeyboardInterrupt:
            print >> sys.stderr, 'Keyboard Interrupt - exiting'
            pool.terminate()
            pool.join()
            sys.exit(1)
    else:        
        (disp_raw, disp_raw_conv, _) = estimate_dispersion_chunk(gene_counts, matrix, sf, CFG, sp.arange(gene_counts.shape[0]), log=CFG['verbose'])

    if CFG['diagnose_plots']:
        plot.mean_variance_plot(counts=gene_counts,
                                disp=disp_raw,
                                matrix=matrix,
                                figtitle='Raw Dispersion Estimate',
                                filename=os.path.join(CFG['plot_dir'], 'dispersion_raw_%s.png' % event_type),
                                CFG=CFG)

    return (disp_raw, disp_raw_conv)


def fit_dispersion(counts, disp_raw, disp_conv, sf, CFG, dmatrix1, event_type):

    mean_count = sp.mean(counts / sf, axis=1)[:, sp.newaxis]
    index = sp.where(disp_conv)[0]

    lowerBound = sp.percentile(sp.unique(disp_raw[index]), 1)
    upperBound = sp.percentile(sp.unique(disp_raw[index]), 99)

    idx = sp.where((disp_raw > lowerBound) & (disp_raw < upperBound))[0]

    matrix = sp.ones((idx.shape[0], 2), dtype='float')
    matrix[:, 0] /= mean_count[idx].ravel()

    modGamma = sm.GLM(disp_raw[idx], matrix, family=sm.families.Gamma(sm.families.links.identity))
    res = modGamma.fit()
    Lambda = res.params

    disp_fitted = disp_raw.copy()
    ok_idx = sp.where(~sp.isnan(disp_fitted))[0]
    disp_fitted[ok_idx] = Lambda[0] / mean_count[ok_idx] + Lambda[1]

    if sp.sum(disp_fitted > 0) > 0:
        print "Found dispersion fit"

    if CFG['diagnose_plots']:
        plot.mean_variance_plot(counts=counts,
                                disp=disp_fitted,
                                matrix=dmatrix1,
                                figtitle='Fitted Dispersion Estimate',
                                filename=os.path.join(CFG['plot_dir'], 'dispersion_fitted_%s.png' % event_type),
                                CFG=CFG)

    return (disp_fitted, Lambda, idx)


def adj_loglikelihood_shrink_scalar_onedisper(disp, explanatory, response, yhat, dispFitted, varPrior, sign):
    """
    """

    loglik_adj = adj_loglikelihood_scalar(disp, explanatory, response, yhat, 1.0)
    logprior = (sp.log(disp) - sp.log(dispFitted)) ** 2 / (2 * varPrior ** 2)
    loglik_adj_shrk = loglik_adj - logprior

    return loglik_adj_shrk * sign


def adj_loglikelihood_scalar(disp, X, y, mu, sign):

    n = 1 / disp
    p = n / (n + mu)
    loglik = sum(nbinom.logpmf(y, n, p))

    diagVec = mu / (1 + mu * disp)
    diagWM = sp.diag(diagVec)
    xtwx = sp.dot(sp.dot(X.T, diagWM), X)
    coxreid = 0.5 * sp.log(sp.linalg.det(xtwx))

    return (loglik - coxreid) * sign


def adjust_dispersion_chunk(counts, dmatrix1, disp_raw, disp_fitted, varPrior, sf, CFG, idx, log=False):

    disp_adj = sp.empty((counts.shape[0], 1))
    disp_adj.fill(sp.nan)
    disp_adj_conv = sp.zeros_like(disp_adj, dtype='bool')
    error_cnt = 0

    for i in range(idx.shape[0]):

        if log:
            log_progress(i, idx.shape[0])

        if not sp.isnan(disp_raw[i]):

            ### init dispersion and response
            disp = 0.1
            resp = counts[i, :].astype('int')

            ### run for max 10 iterations
            for j in range(10):
                modNB = sm.GLM(resp, dmatrix1, family=sm.families.NegativeBinomial(alpha=disp), offset=sp.log(sf))
                result = modNB.fit()

                dispBef = disp
                yhat = result.mu
                sign = -1.0
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    try:
                        res = minimize_scalar(adj_loglikelihood_shrink_scalar_onedisper, args=(dmatrix1, resp, yhat, disp_fitted[i], varPrior, sign), method='Bounded', bounds=(0, 10.0), tol=1e-5)
                    except TypeError:
                        disp_adj[i] = disp 
                        disp_adj_conv[i] = False
                        error_cnt += 1
                        break
                disp = res.x

                if abs(sp.log(disp) - sp.log(dispBef)) < 1e-4:
                    disp_adj[i] = disp
                    disp_adj_conv[i] = True
                    break
            else:
                disp_adj[i] = disp
                disp_adj_conv[i] = False
    if log:
        log_progress(idx.shape[0], idx.shape[0])
        print ''

    if error_cnt > 0:
        print 'Warning: %i events did not fit due to a TypeError' % error_cnt

    return (disp_adj, disp_adj_conv, idx)


def adjust_dispersion(counts, dmatrix1, disp_raw, disp_fitted, idx, sf, CFG, event_type):

    if CFG['verbose']:
        print 'Estimating adjusted dispersions.'

    varLogDispSamp = polygamma(1, (dmatrix1.shape[0] - dmatrix1.shape[1] ) / 2) ## number of samples - number of coefficients
    varPrior = calculate_varPrior(disp_raw, disp_fitted, idx, varLogDispSamp)

    if CFG['parallel'] > 1:
        disp_adj = sp.empty((counts.shape[0], 1))
        disp_adj.fill(sp.nan)
        disp_adj_conv = sp.zeros_like(disp_adj, dtype='bool')

        pool = mp.Pool(processes=CFG['parallel'], initializer=lambda: sig.signal(sig.SIGINT, sig.SIG_IGN))
        binsize = 30
        idx_chunks = [sp.arange(x, min(x + binsize, counts.shape[0])) for x in range(0, counts.shape[0], binsize)]

        try:
            result = [pool.apply_async(adjust_dispersion_chunk, args=(counts[cidx, :], dmatrix1, disp_raw[cidx], disp_fitted[cidx], varPrior, sf, CFG, cidx,)) for cidx in idx_chunks]
            res_cnt = 0
            while result:
                tmp = result.pop(0).get()
                for i, j in enumerate(tmp[2]):
                    if CFG['verbose']:
                        log_progress(res_cnt, counts.shape[0])
                        res_cnt += 1
                    disp_adj[j] = tmp[0][i]
                    disp_adj_conv[j] = tmp[1][i]
            if CFG['verbose']:
                log_progress(counts.shape[0], counts.shape[0])
                print ''
            pool.terminate()
            pool.join()
        except KeyboardInterrupt:
            print >> sys.stderr, 'Keyboard Interrupt - exiting'
            pool.terminate()
            pool.join()
            sys.exit(1)
    else:        
        (disp_adj, disp_adj_conv, _) = adjust_dispersion_chunk(counts, dmatrix1, disp_raw, disp_fitted, varPrior, sf, CFG, sp.arange(counts.shape[0]), log=CFG['verbose'])

    if CFG['diagnose_plots']:
        plot.mean_variance_plot(counts=counts,
                           disp=disp_adj,
                           matrix=dmatrix1,
                           figtitle='Adjusted Dispersion Estimate',
                           filename=os.path.join(CFG['plot_dir'], 'dispersion_adjusted_%s.png' % event_type),
                           CFG=CFG)

    return (disp_adj, disp_adj_conv)


def test_count_chunk(gene_counts, disp_adj, sf, dmatrix0, dmatrix1, CFG, idx, log=False):

    pval = sp.zeros((gene_counts.shape[0], 1), dtype='float')
    pval.fill(sp.nan)

    for i in xrange(idx.shape[0]):

        if log:
            log_progress(i, idx.shape[0])

        if sp.isnan(disp_adj[i]):
            continue

        response = gene_counts[i, :].astype('int')

        if sp.sum(response[:response.shape[0] / 2] == 0) > CFG['max_0_frac'] * response.shape[0] / 2:
            continue
        modNB0 = sm.GLM(response, dmatrix0, family=sm.families.NegativeBinomial(alpha=disp_adj[i]), offset=sp.log(sf))
        modNB1 = sm.GLM(response, dmatrix1, family=sm.families.NegativeBinomial(alpha=disp_adj[i]), offset=sp.log(sf))
        try:
            result0 = modNB0.fit()
            result1 = modNB1.fit()
        except:
            print >> sys.stderr, '\nWARNING: SVD did not converge - skipping'
            #traceback.print_exc(file=sys.stderr)
            continue

        pval[i] = 1 - chi2.cdf(result0.deviance - result1.deviance, dmatrix1.shape[1] - dmatrix0.shape[1])

    if log:
        log_progress(idx.shape[0], idx.shape[0])
        print ''

    return (pval, idx)


def test_count(gene_counts, disp_adj, sf, dmatrix0, dmatrix1, CFG):

    if CFG['verbose']:
        print 'Running the statistical test.'

    if CFG['parallel'] > 1:
        pval = sp.zeros((gene_counts.shape[0], 1), dtype='float')
        pval.fill(sp.nan)

        pool = mp.Pool(processes=CFG['parallel'], initializer=lambda: sig.signal(sig.SIGINT, sig.SIG_IGN))
        binsize = 30
        idx_chunks = [sp.arange(x, min(x + binsize, gene_counts.shape[0])) for x in range(0, gene_counts.shape[0], binsize)]

        try:
            result = [pool.apply_async(test_count_chunk, args=(gene_counts[cidx, :], disp_adj[cidx], sf, dmatrix0, dmatrix1, CFG, cidx,)) for cidx in idx_chunks]
            res_cnt = 0
            while result:
                tmp = result.pop(0).get()
                for i, j in enumerate(tmp[1]):
                    if CFG['verbose']:
                        log_progress(res_cnt, gene_counts.shape[0])
                        res_cnt += 1
                    pval[j] = tmp[0][i]
            if CFG['verbose']:
                log_progress(gene_counts.shape[0], gene_counts.shape[0])
                print ''
            pool.terminate()
            pool.join()
        except KeyboardInterrupt:
            print >> sys.stderr, 'Keyboard Interrupt - exiting'
            pool.terminate()
            pool.join()
            sys.exit(1)
    else:        
        (pval, _) = test_count_chunk(gene_counts, disp_adj, sf, dmatrix0, dmatrix1, CFG, sp.arange(gene_counts.shape[0]), log=CFG['verbose'])

    if CFG['verbose']:
        print ''

    return pval


def adj_pval(pvals, CFG):
    """
    Perform multiple testing correction.
    """

    pvals_adj = pvals.copy()
    idx = ~sp.isnan(pvals)

    if CFG['multiTest'] == 'BH':
        method = 'fdr_bh'
    elif CFG['multiTest'] == 'Bonferroni':
        method = 'bonferroni'
    elif CFG['multiTest'] == 'Holm':
        method = 'holm'
    elif CFG['multiTest'] == 'Hochberg':
        method = 'simes-hochberg'
    elif CFG['multiTest'] == 'Hommel':
        method = 'hommel'
    elif CFG['multiTest'] == 'BY':
        method = 'fdr_by'
    elif CFG['multiTest'] == 'TSBH':
        method = 'tsbh'
    else:
        sys.stderr.write('ERROR: The methods for multiple test correction can only accept \'Bonferroni\', \'Holm\', \'Hochberg\', \'Hommel\', \'BH\', \'BY\' or \'TSBH\' as its input.\n')
        sys.exit()

    mtc = sms.stats.multicomp.multipletests(pvals[idx], alpha=0.1, method=method, returnsorted=False)

    pvals_adj[idx] = mtc[1]

    return pvals_adj


def calculate_varPrior(disp_raw, disp_fitted, idx, varLogDispSamp):

    logRes = sp.log(disp_raw[idx]) - sp.log(disp_fitted[idx])
    stdLogRes = sp.median(abs(logRes - sp.median(logRes))) * 1.4826

    varLogRes = stdLogRes ** 2
    varPrior = varLogRes - varLogDispSamp

    return max(varPrior, 0.1)


def run_testing(cov, dmatrix0, dmatrix1, sf, CFG, event_type, r_idx=None):

    ### estimate dispersion
    (disp_raw, disp_raw_conv) = estimate_dispersion(cov, dmatrix1, sf, CFG, event_type)

    ### fit dispersion
    (disp_fitted, Lambda, disp_idx) = fit_dispersion(cov, disp_raw, disp_raw_conv, sf, CFG, dmatrix1, event_type)

    ### adjust dispersion estimates
    (disp_adj, disp_adj_conv) = adjust_dispersion(cov, dmatrix1, disp_raw, disp_fitted, disp_idx, sf, CFG, event_type)

    ### do test 
    pvals = test_count(cov, disp_adj, sf, dmatrix0, dmatrix1, CFG)

    ### revert from unique
    if r_idx is not None:
        pvals = pvals[r_idx]

    ### reshape and adjust p-values
    pvals = pvals.reshape((2, pvals.shape[0] / 2)).T
    m_idx = sp.zeros(shape=(pvals.shape[0],), dtype='int')
    for i in range(pvals.shape[0]):
        if sp.all(sp.isnan(pvals[i, :])):
            continue
        elif sp.isnan(pvals[i, 0]):
            m_idx[i] = 1
        elif sp.isnan(pvals[i, 1]):
            m_idx[i] = 0
        else:
            m_idx[i] = sp.argmin(pvals[i, :])
    #m_idx = sp.argmin(pvals, axis=1)
    pvals = 2 * sp.array([pvals[i, m_idx[i]] for i in range(pvals.shape[0])], dtype='float')
    offset = cov.shape[0] / 2
    cov_used = sp.array([cov[i, :] if m_idx[i] == 0 else cov[i + offset, :] for i in range(pvals.shape[0])], dtype=cov.dtype)
    disp_raw_used = sp.array([disp_raw[i] if m_idx[i] == 0 else disp_raw[i + offset] for i in range(pvals.shape[0])], dtype=disp_raw.dtype)
    disp_adj_used = sp.array([disp_adj[i] if m_idx[i] == 0 else disp_adj[i + offset] for i in range(pvals.shape[0])], dtype=disp_adj.dtype)
    pvals[pvals > 1] = 1

    return (pvals, cov_used, disp_raw_used, disp_adj_used)


def main():

    ### get command line options
    options = parse_options(sys.argv)

    ### parse parameters from options object
    CFG = settings.parse_args(options, identity='test')
    CFG['use_exon_counts'] = False

    ### generate output directory
    outdir = os.path.join(options.outdir, 'testing')
    if options.timestamp == 'y':
        outdir = '%s_%s' % (outdir, str(datetime.datetime.now()).replace(' ', '_'))
    if options.labelA != 'condA' and options.labelB != 'condB':
        outdir = '%s_%s_vs_%s' % (outdir, options.labelA, options.labelB)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if CFG['diagnose_plots']:
        CFG['plot_dir'] = os.path.join(outdir, 'plots')
        if not os.path.exists(CFG['plot_dir']):
            os.makedirs(CFG['plot_dir'])

    if CFG['debug']:

        print "Generating simulated dataset"

        npr.seed(23)
        CFG['is_matlab'] = False
        #cov = npr.permutation(20000-20).astype('float').reshape(999, 20)
        #cov = sp.r_[cov, sp.c_[sp.ones((1, 10)) *10, sp.ones((1, 10)) * 500000] + npr.normal(10, 1, 20)]
        #sf = sp.ones((cov.shape[1], ), dtype='float')

        setsize = 50
        ### diff event counts
        cov = sp.zeros((500, 2 * setsize), dtype='int')
        for i in range(10):
            cov[i, :setsize] = nbinom.rvs(30, 0.8, size=setsize)
            cov[i, setsize:] = nbinom.rvs(10, 0.8, size=setsize)
        for i in range(10, cov.shape[0]):
            cov[i, :] = nbinom.rvs(30, 0.8, size=2*setsize)

        ### diff gene expression
        cov2 = sp.zeros((500, 2 * setsize), dtype='int')
        for i in range(20):
            cov2[i, :setsize] = nbinom.rvs(2000, 0.2, size=setsize)
            cov2[i, setsize:] = nbinom.rvs(2000, 0.3, size=setsize)
        for i in range(20, cov2.shape[0]):
            cov2[i, :] = nbinom.rvs(2000, 0.3, size=2*setsize)

        cov = sp.c_[cov, cov2] * 10000

        tidx = sp.arange(setsize)

        sf = npr.uniform(0, 5, 2*setsize)
        sf = sp.r_[sf, sf]

        #dmatrix0 = sp.ones((cov.shape[1], 3), dtype='bool')
        dmatrix1 = sp.zeros((cov.shape[1], 4), dtype='float')
        dmatrix1[:, 0] = 1
        dmatrix1[tidx, 1] = 1
        #dmatrix1[tidx, 2] = 1
        dmatrix1[tidx + (2*setsize), 2] = 1
        dmatrix1[(2*setsize):, 3] = 1
        #dmatrix1[:, 4] = sp.log(sf)
        dmatrix0 = dmatrix1[:, [0, 2, 3]]

        cov = cov * sf
        #sf = sp.ones((cov.shape[1], ), dtype='float')

        pvals = run_testing(cov, dmatrix0, dmatrix1, sf, CFG, 'debug')
        pvals_adj = adj_pval(pvals, CFG) 
        pdb.set_trace()
    else:
        val_tag = ''
        if CFG['validate_splicegraphs']:
            val_tag = '.validated'

        if CFG['is_matlab']:
            CFG['fname_genes'] = os.path.join(CFG['out_dirname'], 'spladder', 'genes_graph_conf%i.%s%s.mat' % (CFG['confidence_level'], CFG['merge_strategy'], val_tag))
            CFG['fname_count_in'] = os.path.join(CFG['out_dirname'], 'spladder', 'genes_graph_conf%i.%s%s.count.mat' % (CFG['confidence_level'], CFG['merge_strategy'], val_tag))
        else:
            CFG['fname_genes'] = os.path.join(CFG['out_dirname'], 'spladder', 'genes_graph_conf%i.%s%s.pickle' % (CFG['confidence_level'], CFG['merge_strategy'], val_tag))
            CFG['fname_count_in'] = os.path.join(CFG['out_dirname'], 'spladder', 'genes_graph_conf%i.%s%s.count.hdf5' % (CFG['confidence_level'], CFG['merge_strategy'], val_tag))

        CFG['fname_exp_hdf5'] = os.path.join(CFG['out_dirname'], 'spladder', 'genes_graph_conf%i.%s%s.gene_exp.hdf5' % (CFG['confidence_level'], CFG['merge_strategy'], val_tag))
        if os.path.exists(CFG['fname_exp_hdf5']):
            if CFG['verbose']:
                print 'Loading expression counts from %s' % CFG['fname_exp_hdf5']
            IN = h5py.File(CFG['fname_exp_hdf5'], 'r')
            gene_counts = IN['raw_count'][:]
            gene_strains = IN['strains'][:]
            gene_ids = IN['genes'][:]
            IN.close()
        else:
            condition_strains = None
            if options.subset_samples == 'y':
                condition_strains = sp.unique(sp.r_[sp.array(CFG['conditionA']), sp.array(CFG['conditionB'])])
                CFG['fname_exp_hdf5'] = os.path.join(CFG['out_dirname'], 'spladder', 'genes_graph_conf%i.%s%s.gene_exp.%i.hdf5' % (CFG['confidence_level'], CFG['merge_strategy'], val_tag, hash(tuple(sp.unique(condition_strains))) * -1))
            if os.path.exists(CFG['fname_exp_hdf5']):
                if CFG['verbose']:
                    print 'Loading expression counts from %s' % CFG['fname_exp_hdf5']
                IN = h5py.File(CFG['fname_exp_hdf5'], 'r')
                gene_counts = IN['raw_count'][:]
                gene_strains = IN['strains'][:]
                gene_ids = IN['genes'][:]
                IN.close()
            else:
                gene_counts, gene_strains, gene_ids = get_gene_expression(CFG, fn_out=CFG['fname_exp_hdf5'], strain_subset=condition_strains)

        gene_strains = sp.array([x.split(':')[1] if ':' in x else x for x in gene_strains])

        ### estimate size factors for library size normalization
        sf_ge = get_size_factors(gene_counts, CFG)

        ### get index of samples for difftest
        idx1 = sp.where(sp.in1d(gene_strains, CFG['conditionA']))[0]
        idx2 = sp.where(sp.in1d(gene_strains, CFG['conditionB']))[0]

        ### for TESTING
        #setsize = 100
        #idx1 = sp.arange(0, setsize / 2)
        #idx2 = sp.arange(setsize / 2, setsize)

        ### subset expression counts to tested samples
        gene_counts = gene_counts[:, sp.r_[idx1, idx2]]
        gene_strains = gene_strains[sp.r_[idx1, idx2]]
        sf_ge = sf_ge[sp.r_[idx1, idx2]]
        #sf = sp.r_[sf, sf]

        ### test each event type individually
        for event_type in CFG['event_types']:

            if CFG['verbose']:
                print 'Testing %s events' % event_type

            CFG['fname_events'] = os.path.join(CFG['out_dirname'], 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, CFG['confidence_level']))

            ### quantify events
            (cov, gene_idx, event_idx, event_ids, event_strains) = quantify.quantify_from_counted_events(CFG['fname_events'], idx1, idx2, event_type, CFG)

            ### estimate size factors
            sf_ev = get_size_factors(sp.vstack(cov), CFG)

            sf = sp.r_[sf_ev, sf_ge]

            assert(sp.all(gene_strains == event_strains))

            ### map gene expression to event order
            curr_gene_counts = gene_counts[gene_idx, :]

            ### filter for min expression
            if event_type == 'intron_retention':
                k_idx = sp.where((sp.mean(cov[0] == 0, axis=1) < CFG['max_0_frac']) | \
                                 (sp.mean(cov[1] == 0, axis=1) < CFG['max_0_frac']))[0]
            else:
                ### filter all events that have too many 0 counts in inclusion or exclusion 
                ### filter all events where inclusion and exclusion in 
                mean1 = sp.mean(sp.c_[(cov[0] / sf_ev)[:, idx1.shape[0]:], (cov[1] / sf_ev)[:, idx1.shape[0]:]], axis=1)
                mean2 = sp.mean(sp.c_[(cov[0] / sf_ev)[:, :idx1.shape[0]], (cov[1] / sf_ev)[:, :idx1.shape[0]]], axis=1)
                k_idx = sp.where(((sp.mean(cov[0] == 0, axis=1) <= CFG['max_0_frac']) | \
                                  (sp.mean(cov[1] == 0, axis=1) <= CFG['max_0_frac'])) & \
                                 (sp.mean(sp.c_[cov[0][:, :idx1.shape[0]], cov[1][:, :idx1.shape[0]]] == 0, axis=1) <= CFG['max_0_frac']) & \
                                 (sp.mean(sp.c_[cov[0][:, idx1.shape[0]:], cov[1][:, idx1.shape[0]:]] == 0, axis=1) <= CFG['max_0_frac']) & \
                                 ((mean1 >= 10) | (mean2 >= 10)))[0]
                del mean1, mean2

            if CFG['verbose']:
                print 'Exclude %i of %i %s events (%.2f percent) from testing due to low coverage' % (cov[0].shape[0] - k_idx.shape[0], cov[0].shape[0], event_type, (1 - float(k_idx.shape[0]) / cov[0].shape[0]) * 100)
            if k_idx.shape[0] == 0:
                print 'All events of type %s were filtered out due to low coverage. Please try re-running with less stringent filter criteria' % event_type
                continue

           # k_idx = sp.where((sp.mean(sp.c_[cov[0], cov[1]], axis=1) > 2))[0]
           # k_idx = sp.where((sp.mean(cov[0], axis=1) > 2) & (sp.mean(cov[1], axis=1) > 2))[0]
            cov[0] = cov[0][k_idx, :]
            cov[1] = cov[1][k_idx, :]
            curr_gene_counts = curr_gene_counts[k_idx, :]
            event_idx = event_idx[k_idx]
            gene_idx = gene_idx[k_idx]
            if not event_ids is None:
                event_ids = [x[k_idx] for x in event_ids]

            cov[0] = sp.around(sp.hstack([cov[0], curr_gene_counts]))
            cov[1] = sp.around(sp.hstack([cov[1], curr_gene_counts]))
            cov = sp.vstack(cov)
            if not event_ids is None:
                event_ids = sp.hstack(event_ids)

            tidx = sp.arange(idx1.shape[0])

        #if CFG['debug']:
        #    for i in range(cov.shape[0]):
        #        fig = plt.figure(figsize=(8, 6), dpi=100)
        #        ax = fig.add_subplot(111)
        #        ax.hist(cov[i, :] * sf, 50, histtype='bar', rwidth=0.8)
        #        #ax.plot(sp.arange(cov.shape[1]), sorted(cov[i, :]), 'bo')
        #        ax.set_title('Count Distribution - Sample %i' % i )
        #        plt.savefig('count_dist.%i.pdf' % i, format='pdf', bbox_inches='tight')
        #        plt.close(fig)

            ### build design matrix for testing
            dmatrix1 = sp.zeros((cov.shape[1], 4), dtype='bool')
            dmatrix1[:, 0] = 1                      # intercept
            dmatrix1[tidx, 1] = 1                   # delta splice
            dmatrix1[tidx, 2] = 1                   # delta gene exp
            dmatrix1[tidx + (idx1.shape[0] + idx2.shape[0]), 2] = 1         # delta gene exp
            dmatrix1[(idx1.shape[0] + idx2.shape[0]):, 3] = 1         # is gene exp
            #dmatrix1[:(idx1.shape[0] + idx2.shape[0]), 4] = 1         # is splice
            dmatrix0 = dmatrix1[:, [0, 2, 3]]
            #dmatrix0 = dmatrix1[:, [0, 2, 3, 4]]

            if CFG['diagnose_plots']:
                plot.count_histogram(counts=cov,
                                     matrix=dmatrix1,
                                     figtitle='Count Distributions',
                                     filename=os.path.join(CFG['plot_dir'], 'count_distribution.%s.png' % event_type),
                                     CFG=CFG)
            ### make event splice forms unique to prevent unnecessary tests
            if not event_ids is None:
                event_ids, u_idx, r_idx = sp.unique(event_ids, return_index=True, return_inverse=True)
                if CFG['verbose']:
                    print 'Consider %i unique event splice forms for testing' % u_idx.shape[0]

            ### run testing
            #pvals = run_testing(cov[u_idx, :], dmatrix0, dmatrix1, sf, CFG, event_type, r_idx)
            (pvals, cov_used, disp_raw_used, disp_adj_used) = run_testing(cov, dmatrix0, dmatrix1, sf, CFG, event_type)
            pvals_adj = adj_pval(pvals, CFG) 

            if CFG['diagnose_plots']:
                plot.qq_plot(pvals=pvals_adj,
                             figtitle='Quantile-Quantile Plot (Adjusted)',
                             filename=os.path.join(CFG['plot_dir'], 'qq_plot_%s.adj.png' % event_type),
                             CFG=CFG)
                plot.qq_plot(pvals=pvals,
                             figtitle='Quantile-Quantile Plot',
                             filename=os.path.join(CFG['plot_dir'], 'qq_plot_%s.png' % event_type),
                             CFG=CFG)

            ### 
            ### OUTPUT
            ###

            ### write test summary (what has been tested, which bam files, etc. ...)
            cPickle.dump((gene_strains,
                          event_strains,
                          dmatrix0,
                          dmatrix1,
                          event_type,
                          options,
                          CFG),
                         open(os.path.join(outdir, 'test_setup_C%i_%s.pickle' % (options.confidence, event_type)), 'w'), -1)

            ### write test results
            s_idx = sp.argsort(pvals)
            header = sp.array(['event_id', 'gene', 'p_val', 'p_val_adj']) 
            event_ids = sp.array(['%s_%i' % (event_type, i + 1) for i in event_idx], dtype='str')

            out_fname = os.path.join(outdir, 'test_results_C%i_%s.tsv' % (options.confidence, event_type))
            if CFG['verbose']:
                print 'Writing test results to %s' % out_fname
            if CFG['is_matlab']:
                data_out = sp.c_[event_ids[s_idx], gene_ids[gene_idx[s_idx], 0], pvals[s_idx].astype('str'), pvals_adj[s_idx].astype('str')]
            else:
                data_out = sp.c_[event_ids[s_idx], gene_ids[gene_idx[s_idx]], pvals[s_idx].astype('str'), pvals_adj[s_idx].astype('str')]
            sp.savetxt(out_fname, sp.r_[header[sp.newaxis, :], data_out], delimiter='\t', fmt='%s')

            ### write extended output
            out_fname = os.path.join(outdir, 'test_results_extended_C%i_%s.tsv' % (options.confidence, event_type))
            if CFG['verbose']:
                print 'Writing extended test results to %s' % out_fname
            header_long = sp.r_[header, ['mean_event_count_A', 'mean_event_count_B', 'log2FC_event_count', 'mean_gene_exp_A', 'mean_gene_exp_B', 'log2FC_gene_exp'], ['event_count:%s' % x for x in event_strains], ['gene_exp:%s' % x for x in event_strains], ['disp_raw', 'disp_adj']]

            ### compute means
            s = event_strains.shape[0]
            m_ev1 = sp.nanmean(cov_used[:, :idx1.shape[0]] / sf_ev[:idx1.shape[0]], axis=1)
            m_ev2 = sp.nanmean(cov_used[:, idx1.shape[0]:s] / sf_ev[idx1.shape[0]:], axis=1)
            fc_ev = sp.log2(m_ev1) - sp.log2(m_ev2)
            m_ge1 = sp.nanmean(cov_used[:, s:s+idx1.shape[0]] / sf_ge[:idx1.shape[0]], axis=1)
            m_ge2 = sp.nanmean(cov_used[:, s+idx1.shape[0]:] / sf_ge[idx1.shape[0]:], axis=1)
            fc_ge = sp.log2(m_ge1) - sp.log2(m_ge2)
            m_all = sp.c_[m_ev1, m_ev2, fc_ev, m_ge1, m_ge2, fc_ge]

            data_out = sp.c_[data_out, m_all[s_idx, :], (cov_used[s_idx, :] / sf).astype('str'), disp_raw_used.astype('str'), disp_adj_used.astype('str')]
            sp.savetxt(out_fname, sp.r_[header_long[sp.newaxis, :], data_out], delimiter='\t', fmt='%s')

            ### write output unique over genes
            out_fname = os.path.join(outdir, 'test_results_C%i_%s.gene_unique.tsv' % (options.confidence, event_type))
            if CFG['verbose']:
                print 'Writing gene unique test results to %s' % out_fname
            ks_idx = []
            taken = set()
            for i in s_idx:
                if CFG['is_matlab']:
                    gid = gene_ids[gene_idx[i], 0]
                else:
                    gid = gene_ids[gene_idx[i]]
                if gid in taken:
                    continue
                ks_idx.append(i)
                taken.add(gid)
            ks_idx = sp.array(ks_idx)

            if CFG['is_matlab']:
                data_out = sp.c_[event_ids[ks_idx], gene_ids[gene_idx[ks_idx], 0], pvals[ks_idx].astype('str'), pvals_adj[ks_idx].astype('str')]
            else:
                data_out = sp.c_[event_ids[ks_idx], gene_ids[gene_idx[ks_idx]], pvals[ks_idx].astype('str'), pvals_adj[ks_idx].astype('str')]
            data_out = sp.r_[header[sp.newaxis, :], data_out]
            sp.savetxt(out_fname, data_out, delimiter='\t', fmt='%s')

if __name__ == "__main__":
    main()

