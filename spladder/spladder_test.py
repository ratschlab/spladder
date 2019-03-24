import scipy as sp
import scipy.io as scio
import statsmodels.api as sm
import statsmodels.sandbox as sms
import h5py
import sys
import os
import pdb
import pickle
import warnings
import time
import datetime
from scipy.optimize import minimize_scalar
from scipy.special import polygamma
from scipy.stats import chi2,nbinom,scoreatpercentile
import numpy.random as npr

from .alt_splice import quantify
from .testing import likelihood
from . import settings
from .viz import diagnose as plot

import multiprocessing as mp
import signal as sig

from .helpers import log_progress, decodeUTF8, codeUTF8

TIME0 = time.time()

def get_gene_expression(options, fn_out=None, strain_subset=None):

    if options.verbose:
        sys.stdout.write('Quantifying gene expression ...\n')

    ### load gene information
    genes = pickle.load(open(options.fname_genes, 'rb'))[0]
    numgenes = genes.shape[0]

    ### open hdf5 file containing graph count information
    IN = h5py.File(options.fname_count_in, 'r')
    strains = IN['strains'][:].astype('str')
    ### sort by strain ID
    strain_idx_all = sp.argsort(strains)
    if strain_subset is None:
        strain_idx = strain_idx_all.copy()
    else:
        strain_idx = strain_idx_all[sp.in1d(strains[strain_idx_all], strain_subset)]
    gene_counts = sp.zeros((numgenes, strain_idx.shape[0]), dtype='float')
    gene_names = sp.array([x.name for x in genes], dtype='str')

    seg_lens = IN['seg_len'][:]
    gene_ids_segs = IN['gene_ids_segs'][:].astype('int')

    ### no longer assume that the gene_ids_segs are sorted by gene ID
    s_idx = sp.argsort(gene_ids_segs[:, 0], kind='mergesort')
    _, u_idx = sp.unique(gene_ids_segs[s_idx, 0], return_index=True)
    s_idx = s_idx[u_idx]

    ### iterate over genes
    for gidx, iidx in enumerate(s_idx):

        if options.verbose:
            log_progress(gidx, numgenes, 100)

        ### get idx of non alternative segments
        non_alt_idx = genes[gidx].get_non_alt_seg_ids()
        seg_idx = sp.arange(iidx, iidx + genes[gidx].segmentgraph.seg_edges.shape[0])

        gene_idx = gene_ids_segs[seg_idx, 0]
        if len(gene_idx.shape) > 0:
            gene_idx = gene_idx[0]

        assert(decodeUTF8(IN['gene_names'][:][gene_idx]) == genes[gidx].name)
        assert(genes[gidx].name == gene_names[gidx])

        if options.non_alt_norm:
            seg_idx = seg_idx[non_alt_idx]

        ### compute gene expression as the read count over all non alternative segments
        #gene_counts[gidx, :] = sp.dot(IN['segments'][seg_idx, :].T, IN['seg_len'][:][seg_idx]) / sp.sum(IN['seg_len'][:][seg_idx])
        if seg_idx.shape[0] > 1:
            gene_counts[gidx, :] = sp.squeeze(sp.dot(IN['segments'][seg_idx, :][:, strain_idx].T, seg_lens[seg_idx])) / options.readlen
        else:
            gene_counts[gidx, :] = IN['segments'][seg_idx[0], :][strain_idx] * seg_lens[seg_idx] / options.readlen
        #seg_offset += genes[gidx].segmentgraph.seg_edges.shape[0]

    IN.close()

    if options.verbose:
        sys.stdout.write('\n... done.\n')


    ### write results to hdf5
    if fn_out is not None:
        OUT = h5py.File(fn_out, 'w')
        OUT.create_dataset(name='all_strains', data=codeUTF8(strains[strain_idx_all]))
        OUT.create_dataset(name='strains', data=codeUTF8(strains[strain_idx]))
        OUT.create_dataset(name='genes', data=codeUTF8(gene_names))
        OUT.create_dataset(name='raw_count', data=gene_counts, compression="gzip")
        OUT.close()

    return (gene_counts, strains[strain_idx_all], strains[strain_idx], gene_names)


def get_size_factors(gene_counts, options):

    if options.verbose:
        print('Estimating size factors')

    ### take geometric mean of counts
    gmean = sp.exp(sp.mean(sp.log(gene_counts + 1), axis=1))

    size_factors = []
    for i in range(gene_counts.shape[1]):
        idx = gene_counts[:, i] > 0
        size_factors.append(sp.median(gene_counts[idx, i] / gmean[idx]))

    size_factors = sp.array(size_factors, dtype='float')

    return size_factors


def re_quantify_events(options):
    """This is more a legacy function for testing that requantifies events on a given graph"""

    ev = pickle.load(open(options.fname_events, 'rb'))[0]
    cov = quantify.quantify_from_graph(ev, sp.arange(1000), 'exon_skip', options, fn_merge=sys.argv[1])

    return cov

def estimate_dispersion_chunk(gene_counts, matrix, sf, options, test_idx, idx, log=False):

    disp_raw = sp.empty((idx.shape[0], 1), dtype='float')
    disp_raw.fill(sp.nan)
    disp_raw_conv = sp.zeros((idx.shape[0], 1), dtype='bool')

    npr.seed(23)
    for i in range(idx.shape[0]):

        if log:
            log_progress(i, idx.shape[0])

        disp = 0.1
        resp = gene_counts[i, :].astype('int')

        if sum(resp / sf) < options.min_count or sp.mean(resp == 0) > 0.6 or not test_idx[i]:
            continue

        for j in range(10):
            modNB  = sm.GLM(resp, matrix, family=sm.families.NegativeBinomial(alpha=disp), offset=sp.log(sf))
            result = modNB.fit()

            sp.set_printoptions(12)

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


def estimate_dispersion(gene_counts, matrix, sf, options, test_idx, event_type):

    if options.verbose:
        print('Estimating raw dispersions')

    if options.parallel > 1:
        disp_raw = sp.empty((gene_counts.shape[0], 1), dtype='float')
        disp_raw.fill(sp.nan)
        disp_raw_conv = sp.zeros((gene_counts.shape[0], 1), dtype='bool')

        pool = mp.Pool(processes=options.parallel, initializer=lambda: sig.signal(sig.SIGINT, sig.SIG_IGN))
        binsize = 30
        idx_chunks = [sp.arange(x, min(x + binsize, gene_counts.shape[0])) for x in range(0, gene_counts.shape[0], binsize)]

        try:
            result = [pool.apply_async(estimate_dispersion_chunk, args=(gene_counts[idx, :], matrix, sf, options, test_idx[idx], idx,)) for idx in idx_chunks]
            res_cnt = 0
            while result:
                tmp = result.pop(0).get()
                for i, j in enumerate(tmp[2]):
                    if options.verbose:
                        log_progress(res_cnt, gene_counts.shape[0])
                        res_cnt += 1
                    disp_raw[j] = tmp[0][i]
                    disp_raw_conv[j] = tmp[1][i]
            if options.verbose:
                log_progress(gene_counts.shape[0], gene_counts.shape[0])
                print('')
            pool.terminate()
            pool.join()
        except KeyboardInterrupt:
            print('Keyboard Interrupt - exiting', file=sys.stderr)
            pool.terminate()
            pool.join()
            sys.exit(1)
    else:
        (disp_raw, disp_raw_conv, _) = estimate_dispersion_chunk(gene_counts, matrix, sf, options, test_idx, sp.arange(gene_counts.shape[0]), log=options.verbose)

    if sp.sum(disp_raw_conv) == 0:
        print('\nERROR: None of the dispersion estimates converged. Exiting.', file=sys.stderr)
        sys.exit(1)

    if options.diagnose_plots:
        plot.mean_variance_plot(counts=gene_counts,
                                disp=disp_raw,
                                matrix=matrix,
                                figtitle='Raw Dispersion Estimate',
                                filename=os.path.join(options.plot_dir, 'dispersion_raw_%s.%s' % (event_type, options.plot_format)),
                                options=options)

    return (disp_raw, disp_raw_conv)


def fit_dispersion(counts, disp_raw, disp_conv, sf, options, dmatrix1, event_type):

    mean_count = sp.mean(counts / sf, axis=1)[:, sp.newaxis]
    index = sp.where(disp_conv)[0]

    lowerBound = sp.percentile(sp.unique(disp_raw[index]), 1)
    upperBound = sp.percentile(sp.unique(disp_raw[index]), 99)

    idx = sp.where((disp_raw > lowerBound) & (disp_raw < upperBound))[0]

    matrix = sp.ones((idx.shape[0], 2), dtype='float')
    matrix[:, 0] /= mean_count[idx].ravel()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        modGamma = sm.GLM(disp_raw[idx], matrix, family=sm.families.Gamma(sm.families.links.identity()))
    res = modGamma.fit()
    Lambda = res.params

    disp_fitted = disp_raw.copy()
    ok_idx = sp.where(~sp.isnan(disp_fitted))[0]
    disp_fitted[ok_idx] = Lambda[0] / mean_count[ok_idx] + Lambda[1]

    if sp.sum(disp_fitted > 0) > 0:
        print("\nFound dispersion fit")

    if options.diagnose_plots:
        plot.mean_variance_plot(counts=counts,
                                disp=disp_fitted,
                                matrix=dmatrix1,
                                figtitle='Fitted Dispersion Estimate',
                                filename=os.path.join(options.plot_dir, 'dispersion_fitted_%s.%s' % (event_type, options.plot_format)),
                                options=options)

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


def adjust_dispersion_chunk(counts, dmatrix1, disp_raw, disp_fitted, varPrior, sf, options, idx, log=False):

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
        print('')

    #if error_cnt > 0:
    #    print 'Warning: %i events did not fit due to a TypeError' % error_cnt

    return (disp_adj, disp_adj_conv, idx)


def adjust_dispersion(counts, dmatrix1, disp_raw, disp_fitted, idx, sf, options, event_type):

    if options.verbose:
        print('Estimating adjusted dispersions.')

    varLogDispSamp = polygamma(1, (dmatrix1.shape[0] - dmatrix1.shape[1] ) / 2) ## number of samples - number of coefficients
    varPrior = calculate_varPrior(disp_raw, disp_fitted, idx, varLogDispSamp)

    if options.parallel > 1:
        disp_adj = sp.empty((counts.shape[0], 1))
        disp_adj.fill(sp.nan)
        disp_adj_conv = sp.zeros_like(disp_adj, dtype='bool')

        pool = mp.Pool(processes=options.parallel, initializer=lambda: sig.signal(sig.SIGINT, sig.SIG_IGN))
        binsize = 30
        idx_chunks = [sp.arange(x, min(x + binsize, counts.shape[0])) for x in range(0, counts.shape[0], binsize)]

        try:
            result = [pool.apply_async(adjust_dispersion_chunk, args=(counts[cidx, :], dmatrix1, disp_raw[cidx], disp_fitted[cidx], varPrior, sf, options, cidx,)) for cidx in idx_chunks]
            res_cnt = 0
            while result:
                tmp = result.pop(0).get()
                for i, j in enumerate(tmp[2]):
                    if options.verbose:
                        log_progress(res_cnt, counts.shape[0])
                        res_cnt += 1
                    disp_adj[j] = tmp[0][i]
                    disp_adj_conv[j] = tmp[1][i]
            if options.verbose:
                log_progress(counts.shape[0], counts.shape[0])
                print('')
            pool.terminate()
            pool.join()
        except KeyboardInterrupt:
            print('Keyboard Interrupt - exiting', file=sys.stderr)
            pool.terminate()
            pool.join()
            sys.exit(1)
    else:
        (disp_adj, disp_adj_conv, _) = adjust_dispersion_chunk(counts, dmatrix1, disp_raw, disp_fitted, varPrior, sf, options, sp.arange(counts.shape[0]), log=options.verbose)

    if options.diagnose_plots:
        plot.mean_variance_plot(counts=counts,
                           disp=disp_adj,
                           matrix=dmatrix1,
                           figtitle='Adjusted Dispersion Estimate',
                           filename=os.path.join(options.plot_dir, 'dispersion_adjusted_%s.%s' % (event_type, options.plot_format)),
                           options=options)

    return (disp_adj, disp_adj_conv)


def test_count_chunk(gene_counts, disp_adj, sf, dmatrix0, dmatrix1, options, test_idx, idx, log=False):

    pval = sp.zeros((gene_counts.shape[0], 1), dtype='float')
    pval.fill(sp.nan)
    npr.seed(23)

    for i in range(idx.shape[0]):

        if log:
            log_progress(i, idx.shape[0])

        if sp.isnan(disp_adj[i]) or not test_idx[i]:
            continue

        response = gene_counts[i, :].astype('int')

        if sp.sum(response[:int(response.shape[0] / 2)] == 0) > options.max_0_frac * response.shape[0] / 2:
            continue
        modNB0 = sm.GLM(response, dmatrix0, family=sm.families.NegativeBinomial(alpha=disp_adj[i]), offset=sp.log(sf))
        modNB1 = sm.GLM(response, dmatrix1, family=sm.families.NegativeBinomial(alpha=disp_adj[i]), offset=sp.log(sf))
        try:
            result0 = modNB0.fit()
            result1 = modNB1.fit()
        except:
            print('\nWARNING: SVD did not converge - skipping', file=sys.stderr)
            #traceback.print_exc(file=sys.stderr)
            continue

        pval[i] = 1 - chi2.cdf(result0.deviance - result1.deviance, dmatrix1.shape[1] - dmatrix0.shape[1])

    if log:
        log_progress(idx.shape[0], idx.shape[0])
        print('')

    return (pval, idx)


def test_count(gene_counts, disp_adj, sf, dmatrix0, dmatrix1, options, test_idx):

    if options.verbose:
        print('Running the statistical test.')

    if options.parallel > 1:
        pval = sp.zeros((gene_counts.shape[0], 1), dtype='float')
        pval.fill(sp.nan)

        pool = mp.Pool(processes=options.parallel, initializer=lambda: sig.signal(sig.SIGINT, sig.SIG_IGN))
        binsize = 30
        idx_chunks = [sp.arange(x, min(x + binsize, gene_counts.shape[0])) for x in range(0, gene_counts.shape[0], binsize)]

        try:
            result = [pool.apply_async(test_count_chunk, args=(gene_counts[cidx, :], disp_adj[cidx], sf, dmatrix0, dmatrix1, options, test_idx[cidx], cidx)) for cidx in idx_chunks]
            res_cnt = 0
            while result:
                tmp = result.pop(0).get()
                for i, j in enumerate(tmp[1]):
                    if options.verbose:
                        log_progress(res_cnt, gene_counts.shape[0])
                        res_cnt += 1
                    pval[j] = tmp[0][i]
            if options.verbose:
                log_progress(gene_counts.shape[0], gene_counts.shape[0])
                print('')
            pool.terminate()
            pool.join()
        except KeyboardInterrupt:
            print('Keyboard Interrupt - exiting', file=sys.stderr)
            pool.terminate()
            pool.join()
            sys.exit(1)
    else:
        (pval, _) = test_count_chunk(gene_counts, disp_adj, sf, dmatrix0, dmatrix1, options, test_idx, sp.arange(gene_counts.shape[0]), log=options.verbose)

    if options.verbose:
        print('')

    return pval


def adj_pval(pvals, options):
    """
    Perform multiple testing correction.
    """

    pvals_adj = pvals.copy()
    idx = ~sp.isnan(pvals)

    if options.correction == 'BH':
        method = 'fdr_bh'
    elif options.correction == 'Bonferroni':
        method = 'bonferroni'
    elif options.correction == 'Holm':
        method = 'holm'
    elif options.correction == 'Hochberg':
        method = 'simes-hochberg'
    elif options.correction == 'Hommel':
        method = 'hommel'
    elif options.correction == 'BY':
        method = 'fdr_by'
    elif options.correction == 'TSBH':
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


def run_testing(cov, dmatrix0, dmatrix1, sf, options, event_type, test_idx, r_idx=None):

    ### estimate dispersion
    (disp_raw, disp_raw_conv) = estimate_dispersion(cov, dmatrix1, sf, options, test_idx, event_type)

    ### fit dispersion
    (disp_fitted, Lambda, disp_idx) = fit_dispersion(cov, disp_raw, (disp_raw_conv[:, 0] & test_idx)[:, sp.newaxis], sf, options, dmatrix1, event_type)

    ### adjust dispersion estimates
    (disp_adj, disp_adj_conv) = adjust_dispersion(cov, dmatrix1, disp_raw, disp_fitted, disp_idx, sf, options, event_type)

    ### do test
    pvals = test_count(cov, disp_adj, sf, dmatrix0, dmatrix1, options, test_idx)

    ### revert from unique
    if r_idx is not None:
        pvals = pvals[r_idx]

    ### reshape and adjust p-values
    pvals = pvals.reshape((2, int(pvals.shape[0] / 2))).T
    m_idx = sp.zeros(shape=(pvals.shape[0],), dtype='int')
    #for i in range(pvals.shape[0]):
    #    if sp.all(sp.isnan(pvals[i, :])):
    #        continue
    #    elif sp.isnan(pvals[i, 0]):
    #        m_idx[i] = 1
    #    elif sp.isnan(pvals[i, 1]):
    #        m_idx[i] = 0
    #    else:
    #        m_idx[i] = sp.argmin(pvals[i, :])
    #pvals = 2 * sp.array([pvals[i, m_idx[i]] for i in range(pvals.shape[0])], dtype='float')
    for i in range(pvals.shape[0]):
        if sp.all(sp.isnan(pvals[i, :])):
            continue
        elif sp.isnan(pvals[i, 0]):
            #pvals[i, 1] = sp.nan
            m_idx[i] = 1
        elif sp.isnan(pvals[i, 1]):
            #pvals[i, 0] = sp.nan
            m_idx[i] = 0
        else:
            m_idx[i] = sp.argmax(pvals[i, :])
            #m_idx[i] = sp.argmin(pvals[i, :])
            #pvals[i, m_idx[i]] = min(1, 2*pvals[i, m_idx[i]])

    pvals = sp.array([pvals[i, m_idx[i]] for i in range(pvals.shape[0])], dtype='float')

    offset = int(cov.shape[0] / 2)
    cov_used = sp.array([cov[i, :] if m_idx[i] == 0 else cov[i + offset, :] for i in range(pvals.shape[0])], dtype=cov.dtype)
    disp_raw_used = sp.array([disp_raw[i] if m_idx[i] == 0 else disp_raw[i + offset] for i in range(pvals.shape[0])], dtype=disp_raw.dtype)
    disp_adj_used = sp.array([disp_adj[i] if m_idx[i] == 0 else disp_adj[i + offset] for i in range(pvals.shape[0])], dtype=disp_adj.dtype)
    pvals[pvals > 1] = 1

    return (pvals, cov_used, disp_raw_used, disp_adj_used)


def spladder_test(options):

    ### parse parameters from options object
    options = settings.parse_args(options, identity='test')
    options.use_exon_counts = False

    non_alt_tag = ''
    if options.non_alt_norm:
        non_alt_tag = '.non_alt'

    ### generate output directory
    outdir = os.path.join(options.outdir, 'testing%s' % non_alt_tag)
    if options.timestamp == 'y':
        outdir = '%s_%s' % (outdir, str(datetime.datetime.now()).replace(' ', '_'))
    if options.labelA != 'condA' and options.labelB != 'condB':
        outdir = '%s_%s_vs_%s' % (outdir, options.labelA, options.labelB)
    if options.out_tag != '-':
        outdir += '_%s' % options.out_tag
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if options.diagnose_plots:
        options.plot_dir = os.path.join(outdir, 'plots')
        if not os.path.exists(options.plot_dir):
            os.makedirs(options.plot_dir)

    val_tag = ''
    if options.validate_sg:
        val_tag = '.validated'

    options.fname_genes = os.path.join(options.outdir, 'spladder', 'genes_graph_conf%i.%s%s.pickle' % (options.confidence, options.merge, val_tag))
    options.fname_count_in = os.path.join(options.outdir, 'spladder', 'genes_graph_conf%i.%s%s.count.hdf5' % (options.confidence, options.merge, val_tag))
    options.fname_exp_hdf5 = os.path.join(options.outdir, 'spladder', 'genes_graph_conf%i.%s%s.gene_exp%s.hdf5' % (options.confidence, options.merge, val_tag, non_alt_tag))

    condition_strains = None
    if options.subset_samples:
        condition_strains = sp.unique(sp.r_[sp.array(Coptions.conditionA), sp.array(options.conditionB)])
        options.fname_exp_hdf5 = os.path.join(options.outdir, 'spladder', 'genes_graph_conf%i.%s%s.gene_exp%s.%i.hdf5' % (options.confidence_level, options.merge, val_tag, non_alt_tag, hash(tuple(sp.unique(condition_strains))) * -1))
    if os.path.exists(options.fname_exp_hdf5):
        if options.verbose:
            print('Loading expression counts from %s' % options.fname_exp_hdf5)
        IN = h5py.File(options.fname_exp_hdf5, 'r')
        gene_counts = IN['raw_count'][:]
        gene_strains = decodeUTF8(IN['strains'][:])
        gene_strains_all = decodeUTF8(IN['all_strains'][:])
        gene_ids = decodeUTF8(IN['genes'][:])
        IN.close()
    else:
        gene_counts, gene_strains_all, gene_strains, gene_ids = get_gene_expression(options, fn_out=options.fname_exp_hdf5, strain_subset=condition_strains)

    gene_strains = sp.array([x.split(':')[1] if ':' in x else x for x in gene_strains])
    gene_strains_all = sp.array([x.split(':')[1] if ':' in x else x for x in gene_strains_all])

    ### get index of samples for difftest
    idx1 = sp.where(sp.in1d(gene_strains, options.conditionA))[0]
    idx2 = sp.where(sp.in1d(gene_strains, options.conditionB))[0]
    idx1_all = sp.where(sp.in1d(gene_strains_all, options.conditionA))[0]
    idx2_all = sp.where(sp.in1d(gene_strains_all, options.conditionB))[0]

    ### subset expression counts to tested samples
    gene_counts = gene_counts[:, sp.r_[idx1, idx2]]
    gene_strains = gene_strains[sp.r_[idx1, idx2]]

    ### estimate size factors for library size normalization
    sf_ge = get_size_factors(gene_counts, options)

    ### handle outliers (mask with capped value)
    if options.cap_exp_outliers:
        outlier_cnt = 0
        for gidx in range(gene_counts.shape[0]):
            log_counts = sp.log2(gene_counts[gidx, :] / sf_ge + 1)
            p25 = scoreatpercentile(log_counts, 25)
            p75 = scoreatpercentile(log_counts, 75)
            iqr = (p75 - p25)
            if iqr > 0:
                o_idx = sp.where(log_counts > p75+(1.5*iqr))[0]
                if o_idx.shape[0] > 0:
                    cap = 2**(p75+(1.5*iqr))-1
                    gene_counts[gidx, o_idx] = (cap * sf_ge[o_idx])
                    outlier_cnt += o_idx.shape[0]
        if outlier_cnt > 0:
            total_cnt = gene_counts.shape[0] * gene_counts.shape[1]
            sys.stdout.write('\nCapped %i/%i outlier expression counts (%.2f percent)\n' % (outlier_cnt, total_cnt, float(outlier_cnt) / total_cnt * 100))

    ### test each event type individually
    for event_type in options.event_types:

        if options.verbose:
            print('Testing %s events' % event_type)

        options.fname_events = os.path.join(options.outdir, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, options.confidence))

        ### quantify events
        (cov, gene_idx, event_idx, event_ids, event_strains) = quantify.quantify_from_counted_events(options.fname_events, idx1_all, idx2_all, event_type, options, gen_event_ids=False, low_mem=options.low_memory)

        if options.cap_outliers:
            log_counts = sp.log2(cov[0] + 1)
            p25 = scoreatpercentile(log_counts, 25, axis=1)
            p75 = scoreatpercentile(log_counts, 75, axis=1)
            iqr = (p75 - p25)
            cap = 2**(p75 + (3*iqr)) - 1
            for c in sp.where(iqr > 0)[0]:
                cov[0][c, log_counts[c, :] > (p75[c] + 3*iqr[c])] = cap[c]
            log_counts = sp.log2(cov[1] + 1)
            p25 = scoreatpercentile(log_counts, 25, axis=1)
            p75 = scoreatpercentile(log_counts, 75, axis=1)
            iqr = (p75 - p25)
            cap = 2**(p75 + (3*iqr)) - 1
            for c in sp.where(iqr > 0)[0]:
                cov[1][c, log_counts[c, :] > (p75[c] + 3*iqr[c])] = cap[c]

        ### estimate size factors
        sf_ev = get_size_factors(sp.vstack(cov), options)

        sf = sp.r_[sf_ev, sf_ge]

        assert(sp.all(gene_strains == event_strains))

        ### map gene expression to event order
        curr_gene_counts = gene_counts[gene_idx, :]

        ### filter for min expression
        k_idx1 = ((sp.mean(cov[0][:, idx1.shape[0]:] <= 1, axis=1) <= options.max_0_frac) | \
                  (sp.mean(cov[0][:, :idx1.shape[0]] <= 1, axis=1) <= options.max_0_frac))
        k_idx2 = ((sp.mean(cov[1][:, idx1.shape[0]:] <= 1, axis=1) <= options.max_0_frac) | \
                  (sp.mean(cov[1][:, :idx1.shape[0]] <= 1, axis=1) <= options.max_0_frac))

        k_idx = sp.where(k_idx1 | k_idx2)[0]

        if options.verbose:
            print('Exclude %i of %i %s events (%.2f percent) from testing due to low coverage' % (cov[0].shape[0] - k_idx.shape[0], cov[0].shape[0], event_type, (1 - float(k_idx.shape[0]) / cov[0].shape[0]) * 100))
        if k_idx.shape[0] == 0:
            print('All events of type %s were filtered out due to low coverage. Please try re-running with less stringent filter criteria' % event_type)
            continue

        cov[0] = cov[0][k_idx, :]
        cov[1] = cov[1][k_idx, :]
        curr_gene_counts = curr_gene_counts[k_idx, :]
        event_idx = event_idx[k_idx]
        gene_idx = gene_idx[k_idx]
        if not event_ids is None:
            event_ids = [x[k_idx] for x in event_ids]
        k_idx1 = k_idx1[k_idx]
        k_idx2 = k_idx2[k_idx]

        cov[0] = sp.around(sp.hstack([cov[0], curr_gene_counts]))
        cov[1] = sp.around(sp.hstack([cov[1], curr_gene_counts]))
        cov = sp.vstack(cov)
        if not event_ids is None:
            event_ids = sp.hstack(event_ids)
        test_idx = sp.r_[k_idx1, k_idx2]

        tidx = sp.arange(idx1.shape[0])

    #if options.debug:
    #    for i in range(cov.shape[0]):
    #        fig = plt.figure(figsize=(8, 6), dpi=100)
    #        ax = fig.add_subplot(111)
    #        ax.hist(cov[i, :] * sf, 50, histtype='bar', rwidth=0.8)
    #        #ax.plot(sp.arange(cov.shape[1]), sorted(cov[i, :]), 'bo')
    #        ax.set_title('Count Distribution - Sample %i' % i )
    #        plt.savefig('count_dist.%i.pdf' % i, format='pdf', bbox_inches='tight')
    #        plt.close(fig)

        ### build design matrix for testing
        dmatrix1 = sp.zeros((cov.shape[1], 4), dtype='int')
        dmatrix1[:, 0] = 1                      # intercept
        dmatrix1[tidx, 1] = 1                   # delta splice
        dmatrix1[tidx, 2] = 1                   # delta gene exp
        dmatrix1[tidx + (idx1.shape[0] + idx2.shape[0]), 2] = 1         # delta gene exp
        dmatrix1[(idx1.shape[0] + idx2.shape[0]):, 3] = 1         # is gene exp
        #dmatrix1[:(idx1.shape[0] + idx2.shape[0]), 4] = 1         # is splice
        dmatrix0 = dmatrix1[:, [0, 2, 3]]
        #dmatrix0 = dmatrix1[:, [0, 2, 3, 4]]

        if options.diagnose_plots:
            plot.count_histogram(counts=cov[test_idx, :],
                                 matrix=dmatrix1,
                                 figtitle='Count Distributions',
                                 filename=os.path.join(options.plot_dir, 'count_distribution.%s.%s' % (event_type, options.plot_format)),
                                 options=options)
        ### make event splice forms unique to prevent unnecessary tests
        if not event_ids is None:
            event_ids, u_idx, r_idx = sp.unique(event_ids, return_index=True, return_inverse=True)
            test_idx = test_idx[u_idx]
            if options.verbose:
                print('Consider %i unique event splice forms for testing' % u_idx.shape[0])

        ### run testing
        (pvals, cov_used, disp_raw_used, disp_adj_used) = run_testing(cov, dmatrix0, dmatrix1, sf, options, event_type, test_idx)
        pvals_adj = adj_pval(pvals, options)

        ### compute means and fold changes
        s = event_strains.shape[0]
        m_ev1 = sp.nanmean(cov_used[:, :idx1.shape[0]] / sf_ev[:idx1.shape[0]], axis=1)
        m_ev2 = sp.nanmean(cov_used[:, idx1.shape[0]:s] / sf_ev[idx1.shape[0]:], axis=1)
        fc_ev = sp.log2(m_ev1) - sp.log2(m_ev2)
        m_ge1 = sp.nanmean(cov_used[:, s:s+idx1.shape[0]] / sf_ge[:idx1.shape[0]], axis=1)
        m_ge2 = sp.nanmean(cov_used[:, s+idx1.shape[0]:] / sf_ge[idx1.shape[0]:], axis=1)
        fc_ge = sp.log2(m_ge1) - sp.log2(m_ge2)
        m_all = sp.c_[m_ev1, m_ev2, fc_ev, m_ge1, m_ge2, fc_ge]

        if options.diagnose_plots:
            plot.qq_plot(pvals=pvals_adj,
                         figtitle='Quantile-Quantile Plot (Adjusted)',
                         filename=os.path.join(options.plot_dir, 'qq_plot_%s.adj.%s' % (event_type, options.plot_format)),
                         options=options)
            plot.qq_plot(pvals=pvals,
                         figtitle='Quantile-Quantile Plot',
                         filename=os.path.join(options.plot_dir, 'qq_plot_%s.%s' % (event_type, options.plot_format)),
                         options=options)
            plot.ma_plot(pvals=pvals_adj,
                         counts=sp.nanmean(sp.c_[m_ev1, m_ev2], axis=1),
                         fc=(fc_ev - fc_ge),
                         figtitle='MA Plot',
                         filename=os.path.join(options.plot_dir, 'ma_plot_%s.%s' % (event_type, options.plot_format)),
                         options=options)


        ###
        ### OUTPUT
        ###

        ### write test summary (what has been tested, which bam files, etc. ...)
        pickle.dump((gene_strains,
                      event_strains,
                      dmatrix0,
                      dmatrix1,
                      event_type),
                     open(os.path.join(outdir, 'test_setup_C%i_%s.pickle' % (options.confidence, event_type)), 'wb'), -1)

        ### write test results
        s_idx = sp.argsort(pvals)
        header = sp.array(['event_id', 'gene', 'p_val', 'p_val_adj', 'mean_event_count_A', 'mean_event_count_B', 'log2FC_event_count', 'mean_gene_exp_A', 'mean_gene_exp_B', 'log2FC_gene_exp'])
        event_ids = sp.array(['%s_%i' % (event_type, i + 1) for i in event_idx], dtype='str')

        out_fname = os.path.join(outdir, 'test_results_C%i_%s.tsv' % (options.confidence, event_type))
        if options.verbose:
            print('Writing test results to %s' % out_fname)
        data_out = sp.c_[event_ids[s_idx], gene_ids[gene_idx[s_idx]], pvals[s_idx].astype('str'), pvals_adj[s_idx].astype('str'), m_all[s_idx, :]]
        sp.savetxt(out_fname, sp.r_[header[sp.newaxis, :], data_out], delimiter='\t', fmt='%s')

        ### write extended output
        out_fname = os.path.join(outdir, 'test_results_extended_C%i_%s.tsv' % (options.confidence, event_type))
        if options.verbose:
            print('Writing extended test results to %s' % out_fname)
        header_long = sp.r_[header, ['event_count:%s' % x for x in event_strains], ['gene_exp:%s' % x for x in event_strains], ['disp_raw', 'disp_adj']]

        data_out = sp.c_[data_out, (cov_used[s_idx, :] / sf).astype('str'), disp_raw_used[s_idx].astype('str'), disp_adj_used[s_idx].astype('str')]
        sp.savetxt(out_fname, sp.r_[header_long[sp.newaxis, :], data_out], delimiter='\t', fmt='%s')

        ### write output unique over genes
        out_fname = os.path.join(outdir, 'test_results_C%i_%s.gene_unique.tsv' % (options.confidence, event_type))
        if options.verbose:
            print('Writing gene unique test results to %s' % out_fname)
        ks_idx = []
        taken = set()
        for i in s_idx:
            gid = gene_ids[gene_idx[i]]
            if gid in taken:
                continue
            ks_idx.append(i)
            taken.add(gid)
        ks_idx = sp.array(ks_idx)

        data_out = sp.c_[event_ids[ks_idx], gene_ids[gene_idx[ks_idx]], pvals[ks_idx].astype('str'), pvals_adj[ks_idx].astype('str'), m_all[ks_idx, :]]
        data_out = sp.r_[header[sp.newaxis, :], data_out]
        sp.savetxt(out_fname, data_out, delimiter='\t', fmt='%s')

