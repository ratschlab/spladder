"""
Negative binomial likelihood function.

authors: Yi Zhong
"""

import scipy as sp
from scipy.stats import nbinom
from scipy.special import digamma

def adj_loglikelihood_scalar(disp, X, y, mu, sign):

    n = 1 / disp
    p = n / (n + mu)
    loglik = sum(nbinom.logpmf(y, n, p))

    diagVec = mu / (1 + mu * disp)
    diagWM = sp.diag(diagVec)
    xtwx = sp.dot(sp.dot(X.T, diagWM), X)
    coxreid = 0.5 * sp.log(sp.linalg.det(xtwx))

    return (loglik - coxreid) * sign
