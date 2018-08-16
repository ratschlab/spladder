"""
Negative binomial likelihood function.

authors: Yi Zhong
"""

import scipy as sp
from scipy.stats import nbinom
from scipy.special import digamma


class complexException(Exception):
    pass

def adj_loglikelihood(xVec, lenSampleRibo, lenSampleRna, X, y, mu, sign):

    disp = sp.hstack([sp.repeat(xVec[0], lenSampleRibo), sp.repeat(xVec[1], lenSampleRna)])
    n = 1 / disp
    p = n / (n + mu)
    loglik = sum(nbinom.logpmf(y, n, p))

    diagVec = mu / (1 + sp.dot(mu.transpose(), disp))
    diagWM = sp.diagflat(diagVec)
    xtwx = sp.dot(sp.dot(sp.transpose(X), diagWM), X)

    coxreid = 0.5 * sp.log(sp.linalg.det(xtwx))
    ret = (loglik - coxreid) * sign
    #print "return value is " + str(ret)
    if isinstance(ret, complex):
        raise complexException()
    return ret


def adj_loglikelihood_gradient(xVec, lenSampleRibo, lenSampleRna, X, y, mu, sign):

    disp = sp.hstack([sp.repeat(xVec[0], lenSampleRibo), sp.repeat(xVec[1], lenSampleRna)])
    Gradient = sp.zeros_like(xVec)

    for i in range(len(xVec)):
        f1 = (digamma((1 / xVec[i])) - digamma( (1 / xVec[i]) + y[i])) / (xVec[i] ** 2)
        f2 = -((xVec[i] * mu[i] + (1 + xVec[i] * mu[i]) * sp.log(1 / (1 + xVec[i] * mu[i]))) / ((xVec[i] ** 2) * (1 + xVec[i] * mu[i])))
        f3 = y[i] / (xVec[i] + (xVec[i] ** 2) * mu[i])
        f4 = 0.5 * X.shape[1] * (mu[i] / (1 + sp.dot(mu.transpose(), disp)))
        Gradient[i] = f1 + f2 + f3 + f4

    return Gradient


def adj_loglikelihood_scalar(disp, X, y, mu, sign):

    n = 1 / disp
    p = n / (n + mu)
    loglik = sum(nbinom.logpmf(y, n, p))

    diagVec = mu / (1 + mu * disp)
    diagWM = sp.diag(diagVec)
    xtwx = sp.dot(sp.dot(X.T, diagWM), X)
    coxreid = 0.5 * sp.log(sp.linalg.det(xtwx))

    ret = (loglik - coxreid) * sign
    if isinstance(ret, complex):
        raise complexException()

    return ret
