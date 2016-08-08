import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import scipy as sp
import re

def mean_variance_plot(counts, disp, matrix, figtitle, filename, CFG):

    ### generate 
    fig = plt.figure(figsize=(10, 10), dpi=300)
    fig.suptitle(figtitle, fontsize=12)
    gs = gridspec.GridSpec(2, 2)
    idx = sp.where(~sp.isnan(disp))[0]
    idxa = sp.where(matrix[:, 1] & ~matrix[:, 3])[0]
    idxb = sp.where(~matrix[:, 1] & ~matrix[:, 3])[0]
    idxga = sp.where(matrix[:, 2] & matrix[:, 3])[0]
    idxgb = sp.where(~matrix[:, 2] & matrix[:, 3])[0]

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(sp.mean(sp.log10(counts[:, idxa] + 0.5), axis=1)[idx], sp.sqrt(disp[idx]), 'bo', alpha=0.3, markeredgecolor='none')
    ax.set_title('Events Condition 1')
    ax.set_xlabel('log10(Mean expression + 0.5)')
    ax.set_ylabel('Sqrt. dispersion')

    ax = fig.add_subplot(gs[0, 1])
    ax.plot(sp.mean(sp.log10(counts[:, idxb] + 0.5), axis=1)[idx], sp.sqrt(disp[idx]), 'bo', alpha=0.3, markeredgecolor='none')
    ax.set_title('Events Condition 2')
    ax.set_xlabel('log10(Mean expression + 0.5)')
    ax.set_ylabel('Sqrt. dispersion')

    ax = fig.add_subplot(gs[1, 0])
    ax.plot(sp.mean(sp.log10(counts[:, idxga] + 0.5), axis=1)[idx], sp.sqrt(disp[idx]), 'bo', alpha=0.3, markeredgecolor='none')
    ax.set_title('Gene Expression Condition 1')
    ax.set_xlabel('log10(Mean expression + 0.5)')
    ax.set_ylabel('Sqrt. dispersion')

    ax = fig.add_subplot(gs[1, 1])
    ax.plot(sp.mean(sp.log10(counts[:, idxgb] + 0.5), axis=1)[idx], sp.sqrt(disp[idx]), 'bo', alpha=0.3, markeredgecolor='none')
    ax.set_title('Gene Expression Condition 2')
    ax.set_xlabel('log10(Mean expression + 0.5)')
    ax.set_ylabel('Sqrt. dispersion')

    plt.savefig(filename, format='png', bbox_inches='tight')
    plt.close(fig)


def qq_plot(pvals, figtitle, filename, CFG):
    '''
    create a quantile quantile plot for the given p-values
    '''

    exp = sp.linspace(0, 1, num=pvals.shape[0])

    ### plot no log
    fig = plt.figure(figsize=(10, 10), dpi=100)
    fig.suptitle(figtitle, fontsize=12)
    ax = fig.add_subplot(111)
    ax.set_title(figtitle)
    ax.set_ylabel("Oberserved P-Value")
    ax.set_xlabel("Expected P-Value")
    ax.plot(exp, sp.sort(pvals), 'bo') 
    ax.set_xlim([0, 1.0])
    ax.set_ylim([0, 1.0])
    ax.plot([0, 1.0], [0, 1.0], 'r--')
    plt.savefig(filename, format='png', bbox_inches='tight')
    plt.close(fig)

    ### plot with log
    fig = plt.figure(figsize=(10, 10), dpi=100)
    fig.suptitle(figtitle, fontsize=12)
    ax = fig.add_subplot(111)
    ax.set_title(figtitle)
    ax.set_ylabel("Oberserved P-Value (-log10)")
    ax.set_xlabel("Expected P-Value (-log10)")
    ax.plot(-sp.log10(exp), -sp.log10(sp.sort(pvals)), 'bo') 
    maxlim = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.set_xlim([0, maxlim])
    ax.set_ylim([0, maxlim])
    ax.plot([0, maxlim], [0, maxlim], 'r--')
    plt.savefig(re.sub(r'.png$', '', filename) + '.log10.png', format='png', bbox_inches='tight')
    plt.close(fig)
