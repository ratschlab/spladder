import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import scipy as sp
import re

def mean_variance_plot(counts, disp, matrix, figtitle, filename, options):

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

    plt.savefig(filename, format=options.plot_format, bbox_inches='tight')
    plt.close(fig)


def qq_plot(pvals, figtitle, filename, options):
    '''
    create a quantile quantile plot for the given p-values
    '''

    idx = sp.where(~sp.isnan(pvals))[0]
    exp = sp.linspace(0, 1, num=idx.shape[0])

    ### plot no log
    fig = plt.figure(figsize=(10, 10), dpi=100)
    fig.suptitle(figtitle, fontsize=12)
    ax = fig.add_subplot(111)
    ax.set_title(figtitle)
    ax.set_ylabel("Oberserved P-Value")
    ax.set_xlabel("Expected P-Value")
    ax.plot(exp, sp.sort(pvals[idx]), 'bo') 
    ax.set_xlim([0, 1.0])
    ax.set_ylim([0, 1.0])
    ax.plot([0, 1.0], [0, 1.0], 'r--')
    plt.savefig(filename, format=options.plot_format, bbox_inches='tight')
    plt.close(fig)

    ### plot with log
    fig = plt.figure(figsize=(10, 10), dpi=100)
    fig.suptitle(figtitle, fontsize=12)
    ax = fig.add_subplot(111)
    ax.set_title(figtitle)
    ax.set_ylabel("Oberserved P-Value (-log10)")
    ax.set_xlabel("Expected P-Value (-log10)")
    eps = 10e-5
    ax.plot(-sp.log10(exp + eps), -sp.log10(sp.sort(pvals[idx] + eps)), 'bo') 
    maxlim = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.set_xlim([0, maxlim])
    ax.set_ylim([0, maxlim])
    ax.plot([0, maxlim], [0, maxlim], 'r--')
    plt.savefig(re.sub(r'.%s$' % options.plot_format, '', filename) + '.log10.%s' % options.plot_format, format=options.plot_format, bbox_inches='tight')
    plt.close(fig)


def count_histogram(counts, matrix, figtitle, filename, options):
    '''
    create a histogram plot showing the count distributions
    '''

    if options.verbose:
        print('Plotting count distributions')

    ### generate 
    fig = plt.figure(figsize=(15, 10), dpi=300)
    fig.suptitle(figtitle, fontsize=12)
    gs = gridspec.GridSpec(2, 2)
    idxa = sp.where(matrix[:, 1] & ~matrix[:, 3])[0]
    idxb = sp.where(~matrix[:, 1] & ~matrix[:, 3])[0]
    idxga = sp.where(matrix[:, 2] & matrix[:, 3])[0]
    idxgb = sp.where(~matrix[:, 2] & matrix[:, 3])[0]

    ax = fig.add_subplot(gs[0, 0])
    ax.hist(counts[:, idxa].ravel(), 50, facecolor='blue', alpha=0.7)
    ax.set_title('Event counts Condition 1')
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Expression bin')

    ax = fig.add_subplot(gs[0, 1])
    ax.hist(counts[:, idxb].ravel(), 50, facecolor='blue', alpha=0.7)
    ax.set_title('Event counts Condition 2')
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Expression bin')

    ax = fig.add_subplot(gs[1, 0])
    ax.hist(counts[:, idxga].ravel(), 50, facecolor='blue', alpha=0.7)
    ax.set_title('Gene Expression Counts Condition 1')
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Expression bin')

    ax = fig.add_subplot(gs[1, 1])
    ax.hist(counts[:, idxgb].ravel(), 50, facecolor='blue', alpha=0.7)
    ax.set_title('Gene Expression Counts Condition 2')
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Expression bin')

    plt.savefig(filename, format=options.plot_format, bbox_inches='tight')
    plt.close(fig)

    ### generate 
    fig = plt.figure(figsize=(15, 10), dpi=300)
    fig.suptitle(figtitle, fontsize=12)

    ax = fig.add_subplot(gs[0, 0])
    ax.hist(sp.log10(counts[:, idxa].ravel() + 0.5), 50, facecolor='blue', alpha=0.7)
    ax.set_title('Event counts Condition 1')
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Expression bin (log10)')

    ax = fig.add_subplot(gs[0, 1])
    ax.hist(sp.log10(counts[:, idxb].ravel() + 0.5), 50, facecolor='blue', alpha=0.7)
    ax.set_title('Event counts Condition 2')
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Expression bin (log10)')

    ax = fig.add_subplot(gs[1, 0])
    ax.hist(sp.log10(counts[:, idxga].ravel() + 0.5), 50, facecolor='blue', alpha=0.7)
    ax.set_title('Gene Expression Counts Condition 1')
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Expression bin (log10)')

    ax = fig.add_subplot(gs[1, 1])
    ax.hist(sp.log10(counts[:, idxgb].ravel() + 0.5), 50, facecolor='blue', alpha=0.7)
    ax.set_title('Gene Expression Counts Condition 2')
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Expression bin (log10)')

    plt.savefig(re.sub(r'.%s$' % options.plot_format, '', filename) + '.log10.%s' % options.plot_format, format=options.plot_format, bbox_inches='tight')
    plt.close(fig)

def ma_plot(pvals, counts, fc, figtitle, filename, options, alpha=0.05):
    '''
    create an MA plot summarizing coverage, log fold changes and significant values
    '''

    fig = plt.figure(figsize=(10, 10), dpi=100)
    ax = fig.add_subplot(111)
    ax.set_title(figtitle)
    ax.set_ylabel("Fold change (log2)")
    ax.set_xlabel("Mean normalized counts (log 10)")
    idx = sp.where(pvals > alpha)[0]
    ax.plot(sp.log10(counts[idx] + 1), fc[idx], 'ko')
    idx = sp.where(pvals <= alpha)[0]
    ax.plot(sp.log10(counts[idx] + 1), fc[idx], 'ro')
    plt.savefig(filename, format=options.plot_format, bbox_inches='tight')
    plt.close(fig)


