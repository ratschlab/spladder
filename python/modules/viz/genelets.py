import matplotlib.pyplot as plt
import matplotlib
import scipy as sp

def single(exons, ax=None, x_range=None, count=1, color='blue', grid=False, label=None, padding=None):

    if ax is None:
        fig = plt.figure(figsize = (10, 4))
        ax = fig.add_subplot(111)
    
    ### draw grid
    if grid:
        ax.grid(b=True, which='major', linestyle='--', linewidth=0.2, color='#222222')
        ax.yaxis.grid(False)

    min_ex = None
    max_ex = None
    if len(exons.shape) == 1:
        exons = exons[sp.newaxis, :]

    for e_idx in range(exons.shape[0]):
        exon = exons[e_idx, :]
        if min_ex is None:
            min_ex = exon[0]
            max_ex = exon[1]
        else:
            min_ex = min(min_ex, exon[0])
            max_ex = max(max_ex, exon[1] + 1)
        if not padding is None:
            min_ex = max(min_ex - padding, 0)
            max_ex += padding
        rect = matplotlib.patches.Rectangle((exon[0], -20 - ((count - 1) * 25)), exon[1] - exon[0] + 1, 10, facecolor=color, edgecolor='none', alpha=0.7)
        ax.add_patch(rect)
        if e_idx > 0:
            ax.plot([exons[e_idx - 1, 1] + 1, exons[e_idx, 0]], [-15 - ((count - 1) * 25), -15 - ((count - 1) * 25)], 'b-')

    ### set genelet label
    if not label is None:
        ax.text(exons[0][0] + 10, -5 - ((count - 1) * 25), label)

    ax.set_ylim(((count * -25) - 5, 0))

    if x_range is None:
        ax.set_xlim((min_ex, max_ex))
    else:
        ax.set_xlim(x_range)

    ax.get_yaxis().set_visible(False)

def multiple(exon_set, ax=None, x_range=None, color='blue', labels=None, grid=None, padding=None):

    if labels is not None:
        assert len(labels) == len(exon_set), 'Number of labels and given genelet entities must be identical!\nlen(labels) = %i != len(entities) = %i' % (len(labels), len(exon_set))

    for e_idx, exons in enumerate(exon_set):
        if labels is not None:
            single(exons, ax=ax, x_range=x_range, color=color, count=e_idx + 1, grid=grid, label=labels[e_idx], padding=padding)
        else:
            single(exons, ax=ax, x_range=x_range, color=color, count=e_idx + 1, grid=grid, padding=padding)

