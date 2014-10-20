import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cPickle
import sys
import scipy as sp
import pdb

def plot_graph(splicegraph, ax, xlim=None, highlight=None, highlight_color='magenta'):
    """Takes a splicing graph and visualizes its structure"""

    start = splicegraph.vertices.ravel().min()
    stop = splicegraph.vertices.ravel().max()

    ### draw grid
    ax.grid(b=True, which='major', linestyle='--', linewidth=0.2, color='#222222')
    ax.yaxis.grid(False)

    ### nodes
    nodes = []
    exon_num = sp.zeros((stop - start,))
    exon_loc = sp.zeros((1, stop - start))
    exon_level = sp.zeros((splicegraph.vertices.shape[1], splicegraph.vertices.shape[1]))
    for i in range(splicegraph.vertices.shape[1]):
        cur_vertex = splicegraph.vertices[:, i] - start
        exon_num[cur_vertex[0]:cur_vertex[1]] += 1
        if sp.all(exon_num < 2):
            exon_loc[0, :] = exon_num
            level = 0
        elif exon_num.max() > exon_loc.shape[0]:
            exon_loc = sp.r_[exon_loc, sp.zeros((1, stop - start))]
            exon_loc[-1, cur_vertex[0]:cur_vertex[1]] = 1 
            level = exon_loc.shape[0] - 1
        elif exon_num.max() <= exon_loc.shape[0]:
            idx = sp.where(sp.all(exon_loc[:, cur_vertex[0]:cur_vertex[1]] == 0, 1))[0].min() 
            exon_loc[idx, cur_vertex[0]:cur_vertex[1]] = 1
            level = idx
       
        exon_level[i] = level
        
        nodes.append(matplotlib.patches.Rectangle([cur_vertex[0] + start, 20 + (level * 20)], cur_vertex[1] - cur_vertex[0], 10, facecolor='b', edgecolor='none', alpha=0.7))


    ### edges
    intron_loc = sp.zeros((1, stop - start))
    if splicegraph.edges.shape[0] > 1:
        for i in range(splicegraph.vertices.shape[1]):
            for j in range(i + 1, splicegraph.vertices.shape[1]):
                if splicegraph.edges[i ,j] > 0:
                    if splicegraph.vertices[0,i] < splicegraph.vertices[0,j]:
                        istart = splicegraph.vertices[1, i]
                        istop = splicegraph.vertices[0, j]
                        level1 = exon_level[i]
                        level2 = exon_level[j]
                    else:
                        istart = splicegraph.vertices[1, j]
                        istop = splicegraph.vertices[0, i]
                        level1 = exon_level[j]
                        level2 = exon_level[i]
              
                    cur_intron = [istart - start, istop - start]
                    intron_loc[cur_intron[0]:cur_intron[1]] += 1
                    leveli = [(istart + istop) * 0.5, (level1 + level2) * 0.5]
                    ax.plot([istart, leveli[0]], [25 + (level1 * 20), 32 + (leveli[1] * 20)], '-', color='#999999', linewidth=0.5)
                    ax.plot([leveli[0], istop], [32 + (leveli[1] * 20), 25 + (level2 * 20)], '-', color='#999999', linewidth=0.5)
                    #ax.plot([istart, istop], [20 + (level1 * 20) + 5, 20 + (level2 * 20) + 5], 'b-')

    ### draw nodes 
    for node in nodes:
        ax.add_patch(node)

    ### axes 
    if xlim is not None:
        ax.set_xlim(xlim)
    else:
        ax.set_xlim([start - 100, stop + 100])
    ax.set_ylim([0, 40 + (exon_loc.shape[0] * 20)]) 
    ax.set_yticklabels([])

    ### highlight if requested
    if highlight is not None:
        rect = patches.Rectangle((highlight[0], 0), highlight[1] - highlight[0], ax.get_ylim()[1], facecolor=highlight_color, edgecolor='none', alpha=0.5)
        ax.add_patch(rect)
        

