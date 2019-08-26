#! /usr/bin/env python
import sys
import os
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import pickle
import pdb

from . import settings
from .viz.graph import *
from .viz.coverage import *
from .viz.genelets import *
from .viz import axes as vax
from .helpers_viz import *

### intermediate fix to load pickle files stored under previous version
from .classes import gene as cgene
from .classes import splicegraph as csplicegraph
from .classes import segmentgraph as csegmentgraph
from .classes import datatrack as dt
sys.modules['modules.classes.gene'] = cgene
sys.modules['modules.classes.splicegraph'] = csplicegraph
sys.modules['modules.classes.segmentgraph'] = csegmentgraph


def _add_ax(axes, fig, gs):
    sharex = None if len(axes) == 0 else axes[0]
    axes.append(fig.add_subplot(gs[len(axes), 0], sharex=sharex))
    return axes[-1]

def spladder_viz(options):

    """Main visualization code"""

    ### parse parameters from options object
    options = settings.parse_args(options, identity='viz')

    ### create plot directory if it does not exist yet
    dirname = options.outdir if options.testdir == '-' else options.testdir
    if not os.path.exists(os.path.join(dirname, 'plots')):
        os.mkdir(os.path.join(dirname, 'plots'))

    if options.format == 'd3':
        try:
            import mpld3
            from mpld3 import plugins
        except ImportError:
            sys.stderr.write("ERROR: missing package for output format d3. Package mpld3 required")
            sys.exit(1)

    ### collect elements to be plotted

    ### genes
    ### load gene information
    gene_names = get_gene_names(options)
    gids = []
    for g in options.genes:
        gid = sp.where(gene_names == g)[0]
        if gid.shape[0] == 0:
            sys.stderr.write('ERROR: provided gene ID %s could not be found, please check for correctness\n' % g)
            sys.exit(1)
        assert gid.shape[0] == 1
        gids.append(gid[0])
    genes = load_genes(options, sp.array(gids))
    gene_names = sp.array([_.name for _ in genes])

    ### events
    eids = []
    for e in options.events:
        eid = e.rsplit('.', 1)
        if len(eid) == 1:
            eid = e.rsplit('_', 1)
            if len(eid) == 1:
                sys.stderr.write('ERROR: provided event ID %s could not be parsed, please check for correctness\n' % e)    
                sys.exit(1)
        eid[-1] = int(eid[-1]) - 1
        eids.append(eid)
    events = load_events(options, sp.array(eids))

    ### check that everthing is on the same chromosome
    plotchrm = sp.unique([_.chr for _ in sp.r_[events, genes]])
    if plotchrm.shape[0] > 1:
        sys.stderr.write('ERROR: the provided genes / events are on different chromosomes and canot be plotted jointly\n')
        sys.exit(1)

    ### identify the plotting range
    plotrange = None
    for g in genes:
        if not plotrange:
            plotrange = [g.start, g.stop]
        else:
            plotrange[0] = min(g.start, plotrange[0])
            plotrange[1] = min(g.stop, plotrange[1])
    for e in events:
        if not plotrange:
            plotrange = [e.exons2.min(), e.exons2.max()]
        else:
            plotrange[0] = min(e.exons2.min(), plotrange[0])
            plotrange[1] = max(e.exons2.max(), plotrange[1])

    if plotrange[1] - plotrange[0] > 1000000:
        sys.stderr.write('ERROR: plotting area has a width of more than 1 000 000 positions (%i) - aborting\n' % (plotrange[1] - plotrange[0]))
        sys.exit(1)

    ### data track
    data_tracks = []
    ### parse all elements to be plotted as data tracks
    for data_element in options.data_tracks:
        track_types = data_element[0].split(',')
        for track_type in track_types:
            data_tracks.append(dt.DataTrack(track_type, data_element[1:]))

    ### set color maps
    cmap_cov = plt.get_cmap('jet')
    cmap_edg = plt.get_cmap('jet')

    ### plot log scale?
    log_tag = ''
    if options.log:
        log_tag = '.log'
    event_tag = ''

    gs = gridspec.GridSpec(len(data_tracks), 1)
    fig = plt.figure(figsize = (10, 3 * len(data_tracks)), dpi=200)

    axes = []
    for i,data_track in enumerate(data_tracks):

        ### plot splicing graph
        if data_track.type == 'splicegraph':
            ax = _add_ax(axes, fig, gs)
            for gene in genes:
                gene.from_sparse()
                plot_graph(gene.splicegraph.vertices, gene.splicegraph.edges, ax, xlim=plotrange)
                gene.to_sparse()
            ax.set_title('Splicing graph for %s' % ','.join(gene_names))
        ### plot annotated transcripts
        if data_track.type == 'transcripts':
            ax = _add_ax(axes, fig, gs)
            for gene in genes:
                gene.from_sparse()
                multiple(gene.exons, ax=ax, x_range=plotrange)
                gene.to_sparse()
            ax.set_title('Annotated transcripts for %s' % ','.join(gene_names))
        ### plot coverage tracks
        if data_track.type == 'coverage':
            ax = _add_ax(axes, fig, gs)
            min_sample_size = 5 ### TODO make that an option
            caxes = []
            for gene in genes:
                caxes = []
                labels = []
                for g, bam_group in enumerate(data_track.bam_fnames):
                    label = 'group %i' % (g + 1) if len(data_track.group_labels) == 0 else data_track.group_labels[g]
                    norm = plt.Normalize(0, len(bam_group))
                    caxes.append(cov_from_bam(gene.chr, 
                                              plotrange[0], 
                                              plotrange[1], 
                                              bam_group, 
                                              subsample=min_sample_size, 
                                              ax=ax, 
                                              intron_cnt=True,
                                              log=options.log, 
                                              title='Expression', 
                                              xlim=plotrange, 
                                              color_cov=cmap_cov(norm(g)), 
                                              color_intron_edge= cmap_edg(norm(g)),
                                              grid=True, 
                                              min_intron_cnt=options.mincount, 
                                              return_legend_handle=True, 
                                              label=label))
                    labels.append(label)
            if len(caxes) > 0:
                ax.legend(caxes, labels)
        ### plot segment counts
        if data_track.type == 'segments':
            ax = _add_ax(axes, fig, gs)
            for g, gid in enumerate(gids):
                (segments, edges, edge_idx, strains) = get_seg_counts(gid, options.outdir, options.confidence, options.validate_sg)
                seg_sample_idx = None
                if len(data_track.strains) > 0:
                    seg_sample_idx = []
                    for track_strains in data_track.strains:
                        seg_sample_idx.append(sp.where(sp.in1d(strains, track_strains))[0])
                        cov_from_segments(genes[g], segments, edges, edge_idx, ax, xlim=plotrange, log=options.log, grid=True, order='C', sample_idx=seg_sample_idx)
            ax.set_title('Segment counts')

    ### the plotting happens on the results of spladder test
    ### the user chooses to plot the top k significant events
    ### this requires the event type to be specified
    #if options.test_result > 0:
    #    gene_names = []
    #    for event_type in options.event_types:
    #        ### the testing script should generate a setup file for the test
    #        ### SETUP is structured as follows:
    #        ###  [gene_strains, event_strains, dmatrix0, dmatrix1, event_type, options]
    #        labels = options.test_labels.split(':')
    #        options.labels = labels
    #        if options.testdir != '-':
    #            testdir = dirname
    #        else:
    #            testdir = os.path.join(dirname, 'testing_%s_vs_%s' % (labels[0], labels[1]))
    #        SETUP = pickle.load(open(os.path.join(testdir, 'test_setup_C%i_%s.pickle' % (options.confidence, event_type)), 'rb'))

    #        ### get strains to plot
    #        idx1 = sp.where(sp.in1d(SETUP[0], SETUP[6]['conditionA']))[0]
    #        idx2 = sp.where(sp.in1d(SETUP[0], SETUP[6]['conditionB']))[0]

    #        ### load test results
    #        for l, line in enumerate(open(os.path.join(testdir, 'test_results_C%i_%s.tsv' % (options.confidence, event_type)), 'r')):
    #            if l == 0:
    #                continue
    #            if l > options.test_result:
    #                break
    #            sl = line.strip().split('\t')
    #            gene_names.append([sl[1], sl[0]])
    #    gids = get_gene_ids(options, gene_names)
    #### no gene specified but result provided - plot all genes with confirmed events
    #### if an event_id is provided, only the associated gene will be plotted
    #else:
    #    gids = get_gene_ids(options)


    ### iterate over genes to plot
    #for gid in gids:
    #    ### gather information about the gene we plot
    #    gene = load_genes(options, idx=[gid[0]])[0]
    #    if options.verbose:
    #        print('plotting information for gene %s' % gene.name)
    #    gene.from_sparse()
    #    if gid[1] is not None:
    #        event_info = [x[::-1] for x in re.split(r'[._]', gid[1][::-1], maxsplit=1)[::-1]]
    #        event_info[1] = int(event_info[1]) - 1
    #        event_info = sp.array(event_info, dtype='str')[sp.newaxis, :]
    #        event_tag = '.%s' % gid[1]
    #    ### get all confident events of the current gene
    #    else:
    #        event_info = get_conf_events(options, gid[0])

    #    ### go over different plotting options
    #    axes = []
    #    ### plot result of testing
    #    if options.test_result > 0:
    #        fig = plt.figure(figsize = (9, 5), dpi=200)
    #        gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
    #        _add_ax(fig, axes, gs)
    #        _add_ax(fig, axes, gs)
    #        _plot_event(options, event_info, fig, axes[1], gs, None, padding=100)
    #        start, stop = axes[1].get_xlim()
    #        plot_bam(options, gene, options.bam_fnames, fig, axes[0], gs, None, cmap_cov, cmap_edg, single=False, sharex=axes[1], start=int(start), stop=int(stop))

    #    ### plot custom layout
    #    else:
    #        ### set defaults
    #        if not options.user:
    #            options.splicegraph = True
    #            options.transcripts = True

    #        if options.format == 'd3':
    #            fig = plt.figure(figsize = (12, 2*rows), dpi=100)
    #        else:
    #            fig = plt.figure(figsize = (18, 3*rows), dpi=200)

    #        xlim = None
    #        ### plot structure of a single given event
    #        _plot_event(options, event_info, fig, axes, gs, xlim)

    #    ### we only need to adapt the zoom for one axis object - as we share the x
    #    zoom_x = [float(x) for x in options.zoom_x.split(',')]
    #    xlim = axes[0].get_xlim()
    #    xdiff = xlim[1] - xlim[0]
    #    axes[0].set_xlim([xlim[0] + (zoom_x[0] * xdiff), xlim[0] + (zoom_x[1] * xdiff)])

    for ax in axes:
        vax.clean_axis(ax)

    plt.tight_layout()
    ### save plot into file
    if options.format == 'd3':
        out_fname = os.path.join(dirname, 'plots', 'gene_overview_C%i_%s%s%s.html' % (options.confidence, ''.join(gene_names), event_tag, log_tag))
        plugins.clear(fig)
        plugins.connect(fig, plugins.Zoom(enabled=True))
        mpld3.save_html(fig, open(out_fname, 'w'))
    else:
        if options.test_result > 0:
            out_fname = os.path.join(dirname, 'plots', 'gene_overview_C%i_%s%s%s.%s' % (options.confidence, ''.join(gene_names), event_tag, log_tag, options.format))
        else:
            out_fname = os.path.join(dirname, 'plots', 'gene_overview_C%i_%s%s%s.%s' % (options.confidence, ''.join(gene_names), event_tag, log_tag, options.format))
        plt.savefig(out_fname, format=options.format, bbox_inches='tight')
    plt.close(fig)


def plot_bam(options, gene, samples, labels, fig, axes, gs, xlim, cmap_cov, cmap_edg, single=True, subsample_size=5, x_range=None):

    min_sample_size = min(subsample_size, min([len(x) for x in samples]))
    if x_range is None:
        x_range = [gene.splicegraph.vertices.min(), gene.splicegraph.vertices.max()]

    norm = plt.Normalize(0, len(samples))
    caxes = []
    labels = []

    for s, bams in enumerate(samples):

        if len(labels) > 0:
            label = options.labels[s]
        else:
            label = 'group %i' % (s + 1)

        if single:
            axes.append(fig.add_subplot(gs[len(axes), 0], sharex=sharex))
            title = 'Expression (%s)' % label
            color_cov = cmap_cov(norm(0)) # '#d7191c'
            color_edg = cmap_edg(norm(0)) # '#1a9641'
            ax = axes[-1]
        else:
            if s == 0:
                if hasattr(axes, '__iter__'):
                    axes.append(fig.add_subplot(gs[len(axes), 0], sharex=sharex))
                    ax = axes[-1]
                else:
                    ax = axes
            title = 'Expression for all sample groups'
            color_cov = cmap_cov(norm(s))
            color_edg = cmap_edg(norm(s))

        caxes.append(cov_from_bam(gene.chr, start, stop, bams, subsample=min_sample_size, ax=ax, intron_cnt=True,
                     log=options.log, title=title, xlim=xlim, color_cov=color_cov, color_intron_edge=color_edg,
                     grid=True, min_intron_cnt=options.mincount, return_legend_handle=True, label=label))
        labels.append(label)
        ax.set_xlabel('')

    ax.legend(caxes, labels)


def _plot_event(options, event_info, fig, axes, gs, xlim, padding=None):
    """This function takes the event_id given in the options object and 
    plots it into ax."""

    axes.append(fig.add_subplot(gs[len(axes), 0]))
    event_list = [_ for event in load_events(options, event_info) for _ in [event.exons1, event.exons2]]
    multiple(event_list, ax=axes[-1], x_range=xlim, color='green', padding=padding) 
    #ax.set_title('Alt event structure') # of %s' % options.event_id)
    vax.clean_axis(axes[-1], allx=True)


