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


def get_plot_len(options):
    """Identifies the number of rows we need in our plot"""

    rows = 3 # splicing graph + events + segments
    if len(options.bam_fnames) > 0:
        rows += len(options.bam_fnames)
        if len(options.bam_fnames) > 1:
            rows += 1
    rows += int(options.transcripts)

    return rows


def _add_ax(fig, axes, gs):
    sharex = None if len(axes) == 0 else axes[0]
    axes.append(fig.add_subplot(gs[len(axes), 0], sharex=sharex))

def spladder_viz(options):

    """Main visualization code"""

    ### parse parameters from options object
    options = settings.parse_args(options, identity='viz')

    ### create plot directory if it does not exist yet
    if options.testdir != '-':
        dirname = options.testdir
    else:
        dirname = options.outdir
    if not os.path.exists(os.path.join(dirname, 'plots')):
        os.mkdir(os.path.join(dirname, 'plots'))

    if options.format == 'd3':
        try:
            import mpld3
            from mpld3 import plugins
        except ImportError:
            sys.stderr.write("ERROR: missing package for output format d3. Package mpld3 required")
            sys.exit(1)

    ### load gene information
    gene_names = get_gene_names(options)

    rows = get_plot_len(options)
    gs = gridspec.GridSpec(rows, 1)

    ### set color maps
    cmap_cov = plt.get_cmap('jet')
    cmap_edg = plt.get_cmap('jet')

    ### plot log scale?
    log_tag = ''
    if options.log:
        log_tag = '.log'
    event_tag = ''

    ### did we get any labels?
    if ',' in options.labels:
        options.labels = options.labels.strip(',').split(',')
        assert len(options.labels) == len(options.bam_fnames), "The number of given labels (%i) needs to match the number of given bam file groups (%i)" % (len(options.labels), len(options.bam_fnames))


    # Collect genes to be plotted - there are three cases possible
    # 1) the user provides a gene ID
    # 2) all genes containing a significant event of a given type from differential testing are plotted
    # 3) all genes that contain any event are plotted

    ### the user chose a specific gene for plotting
    ### create pairs of gene ids and an event_id (the latter is None by default)
    if options.gene_name is not None:
        gids = [[sp.where(sp.array(gene_names) == options.gene_name)[0][0], options.event_id]]
        if len(gids) == 0:
            sys.stderr.write('ERROR: provided gene ID %s could not be found, please check for correctness\n' % options.gene_name)
            sys.exit(1)
    ### the plotting happens on the results of spladder test
    ### the user chooses to plot the top k significant events
    ### this requires the event type to be specified
    elif options.test_result > 0:
        gene_names = []
        for event_type in options.event_types:
            ### the testing script should generate a setup file for the test
            ### SETUP is structured as follows:
            ###  [gene_strains, event_strains, dmatrix0, dmatrix1, event_type, options]
            labels = options.test_labels.split(':')
            options.labels = labels
            if options.testdir != '-':
                testdir = dirname
            else:
                testdir = os.path.join(dirname, 'testing_%s_vs_%s' % (labels[0], labels[1]))
            SETUP = pickle.load(open(os.path.join(testdir, 'test_setup_C%i_%s.pickle' % (options.confidence, event_type)), 'rb'))

            ### get strains to plot
            idx1 = sp.where(sp.in1d(SETUP[0], SETUP[6]['conditionA']))[0]
            idx2 = sp.where(sp.in1d(SETUP[0], SETUP[6]['conditionB']))[0]

            ### load test results
            for l, line in enumerate(open(os.path.join(testdir, 'test_results_C%i_%s.tsv' % (options.confidence, event_type)), 'r')):
                if l == 0:
                    continue
                if l > options.test_result:
                    break
                sl = line.strip().split('\t')
                gene_names.append([sl[1], sl[0]])
        gids = get_gene_ids(options, gene_names)
    ### no gene specified but result provided - plot all genes with confirmed events
    ### if an event_id is provided, only the associated gene will be plotted
    else:
        gids = get_gene_ids(options)


    ### iterate over genes to plot
    for gid in gids:
        ### gather information about the gene we plot
        gene = load_genes(options, idx=[gid[0]])[0]
        if options.verbose:
            print('plotting information for gene %s' % gene.name)
        gene.from_sparse()

        ### event to plot is specified with the gene id list
        if gid[1] is not None:
            event_info = [x[::-1] for x in re.split(r'[._]', gid[1][::-1], maxsplit=1)[::-1]]
            event_info[1] = int(event_info[1]) - 1
            event_info = sp.array(event_info, dtype='str')[sp.newaxis, :]
            event_tag = '.%s' % gid[1]
        ### get all confident events of the current gene
        else:
            event_info = get_conf_events(options, gid[0])

        ### go over different plotting options
        axes = []
        ### plot result of testing
        if options.test_result > 0:
            fig = plt.figure(figsize = (9, 5), dpi=200)
            gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
            _add_ax(fig, axes, gs)
            _add_ax(fig, axes, gs)
            _plot_event(options, event_info, fig, axes[1], gs, None, padding=100)
            start, stop = axes[1].get_xlim()
            plot_bam(options, gene, options.bam_fnames, fig, axes[0], gs, None, cmap_cov, cmap_edg, single=False, sharex=axes[1], start=int(start), stop=int(stop))

        ### plot custom layout
        else:
            ### set defaults
            if not options.user:
                options.splicegraph = True
                options.transcripts = True

            if options.format == 'd3':
                fig = plt.figure(figsize = (12, 2*rows), dpi=100)
            else:
                fig = plt.figure(figsize = (18, 3*rows), dpi=200)

            xlim = None
            ### plot splicing graph
            if options.splicegraph:
                _plot_splicegraph(gene, fig, axes, gs)
                xlim = axes[-1].get_xlim()

            ### plot annotated transcripts
            if options.transcripts:
                sharex = None if len(axes) == 0 else axes[0]
                axes.append(fig.add_subplot(gs[len(axes), 0], sharex=sharex))
                multiple(gene.exons, ax=axes[-1], x_range=xlim)
                axes[-1].set_title('Annotated Transcripts')

            ### plot coverage information for a set of given samples
            if len(options.bam_fnames) > 0:
                plot_bam(options, gene, options.bam_fnames, fig, axes, gs, xlim, cmap_cov, cmap_edg)

                ### plot all the samples in a single plot
                if len(options.bam_fnames) > 1:
                    plot_bam(options, gene, options.bam_fnames, fig, axes, gs, xlim, cmap_cov, cmap_edg, single=False)

            ### plot segment counts
            if len(options.bam_fnames) == 0 or False: # add option for segment plots
                if options.test_result > 0:
                    _plot_segments(options, gid, fig, axes, gs, [idx1, idx2])
                else:
                    _plot_segments(options, gid, fig, axes, gs)

            ### plot structure of a single given event
            _plot_event(options, event_info, fig, axes, gs, xlim)

        ### we only need to adapt the xoom for one axis object - as we share the x
        zoom_x = [float(x) for x in options.zoom_x.split(',')]
        xlim = axes[0].get_xlim()
        xdiff = xlim[1] - xlim[0]
        axes[0].set_xlim([xlim[0] + (zoom_x[0] * xdiff), xlim[0] + (zoom_x[1] * xdiff)])

        for ax in axes:
            vax.clean_axis(ax)

        plt.tight_layout()
        ### save plot into file
        if options.format == 'd3':
            out_fname = os.path.join(dirname, 'plots', 'gene_overview_C%i_%s%s%s.html' % (options.confidence, gene.name, event_tag, log_tag))
            plugins.clear(fig)
            plugins.connect(fig, plugins.Zoom(enabled=True))
            mpld3.save_html(fig, open(out_fname, 'w'))
        else:
            if options.test_result > 0:
                out_fname = os.path.join(dirname, 'plots', 'gene_overview_C%i_%s%s%s.%s' % (options.confidence, gene.name, event_tag, log_tag, options.format))
            else:
                out_fname = os.path.join(dirname, 'plots', 'gene_overview_C%i_%s%s%s.%s' % (options.confidence, gene.name, event_tag, log_tag, options.format))
            plt.savefig(out_fname, format=options.format, bbox_inches='tight')
        plt.close(fig)


def plot_bam(options, gene, samples, fig, axes, gs, xlim, cmap_cov, cmap_edg, single=True, subsample_size=5, sharex=None, start=None, stop=None):

    min_sample_size = min(subsample_size, min([len(x) for x in samples]))
    if start is None:
        start = gene.splicegraph.vertices.min()
    if stop is None:
        stop = gene.splicegraph.vertices.max()

    norm = plt.Normalize(0, len(samples))
    caxes = []
    labels = []

    if sharex is None:
        sharex = None if len(axes) == 0 else axes[0]
    for s, bams in enumerate(samples):

        if options.labels != '':
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

    if not single:
        ax.legend(caxes, labels)


def _plot_event(options, event_info, fig, axes, gs, xlim, padding=None):
    """This function takes the event_id given in the options object and 
    plots it into ax."""

    axes.append(fig.add_subplot(gs[len(axes), 0]))
    event_list = [_ for event in load_events(options, event_info) for _ in [event.exons1, event.exons2]]
    multiple(event_list, ax=axes[-1], x_range=xlim, color='green', padding=padding) 
    #ax.set_title('Alt event structure') # of %s' % options.event_id)
    vax.clean_axis(axes[-1], allx=True)


def _plot_splicegraph(gene, fig, axes, gs):
    """Append a new object to the axes list, and into the gridspec
       at the current position and plot a splicing graph"""

    axes.append(fig.add_subplot(gs[len(axes), 0]))
    plot_graph(gene.splicegraph.vertices, gene.splicegraph.edges, axes[-1])
    axes[-1].set_title('Splicing graph for %s' % gene.name)


def _plot_segments(options, gid, fig, axes, gs, seg_sample_idx=None):

    print('get segment counts')
    (segments, edges, edge_idx, strains) = get_seg_counts(options, gid[0])
    seg_sample_idx = None
    if len(options.strains) > 0:
        seg_sample_idx = []
        for group in options.strains:
            seg_sample_idx.append(sp.where(sp.in1d(strains, group))[0])
    axes.append(fig.add_subplot(gs[len(axes), 0], sharex=axes[0]))
    print('plot segment counts')
    cov_from_segments(gene, segments, edges, edge_idx, axes[-1], xlim=xlim, log=options.log, grid=True, order='C', sample_idx=seg_sample_idx)
    axes[-1].set_title('Segment counts')

