#! /usr/bin/env python
import sys
import os
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('seaborn')
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import pickle
import copy
from collections import namedtuple

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

EVENT_TYPES = ['exon_skip', 'intron_retention', 'alt_3prime', 'alt_5prime', 'mult_exon_skip', 'mutex_exons']
def _add_gene_events(id_array, gene_id, event_type, outdir, confidence):
    
    if event_type == 'any':
        curr_event_types = EVENT_TYPES
    elif not event_type in EVENT_TYPES:
        sys.stderr.write('ERROR: Given event type (%s) is not in the list of allowed event types (%s)\n' % (event_type, ','.join(EVENT_TYPES)))
        sys.exit(1)
    else:
        curr_event_types = [event_type]

    for t in curr_event_types:
        id_array.extend([[t, _] for _ in get_event_ids_from_gene(gene_id, t, outdir, confidence)])

def _add_ax(axes, fig, gs):
    sharex = None if len(axes) == 0 else axes[0]
    axes.append(fig.add_subplot(gs[len(axes), 0])) 
    return axes[-1]

def _parse_event_info(id_array, gids, event_info, events, outdir, confidence, verbose):

    for e in event_info:
        ### if there is no number in the event name, we assume all events of that type
        if len([_ for _ in re.split(r'[._]', e) if _.isdigit()]) == 0:
            ### get all relevant events from each gene
            for gene_id in gids:
                _add_gene_events(id_array, gene_id, e, outdir, confidence)
        else:
            ### normally parse event ID
            eid = e.rsplit('.', 1)
            if len(eid) == 1:
                eid = e.rsplit('_', 1)
                if len(eid) == 1:
                    sys.stderr.write('ERROR: provided event ID "%s" could not be parsed, please check for correctness\n' % e)    
                    sys.exit(1)
            eid[-1] = int(eid[-1]) - 1
            id_array.append(eid)

    events.extend(load_events(sp.array(id_array), outdir, confidence, verbose))

def _parse_gene_info(track_info, genes, gene_names, gids, all_gene_names, outdir, confidence, validate_sg, verbose):
    for g in track_info:
        gid = sp.where(all_gene_names == g)[0]
        if gid.shape[0] == 0:
            sys.stderr.write('ERROR: provided gene ID "%s" could not be found, please check for correctness\n' % g)
            sys.exit(1)
        assert gid.shape[0] == 1
        gids.append(gid[0])
    genes.extend(load_genes(outdir, confidence, validate_sg, verbose, sp.array(gids)))
    gene_names.extend([_.name for _ in genes])

def _parse_test_info(test_info, outdir, confidence, data_tracks, test_tag):

    ### get correct testing directory
    if len(test_info) == 0 or test_info[0] == 'default':
        testdir = os.path.join(outdir, 'testing')
    else:
        testdir = os.path.join(outdir, test_info[0])
    assert os.path.exists(testdir), '\nERROR: directory containing test results (%s) could not be found.' % testdir

    ### get event types to be plotted
    if len(test_info) < 2:
        event_types = EVENT_TYPES 
    elif test_info[1] == 'any':
        event_types = EVENT_TYPES
    else:
        event_types = test_info[1].split(',')

    ### get number of top events to be plotted
    if len(test_info) < 3:
        topk = 1
    else:
        topk = int(test_info[2])

    for event_type in event_types:
        
        fname_pickle = os.path.join(testdir, 'test_setup_C%i_%s.pickle' % (confidence, event_type))
        fname_result = os.path.join(testdir, 'test_results_C%i_%s.tsv' % (confidence, event_type))
        if not os.path.exists(fname_pickle) or not os.path.exists(fname_result):
            if test_info[0] != 'any':
                sys.stderr.write('WARNING: no test information found for event type %s.\nPlease make sure that test results for this event type are available. If this is the case, try to re-run the test using the latest version of SplAdder.\n' % (event_type))
            continue

        ### load test setup
        test_setup = pickle.load(open(os.path.join(testdir, 'test_setup_C%i_%s.pickle' % (confidence, event_type)), 'rb'), encoding='latin1')
        halfsize = test_setup[3].shape[0] // 2
        idxA = sp.where(test_setup[3][:halfsize, 2] == 0)[0]
        idxB = sp.where(test_setup[3][:halfsize, 2] == 1)[0]
        ### iterate over top k test results
        for k, line in enumerate(open(os.path.join(testdir, 'test_results_C%i_%s.tsv' % (confidence, event_type)), 'r')):
            if k == 0:
                continue
            if k > topk:
                break
            sl = line.strip().split('\t')
            data_tracks.append([])
            data_tracks[-1].append(['segments', ','.join(test_setup[1][idxA]), ','.join(test_setup[1][idxB])])
            data_tracks[-1].append(['event', sl[0]])
            data_tracks[-1].append(['gene', sl[1]]) 
            test_tag.append('.difftest_%s' % sl[0])


def _set_plotting_range(genes, events, coords):
    plotrange = None
    for g in genes:
        if not plotrange:
            plotrange = [g.start, g.stop]
        else:
            plotrange[0] = min(g.start, plotrange[0])
            plotrange[1] = max(g.stop, plotrange[1])
    for e in events:
        if not plotrange:
            plotrange = [e.exons2.min(), e.exons2.max()]
        else:
            plotrange[0] = min(e.exons2.min(), plotrange[0])
            plotrange[1] = max(e.exons2.max(), plotrange[1])
    for c in coords:
        if not plotrange:
            plotrange = [c.start, c.stop]
        else:
            plotrange[0] = min(c.start, plotrange[0])
            plotrange[1] = max(c.stop, plotrange[1])

    if plotrange[1] - plotrange[0] > 1000000:
        sys.stderr.write('ERROR: plotting range has a width of more than 1 000 000 positions (%i) - aborting\n' % (plotrange[1] - plotrange[0]))
        sys.exit(1)

    ### slightly increase plotrange to get margins to left and right
    plotrange_orig = copy.copy(plotrange)
    delta = max(10, int((plotrange[1] - plotrange[0]) * 0.05))
    plotrange[0] -= delta
    plotrange[1] += delta

    return plotrange, plotrange_orig


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

    ### load gene information
    all_gene_names = get_gene_names(options.outdir, options.confidence, options.validate_sg, options.verbose)

    ### set color maps
    cmap_cov = plt.get_cmap('jet')
    cmap_edg = plt.get_cmap('jet')

    ### Data Range
    RangeData = namedtuple('RangeData', ['chr', 'start', 'stop'])

    plots = []
    test_tag = []
    ### if the --test option was given, populate the options object with additional track fields
    ### otherwise just create a single plot
    if options.test is None:
        plots.append(options.data_tracks)
        test_tag.append('')
    else:
        for test_info in options.test:
            _parse_test_info(test_info, options.outdir, options.confidence, plots, test_tag)
        for plot_data in plots:
            if len(options.data_tracks) > 0:
                plot_data.append(options.data_tracks)

    ### generate all plots to be completed
    for p, plot_data in enumerate(plots):
    
        ### get range information
        genes, gene_names, gids = [], [], []
        events, eids = [], []
        coords = []
        for range_info in options.range:
            ### genes
            if range_info[0] == 'gene':
                _parse_gene_info(range_info[1:], genes, gene_names, gids, all_gene_names, options.outdir, options.confidence, options.validate_sg, options.verbose)
            ### events
            elif range_info[0] == 'event':
                _parse_event_info(eids, gids, range_info[1:], events, options.outdir, options.confidence, options.verbose)
            ### coordinate ranges
            elif range_info[0] == 'coordinate':
                coords.append(RangeData._make(range_info[1:4]))

        ### check data tracks for additional range information if no --range was given
        if len(genes) + len(events) + len(coords) == 0:
            for range_info in plot_data:
                ### splicegraph / transcripts
                if range_info[0] in ['splicegraph', 'transcript']:
                    _parse_gene_info(range_info[1:], genes, gene_names, gids, all_gene_names, options.outdir, options.confidence, options.validate_sg, options.verbose)
                ### events
                elif range_info[0] == 'event':
                    _parse_event_info(eids, gids, range_info[1:], events, options.outdir, options.confidence, options.verbose)
                ### gene (this is only needed internally for --test mode)
                elif range_info[0] == 'gene':
                    _parse_gene_info(range_info[1:], genes, gene_names, gids, all_gene_names, options.outdir, options.confidence, options.validate_sg, options.verbose)
                    

        ### check that everything is on the same chromosome
        plotchrm = sp.unique([_.chr for _ in sp.r_[events, genes, coords]])
        if plotchrm.shape[0] > 1:
            sys.stderr.write('ERROR: the provided gene/event/coordinate ranges are on different chromosomes and canot be plotted jointly\n')
            sys.exit(1)

        ### identify the plotting range
        plotrange, plotrange_orig = _set_plotting_range(genes, events, coords)

        ### parse all elements to be plotted as data tracks
        data_tracks = []
        for data_element in plot_data:
            track_types = data_element[0].split(',')
            for track_type in track_types:
                data_tracks.append(dt.DataTrack(track_type, data_element[1:]))

        ### set up figure and layout
        gs = gridspec.GridSpec(len(data_tracks), 1, hspace=0.05)
        fig = plt.figure(figsize = (15, 3 * len(data_tracks)), dpi=200)

        axes = []
        for i,data_track in enumerate(data_tracks):

            ### plot splicing graph
            if data_track.type == 'splicegraph':
                ax = _add_ax(axes, fig, gs)
                if len(data_track.event_info) == 0:
                    _genes = genes
                    _gene_names = gene_names
                else:
                    _genes, _gene_names, _gids = [], [], []
                    _parse_gene_info(range_info[1:], _genes, _gene_names, _gids, all_gene_names, options.outdir, options.confidence, options.validate_sg, options.verbose)
                for gene in _genes:
                    gene.from_sparse()
                    plot_graph(gene.splicegraph.vertices, gene.splicegraph.edges, ax, xlim=plotrange, label=gene.name)
                    gene.to_sparse()
                ax.set_ylabel('splicing graph')
                ax.get_yaxis().set_label_coords(1.02,0.5)
            ### plot annotated transcripts
            if data_track.type == 'transcript':
                ax = _add_ax(axes, fig, gs)
                if len(data_track.event_info) == 0:
                    _genes = genes
                    _gene_names = gene_names
                else:
                    _genes, _gene_names, _gids = [], [], []
                    _parse_gene_info(range_info[1:], _genes, _gene_names, _gids, all_gene_names, options.outdir, options.confidence, options.validate_sg, options.verbose)
                for gene in _genes:
                    gene.from_sparse()
                    multiple(gene.exons, ax=ax, x_range=plotrange, labels=gene.transcripts, grid=True)
                    gene.to_sparse()
                ax.set_ylabel('transcripts')
                ax.get_yaxis().set_label_coords(1.02,0.5)
            ### plot events
            if data_track.type == 'event':
                ax = _add_ax(axes, fig, gs)
                ### no event ids given - plot the one from range
                if len(data_track.event_info) == 0:
                    _events = events
                ### events are given in the track - plot those instead
                else:
                    _eids, _events = [], []
                    _parse_event_info(_eids, gids, data_track.event_info, _events, options.outdir, options.confidence, options.verbose)
                event_list = [[event.exons1, event.exons2] for event in _events]
                labels = ['_'.join([_.event_type, str(_.id)]) for _ in _events]
                multiple(event_list, ax=ax, x_range=plotrange, color='green', padding=None, grid=True, labels=labels) 
                vax.clean_axis(ax, allx=True)
                ax.set_ylabel('events')
                ax.get_yaxis().set_label_coords(1.02,0.5)
            ### plot coverage tracks
            if data_track.type == 'coverage':
                ax = _add_ax(axes, fig, gs)
                min_sample_size = 5 ### TODO make that an option
                caxes = []
                for gene in genes:
                    caxes = []
                    labels = []
                    norm = plt.Normalize(0, len(data_track.bam_fnames))
                    for g, bam_group in enumerate(data_track.bam_fnames):
                        label = 'group %i' % (g + 1) if len(data_track.group_labels) == 0 else data_track.group_labels[g]
                        caxes.append(cov_from_bam(gene.chr, 
                                                  plotrange[0], 
                                                  plotrange[1], 
                                                  bam_group, 
                                                  subsample=min_sample_size, 
                                                  ax=ax, 
                                                  intron_cnt=True,
                                                  log=options.log, 
                                                  xlim=plotrange, 
                                                  color_cov=cmap_cov(norm(g)), 
                                                  color_intron_edge= cmap_edg(norm(g)),
                                                  grid=True, 
                                                  min_intron_cnt=options.mincount, 
                                                  return_legend_handle=True, 
                                                  label=label))
                        labels.append(label)
                ax.get_yaxis().set_label_coords(1.02,0.5)
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
                ax.get_yaxis().set_label_coords(1.02,0.5)
                ax.set_ylabel('segment counts')

        ### set x axis ticks
        for i, ax in enumerate(axes):
            if i == len(axes) - 1:
                ax.set_xticks(sp.around(sp.linspace(plotrange_orig[0], plotrange_orig[1], 10)))
                ax.set_xlabel('chromosome ' + plotchrm[0])
            else:
                ax.set_xticks([])

        ### save plot into file
        if options.format == 'd3':
            out_fname = os.path.join(dirname, 'plots', '%s%s.html' % (options.outbase, test_tag[p]))
            plugins.clear(fig)
            plugins.connect(fig, plugins.Zoom(enabled=True))
            mpld3.save_html(fig, open(out_fname, 'w'))
        else:
            out_fname = os.path.join(dirname, 'plots', '%s%s.%s' % (options.outbase, test_tag[p], options.format))
            plt.savefig(out_fname, format=options.format, bbox_inches='tight')
        plt.close(fig)

