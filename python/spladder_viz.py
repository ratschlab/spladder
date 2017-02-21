#! /usr/bin/env python 
import sys
import os
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import cPickle
import pdb

import modules.settings as settings
from modules.classes.gene import Gene
from modules.viz.graph import *
from modules.viz.coverage import *
from modules.viz.genelets import *
from modules.identity import *

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'MANDATORY')
    required.add_option('-o', '--outdir', dest='outdir', metavar='DIR', help='spladder directory containing the spladder results', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-b', '--bams', dest='bams', metavar='FILE1A,FILE2A:FILE1B,FILE2B,,...', help='alignment files in BAM format (comma separated list,colon separated groups)', default='-')
    optional.add_option('-L', '--labels', dest='labels', metavar='LABEL_A,LABEL_B,...', help='group labels for alignment files groups (comma separated list)', default='-')
    optional.add_option('-c', '--confidence', dest='confidence', metavar='INT', type='int', help='confidence level (0 lowest to 3 highest) [3]', default=3)
    optional.add_option('-m', '--mincount', dest='mincount', metavar='INT', type='int', help='minimum count of introns to be displayed in coverage plot [0]', default=0)
    optional.add_option('-l', '--log', dest='log', action='store_true', help='plot coverage information in log scale [off]', default=False)
    optional.add_option('-g', '--gene_name', dest='gene_name', metavar='STR', help='gene_name to be plotted', default=None)
    optional.add_option('-e', '--event_id', dest='event_id', metavar='STR', help='event to be plotted', default=None)
    optional.add_option('-f', '--format', dest='format', metavar='STR', help='plot file format [pdf, png, d3]', default='pdf')
    optional.add_option('', '--zoom_x', dest='zoom_x', metavar='percent_left,percent_right', help='zoom x axis from percent_left to percent_right [0.0,1.0]', default='0.0,1.0')
    optional.add_option('-V', '--validate_sg', dest='validate_sg', metavar='y|n', help='validate splice graph [n]', default='n')
    optional.add_option('-T', '--transcripts', dest='transcripts', metavar='y|n', help='plot annotated transcripts', default='n')
    optional.add_option('--test-result', dest='test_result', metavar='INT', type='int', help='plot top k significant events from test', default=0)
    optional.add_option('--test-labels', dest='test_labels', metavar='STR', type='str', help='labels used for the groups in the test (order matters) [condA:condB]', default='condA:condB')
    optional.add_option('-t', '--event_types', dest='event_types', metavar='EVENT1,EVENT2,...', help='list of alternative splicing events to extract [exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons]', default='exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons')
    optional.add_option('-v', '--verbose', dest='verbose', metavar='y|n', help='verbosity', default='n')
    optional.add_option('-d', '--debug', dest='debug', metavar='y|n', help='use debug mode [n]', default='n')
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()
    #options.event_types = options.event_types.strip(',').split(',')

    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)

    options.parser = parser
    return options

def get_plot_len(CFG):
    """Identifies the number of rows we need in our plot"""

    rows = 3 # splicing graph + events + segments
    if len(CFG['bam_fnames']) > 0:
        rows += len(CFG['bam_fnames'])
        if len(CFG['bam_fnames']) > 1:
            rows += 1
    rows += int(CFG['plot_transcripts'])

    return rows


def spladder_viz():

    """Main visualization code"""
    
    ### parse command line parameters
    options = parse_options(sys.argv)

    ### parse parameters from options object
    CFG = settings.parse_args(options, identity='viz')
   
    ### create plot directory if it does not exist yet
    if not os.path.exists(os.path.join(options.outdir, 'plots')):
        os.mkdir(os.path.join(options.outdir, 'plots'))

    if options.format == 'd3':
        try:
            import mpld3
            from mpld3 import plugins
        except ImportError:
            sys.stderr.write("ERROR: missing package for output format d3. Package mpld3 required")
            sys.exit(1)

    ### load gene information
    genes = load_genes(CFG)

    rows = get_plot_len(CFG)
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
    if options.labels != '-':
        options.labels = options.labels.strip(',').split(',')
        assert len(options.labels) == len(CFG['bam_fnames']), "The number of given labels (%i) needs to match the number of given bam file groups (%i)" % (len(options.labels), len(CFG['bam_fnames']))

    ### the user chose a specific gene for plotting
    if options.gene_name is not None:
        gid = sp.where(sp.array([x.name.split('.')[0] for x in genes]) == options.gene_name.split('.')[0])[0]
        if gid.shape[0] == 0:
            sys.stderr.write('ERROR: provided gene ID %s could not be found, please check for correctness\n' % options.gene_name)
            sys.exit(1)
        gids = [[sp.where(sp.array([x.name for x in genes]) == options.gene_name)[0][0], options.event_id]]
    ### the plotting happens on the results of spladder test
    ### the user chooses to plot the top k significant events
    ### this requires the event type to be specified
    elif options.test_result > 0:
        gene_names = []
        for event_type in CFG['event_types']:
            ### the testing script should generate a setup file for the test
            ### SETUP is structured as follows:
            ###  [gene_strains, event_strains, dmatrix0, dmatrix1, event_type, options, CFG]
            labels = options.test_labels.split(':')
            SETUP = cPickle.load(open(os.path.join(CFG['out_dirname'], 'testing_%s_vs_%s' % (labels[0], labels[1]), 'test_setup_C%i_%s.pickle' % (CFG['confidence_level'], event_type)), 'r'))

            ### get strains to plot
            idx1 = sp.where(sp.in1d(SETUP[0], SETUP[6]['conditionA']))[0]
            idx2 = sp.where(sp.in1d(SETUP[0], SETUP[6]['conditionB']))[0]
    
            ### load test results
            for l, line in enumerate(open(os.path.join(CFG['out_dirname'], 'testing_%s_vs_%s' % (labels[0], labels[1]), 'test_results_C%i_%s.tsv' % (CFG['confidence_level'], event_type)), 'r')):
                if l == 0:
                    continue
                if l > options.test_result:
                    break
                sl = line.strip().split('\t')
                gene_names.append([sl[1], sl[0]])
        gids = get_gene_ids(CFG, gene_names)
    ### no gene specified but result provided - plot all genes with confirmed events
    ### if an event_id is provided, only the associated gene will be plotted
    else:
        gids = get_gene_ids(CFG)
        
    ### iterate over genes to plot
    for gid in gids:
        if options.format == 'd3':
            fig = plt.figure(figsize = (12, 2*rows), dpi=100)
        else:
            fig = plt.figure(figsize = (18, 3*rows), dpi=200)
        axes = []

        ### gather information about the gene we plot
        gene = get_gene(genes[gid[0]])
        gene.from_sparse()
        if CFG['verbose']:
            print 'plotting information for gene %s' % gene.name

        ### plot splicing graph
        axes.append(fig.add_subplot(gs[len(axes), 0]))
        plot_graph(gene.splicegraph.vertices, gene.splicegraph.edges, axes[-1])
        xlim = axes[-1].get_xlim()
        axes[-1].set_title('Splicing graph for %s' % gene.name)

        ### plot annotated transcripts
        if CFG['plot_transcripts']:
            axes.append(fig.add_subplot(gs[len(axes), 0], sharex=axes[0]))
            multiple(gene.exons, ax=axes[-1], x_range=xlim)
            axes[-1].set_title('Annotated Transcripts')

        ### plot coverage information for a set of given samples
        if len(CFG['bam_fnames']) > 0:
            plot_bam(options, gene, CFG['bam_fnames'], fig, axes, gs, xlim, cmap_cov, cmap_edg)
           
            ### plot all the samples in a single plot
            if len(CFG['bam_fnames']) > 1:
                plot_bam(options, gene, CFG['bam_fnames'], fig, axes, gs, xlim, cmap_cov, cmap_edg, single=False)

        ### plot segment counts
        print 'get segment counts'
        (segments, edges, edge_idx, strains) = get_seg_counts(CFG, gid[0])
        seg_sample_idx = None
        if len(CFG['strains']) > 0:
            seg_sample_idx = []
            for group in CFG['strains']:
                seg_sample_idx.append(sp.where(sp.in1d(strains, group))[0])
        if options.test_result > 0:
            seg_sample_idx = [idx1, idx2]
        axes.append(fig.add_subplot(gs[len(axes), 0], sharex=axes[0]))
        print 'plot segment counts'
        if identity() == 'matlab':
            cov_from_segments(gene, segments, edges, edge_idx, axes[-1], xlim=xlim, log=options.log, grid=True, order='F')
        else:
            cov_from_segments(gene, segments, edges, edge_idx, axes[-1], xlim=xlim, log=options.log, grid=True, order='C', sample_idx=seg_sample_idx)
        axes[-1].set_title('Segment counts')


        ### plot structure of a single given event
        #if options.event_id is not None:
        axes.append(fig.add_subplot(gs[len(axes), 0], sharex=axes[0]))

        ### event to plot is specified with the gene id list
        if gid[1] is not None:
            event_info = [x[::-1] for x in re.split(r'[._]', gid[1][::-1], maxsplit=1)[::-1]]
            event_info[1] = int(event_info[1]) - 1
            event_info = sp.array(event_info, dtype='str')[sp.newaxis, :]
            event_tag = '.%s' % gid[1]
        #if options.event_id is not None:
        #    event_info = [x[::-1] for x in re.split(r'[._]', options.event_id[::-1], maxsplit=1)[::-1]]
        #    event_info[1] = int(event_info[1]) - 1
        #    event_info = sp.array(event_info, dtype='str')[sp.newaxis, :]
        #    event_tag = '.%s' % options.event_id
        ### get all significant events of the current gene
        else:
            event_info = get_conf_events(CFG, gid[0])
        
        plot_event(CFG, event_info, axes[-1], xlim)

        ### we only need to adapt the xoom for one axis object - as we share the x
        zoom_x = [float(x) for x in options.zoom_x.split(',')]
        xlim = axes[0].get_xlim()
        xdiff = xlim[1] - xlim[0]
        axes[0].set_xlim([xlim[0] + (zoom_x[0] * xdiff), xlim[0] + (zoom_x[1] * xdiff)])

        plt.tight_layout()
        ### save plot into file
        if options.format == 'd3':
            out_fname = os.path.join(options.outdir, 'plots', 'gene_overview_C%i_%s%s%s.html' % (options.confidence, gene.name, event_tag, log_tag))
            plugins.clear(fig)
            plugins.connect(fig, plugins.Zoom(enabled=True))
            mpld3.save_html(fig, open(out_fname, 'w'))
        else:
            if options.test_result > 0:
                out_fname = os.path.join(options.outdir, 'plots', 'gene_overview_C%i_%s%s%s.%s' % (options.confidence, gene.name, event_tag, log_tag, options.format))
            else:
                out_fname = os.path.join(options.outdir, 'plots', 'gene_overview_C%i_%s%s%s.%s' % (options.confidence, gene.name, event_tag, log_tag, options.format))
            plt.savefig(out_fname, format=options.format, bbox_inches='tight')
        plt.close(fig)


def plot_bam(options, gene, samples, fig, axes, gs, xlim, cmap_cov, cmap_edg, single=True):

    min_sample_size = min(20, min([len(x) for x in samples]))
    start = gene.splicegraph.vertices.min()
    stop = gene.splicegraph.vertices.max()

    norm = plt.Normalize(0, len(samples))
    caxes = []
    labels = []

    for s, bams in enumerate(samples):

        if options.labels != '-':
            label = options.labels[s]
        else:
            label = 'group %i' % (s + 1)
            
        if single:
            axes.append(fig.add_subplot(gs[len(axes), 0], sharex=axes[0]))
            title = 'Expression (%s)' % label
            color_cov = cmap_cov(norm(0)) # '#d7191c'
            color_edg = cmap_edg(norm(0)) # '#1a9641'
        else:
            if s == 0:
                axes.append(fig.add_subplot(gs[len(axes), 0], sharex=axes[0]))
            title = 'Expression all Sample Groups'
            color_cov = cmap_cov(norm(s))
            color_edg = cmap_edg(norm(s))

        caxes.append(cov_from_bam(gene.chr, start, stop, bams, subsample=min_sample_size, ax=axes[-1], intron_cnt=True, 
                     log=options.log, title=title, xlim=xlim, color_cov=color_cov, color_intron_edge=color_edg,
                     grid=True, min_intron_cnt=options.mincount, return_legend_handle=True, label=label))
        labels.append(label)
        axes[-1].set_xlabel('')

    if not single:
        plt.legend(caxes, labels)
 

def plot_event(CFG, event_info, ax, xlim):
    """This function takes the event_id given in the CFG object and 
    plots it into ax."""

    ### plot it
    event_list = [_ for event in load_events(CFG, event_info) for _ in [event.exons1, event.exons2]]
    multiple(event_list, ax=ax, x_range=xlim, color='green') 
    ax.set_title('Alt event structure') # of %s' % options.event_id)

if __name__ == "__main__":
    spladder_viz()

