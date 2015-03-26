import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import h5py
import cPickle
import pdb

from modules.classes.gene import Gene
from modules.viz.graph import *
from modules.viz.coverage import *
from modules.viz.genelets import *

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'MANDATORY')
    required.add_option('-o', '--outdir', dest='outdir', metavar='DIR', help='spladder directory containing the spladder results', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-b', '--bams', dest='bams', metavar='FILE1,FILE2,...', help='alignment files in BAM format (comma separated list)', default='-')
    optional.add_option('-c', '--confidence', dest='confidence', metavar='INT', type='int', help='confidence level (0 lowest to 3 highest) [3]', default=3)
    optional.add_option('-l', '--log', dest='log', action='store_true', help='plot coverage information in log scale [off]', default=False)
    optional.add_option('-g', '--gene_name', dest='gene_name', metavar='STR', help='gene_name to be plotted', default=None)
    optional.add_option('-f', '--format', dest='format', metavar='STR', help='plot file format [pdf]', default='pdf')
    optional.add_option('-t', '--transcripts', dest='transcripts', action='store_true', help='plot annotated transcripts', default=False)
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default=False)
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()

    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)

    options.parser = parser
    return options

def get_plot_len(options):
    """Identifies the number of rows we need in our plot"""

    rows = 0
    if options.gene_name is not None:
        rows += 1
    if options.bams != '-':
        samples = options.bams.strip(':').split(':')
        rows += len(samples)
        if len(samples) > 1:
            rows += 1
    if options.transcripts:
        rows += 1

    return rows

def spladder_viz():

    """Main visualization code"""
    
    ### parse command line parameters
    options = parse_options(sys.argv)

    ### create plot directory if it does not exist yet
    if not os.path.exists(os.path.join(options.outdir, 'plots')):
        os.mkdir(os.path.join(options.outdir, 'plots'))

    ### load gene information
    (genes, events) = cPickle.load(open(os.path.join(options.outdir, 'spladder', 'genes_graph_conf%s.merge_graphs.pickle' % options.confidence), 'r'))

    rows = get_plot_len(options)
    fig = plt.figure(figsize = (18, 3*rows), dpi=200)
    gs = gridspec.GridSpec(rows, 1)
    axes = []
    xlim = None

    ### get coloring
    cmap_cov = plt.get_cmap('jet')
    cmap_edg = plt.get_cmap('jet')

    ### plot splicing graph
    if options.gene_name is not None:
        axes.append(fig.add_subplot(gs[len(axes), 0]))
        i = sp.where(sp.array([x.name for x in genes]) == options.gene_name)[0]
        if i.shape[0] > 0:
            plot_graph(genes[i[0]].splicegraph.vertices, genes[i[0]].splicegraph.edges, axes[-1])
            xlim = axes[-1].get_xlim()
            start = genes[i[0]].splicegraph.vertices.min()
            stop = genes[i[0]].splicegraph.vertices.max()
            axes[-1].set_title('Splicing graph for %s' % options.gene_name)

        if options.transcripts:
            ### plot annotated transcripts
            axes.append(fig.add_subplot(gs[len(axes), 0]))
            #multiple([x for x in anno['genes'][0, a_idx]['exons'][0, :]], ax=ax, x_range=xlim)                                                                                                                                                                                                                                                                                                                             
            multiple(genes[i[0]].exons, ax=axes[-1], x_range=xlim)                                                                                                                                                                                                                                                                                                                             
            axes[-1].set_title('Annotated Transcripts (TAIR 10)')

        ### plot coverage information for a set of samples
        if options.bams != '-':
            samples = options.bams.split(':')
            for s, sample in enumerate(samples):
                bams = sample.split(',')
                axes.append(fig.add_subplot(gs[len(axes), 0]))
                cov_from_bam(genes[i[0]].chr, start, stop, bams, subsample=20, ax=axes[-1], intron_cnt=True, log=options.log, title='Expression (Sample %i)' % (s+1), xlim=xlim, color_cov='#d7191c', color_intron_edge='#1a9641', grid=True)
                xlim = axes[-1].get_xlim()
                axes[-1].set_xlabel('')
            
            ### plot all the samples in a single plot
            if len(samples) > 1:
                axes.append(fig.add_subplot(gs[len(axes), 0]))
                norm = plt.Normalize(0, len(samples))
                for s, sample in enumerate(samples):
                    bams = sample.split(',')
                    cov_from_bam(genes[i[0]].chr, start, stop, bams, subsample=20, ax=axes[-1], intron_cnt=True, log=options.log, title='Expression all Samples', xlim=xlim, color_cov=cmap_cov(norm(s)), color_intron_edge=cmap_edg(norm(s)), grid=True)
                axes[-1].set_xlabel('')

        ### plot the identified events for this gene

        plt.savefig('test_%s.%s' % (options.gene_name, options.format), format=options.format, bbox_inches='tight')
    plt.close(fig)




if __name__ == "__main__":
    spladder_viz()

