import sys
import os
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import h5py
import cPickle
import pdb
import copy
import scipy.io as scio

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
    optional.add_option('-b', '--bams', dest='bams', metavar='FILE1A,FILE2A:FILE1B,FILE2B,,...', help='alignment files in BAM format (comma separated list,colon separated groups)', default='-')
    optional.add_option('-L', '--labels', dest='labels', metavar='LABEL_A,LABEL_B,...', help='group labels for alignment files groups (comma separated list)', default='-')
    optional.add_option('-c', '--confidence', dest='confidence', metavar='INT', type='int', help='confidence level (0 lowest to 3 highest) [3]', default=3)
    optional.add_option('-m', '--mincount', dest='mincount', metavar='INT', type='int', help='minimum count of introns to be displayed in coverage plot [0]', default=0)
    optional.add_option('-l', '--log', dest='log', action='store_true', help='plot coverage information in log scale [off]', default=False)
    optional.add_option('-g', '--gene_name', dest='gene_name', metavar='STR', help='gene_name to be plotted', default=None)
    optional.add_option('-e', '--event_id', dest='event_id', metavar='STR', help='event to be plotted', default=None)
    optional.add_option('-f', '--format', dest='format', metavar='STR', help='plot file format [pdf]', default='pdf')
    optional.add_option('-V', '--validate_sg', dest='validate_sg', action='store_true', help='use validated splice graph [off]', default=False)
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
    if options.event_id is not None:
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
    if options.validate_sg:
        (genes, events) = cPickle.load(open(os.path.join(options.outdir, 'spladder', 'genes_graph_conf%s.merge_graphs.validated.pickle' % options.confidence), 'r'))
    else:
        (genes, events) = cPickle.load(open(os.path.join(options.outdir, 'spladder', 'genes_graph_conf%s.merge_graphs.pickle' % options.confidence), 'r'))

    rows = get_plot_len(options)
    fig = plt.figure(figsize = (18, 3*rows), dpi=200)
    gs = gridspec.GridSpec(rows, 1)
    axes = []
    xlim = None

    ### get coloring
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
        assert len(options.labels) == len(options.bams.strip(':').split(':')), "The number of given labels (%i) needs to match the number of given bam file groups (%i)" % (len(options.labels), len(options.bams.strip(':').split(':')))

    ### plot splicing graph
    if options.gene_name is not None:
        axes.append(fig.add_subplot(gs[len(axes), 0]))
        gid = sp.where(sp.array([x.name for x in genes]) == options.gene_name)[0]
        min_sample_size = min(20, min([len(x.split(',')) for x in options.bams.strip(':').split(':')]))
        if gid.shape[0] > 0:
            gene = copy.deepcopy(genes[gid[0]])
            del gid
            del genes
            plot_graph(gene.splicegraph.vertices, gene.splicegraph.edges, axes[-1])
            xlim = axes[-1].get_xlim()
            start = gene.splicegraph.vertices.min()
            stop = gene.splicegraph.vertices.max()
            axes[-1].set_title('Splicing graph for %s' % options.gene_name)

            if options.transcripts:
                ### plot annotated transcripts
                axes.append(fig.add_subplot(gs[len(axes), 0]))

                multiple(gene.exons, ax=axes[-1], x_range=xlim)                                                                                                                                                                                                                                                                                                                             
                axes[-1].set_title('Annotated Transcripts')

            ### plot coverage information for a set of samples
            if options.bams != '-':
                samples = options.bams.strip(':').split(':')
                for s, sample in enumerate(samples):
                    bams = sample.split(',')
                    axes.append(fig.add_subplot(gs[len(axes), 0]))
                    if options.labels != '-':
                        title = 'Expression (%s)' % options.labels[s]
                    else:
                        title = 'Expression (sample %i)' % (s + 1)
                    cov_from_bam(gene.chr, start, stop, bams, subsample=min_sample_size, ax=axes[-1], intron_cnt=True, log=options.log, title=title, xlim=xlim, color_cov='#d7191c', color_intron_edge='#1a9641', grid=True, min_intron_cnt=options.mincount)
                    xlim = axes[-1].get_xlim()
                    axes[-1].set_xlabel('')
                
                ### plot all the samples in a single plot
                if len(samples) > 1:
                    axes.append(fig.add_subplot(gs[len(axes), 0]))
                    norm = plt.Normalize(0, len(samples))
                    caxes = []
                    labels = []
                    for s, sample in enumerate(samples):
                        bams = sample.split(',')
                        if options.labels != '-':
                            caxes.append(cov_from_bam(gene.chr, start, stop, bams, subsample=min_sample_size, ax=axes[-1], intron_cnt=True, log=options.log, title='Expression all Samples', xlim=xlim, color_cov=cmap_cov(norm(s)), color_intron_edge=cmap_edg(norm(s)), grid=True, min_intron_cnt=options.mincount, return_legend_handle=True, label=options.labels[s]))
                            labels.append(options.labels[s])
                        else:
                            caxes.append(cov_from_bam(gene.chr, start, stop, bams, subsample=min_sample_size, ax=axes[-1], intron_cnt=True, log=options.log, title='Expression all Samples', xlim=xlim, color_cov=cmap_cov(norm(s)), color_intron_edge=cmap_edg(norm(s)), grid=True, min_intron_cnt=options.mincount, return_legend_handle=True, label='sample %i' % (s + 1)))
                            labels.append('sample %i' % (s + 1))

                    plt.legend(caxes, labels)
                    axes[-1].set_xlabel('')

            ### plot structure of a single given event
            if options.event_id is not None:
                axes.append(fig.add_subplot(gs[len(axes), 0]))
                event_info = [x[::-1] for x in re.split(r'[._]', options.event_id[::-1], maxsplit=1)[::-1]]
                event = scio.loadmat(os.path.join(options.outdir, 'merge_graphs_%s_C%s.mat' % (event_info[0], options.confidence)), struct_as_record=False)['events_all'][0, int(event_info[1]) - 1]
                if event.event_type[0] == 'exon_skip':
                    exons = [sp.r_[event.exon_pre, event.exon_aft], sp.r_[event.exon_pre, event.exon, event.exon_aft]]
                elif event.event_type[0] == 'intron_retention':
                    exons = [sp.r_[event.exon1, event.exon2], sp.array([event.exon1[0, 0], event.exon2[0, 1]])]
                elif event.event_type[0] in ['alt_3prime', 'alt_5prime']:
                    exons = [sp.r_[event.exon_const, event.exon_alt1], sp.r_[event.exon_const, event.exon_alt2]]
                    for e, ex in enumerate(exons):
                        s_idx = sp.argsort(ex[:, 0])
                        exons[e] = ex[s_idx, :]
                elif event.event_type[0] == 'mutex_exons':
                    exons = [sp.r_[event.exon_pre, event.exon1, event.exon_aft], sp.r_[event.exon_pre, event.exon2, event.exon_aft]]
                elif event.event_type[0] == 'mult_exon_skip':
                    exons = [sp.r_[event.exon_pre, event.exon_aft], sp.r_[event.exon_pre, event.exons.reshape(event.exons.shape[1] / 2, 2), event.exon_aft]]
                multiple(exons, ax=axes[-1], x_range=xlim, color='green') 
                axes[-1].set_title('Event structure of %s' % options.event_id)
                event_tag = '.%s' % options.event_id

        ### plot the identified events for this gene
        out_fname = os.path.join(options.outdir, 'plots', 'gene_overview_%s%s%s.%s' % (options.gene_name, event_tag, log_tag, options.format))
        plt.savefig(out_fname, format=options.format, bbox_inches='tight')
    plt.close(fig)




if __name__ == "__main__":
    spladder_viz()

