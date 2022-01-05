import sys

### check which python version we are running
if sys.version_info[0] < 3:
    sys.stderr.write('\nERROR: SplAdder requires Python 3 to run. You are currently running Python %i.%i.%i\n' % (sys.version_info[0], sys.version_info[1], sys.version_info[2]))
    sys.exit(1)

from .spladder_build import spladder
from .spladder_viz import spladder_viz
from .spladder_test import spladder_test
from .spladder_prep import spladder_prep

def parse_options(argv):

    """Parses options from the command line """

    from argparse import ArgumentParser

    parser = ArgumentParser(prog='spladder')
    subparsers = parser.add_subparsers(help='Running modes', metavar='{prep, build, test, viz}')

    ### RUN MODE "BUILD"
    parser_build = subparsers.add_parser('build', help='run mode to build and construct graphs')
    required = parser_build.add_argument_group('MANDATORY')
    required.add_argument('-b', '--bams', dest='bams', metavar='FILE1,FILE2,...', help='alignment files in BAM format (comma separated list)', default='-', required=True)
    required.add_argument('-o', '--outdir', dest='outdir', metavar='DIR', help='output directory', default='-', required=True)
    required.add_argument('-a', '--annotation', dest='annotation', metavar='FILE', help='file name for annotation in GTF/GFF3 or format', default='-', required=True)

    general = parser_build.add_argument_group('GENERAL')
    general.add_argument('--parallel', dest='parallel', metavar='<INT>', type=int, help='use INT processors [1]', default=1)
    general.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='use verbose output mode [off]', default=False)
    general.add_argument('-d', '--debug', dest='debug', action='store_true', help='use debug output mode [off]', default=False)

    inputs = parser_build.add_argument_group('INPUT OPTIONS')
    inputs.add_argument('-n', '--readlen', dest='readlen', metavar='INT', type=int, help='read length (used for automatic confidence levele settings) [36]', default=50)
    inputs.add_argument('--primary-only', dest='primary_only', action='store_true', help='only use primary alignments [on]', default=True)
    inputs.add_argument('--no-primary-only', dest='primary_only', action='store_false', default=True)
    inputs.add_argument('--var-aware', dest='var_aware', action='store_true', help='alignment files are variation aware (presence of XM and XG tags) [off]', default=False)
    inputs.add_argument('--no-var-aware', dest='var_aware', action='store_false', default=False)
    inputs.add_argument('--set-mm-tag', dest='mm_tag', help='sets the sequence of the mismatch tag used in alignments [NM]', default='NM')
    inputs.add_argument('--labels', dest='labels', metavar='STRING', help='use labels instead of bam file names (comma separated list) [-]', default='-')
    inputs.add_argument('--filter-overlap-genes', dest='filter_overlap_genes', action='store_true', help='remove genes from annotation that overlap each other [off]', default=False)
    inputs.add_argument('--filter-overlap-exons', dest='filter_overlap_exons', action='store_true', help='remove exons from annotation that overlap each other [off]', default=False)
    inputs.add_argument('--filter-overlap-transcripts', dest='filter_overlap_transcripts', action='store_true', help='remove transcripts from annotation that overlap each other [off]', default=False)
    inputs.add_argument('--filter-consensus', dest='filter_consensus', metavar='STRING', help='require new junctions to have consensus (needs ref genome) [off]; choices: strict (GT/AG), lenient (G[TC]/AG)', default='')
    inputs.add_argument('--ignore-mismatches', dest='ignore_mismatches', action='store_true', help='does not filter by edit operations - does not require NM in BAM [off]', default=False)
    inputs.add_argument('--reference', dest='ref_genome', metavar='FILE', help='reference genome (only needed for CRAM file de-compression or consensus filtering)', default='')

    outputs = parser_build.add_argument_group('OUTPUT OPTIONS')
    outputs.add_argument('-l', '--logfile', dest='logfile', metavar='FILE', help='log file name [stdout]', default='-')
    outputs.add_argument('--output-txt', dest='output_txt', action='store_true', help='outputs all events in txt format (can be big) [off]', default=False)
    outputs.add_argument('--output-txt-conf', dest='output_confirmed_txt', action='store_true', help='outputs confirmed events in txt format [on]', default=True)
    outputs.add_argument('--no-output-txt-conf', dest='output_confirmed_txt', action='store_false', default=True)
    outputs.add_argument('--output-gff3', dest='output_gff3', action='store_true', help='outputs all events in GFF3 format [off]', default=False)
    outputs.add_argument('--output-gff3-conf', dest='output_confirmed_gff3', action='store_true', help='outputs confirmed events in GFF3 format [on]', default=True)
    outputs.add_argument('--no-output-gff3-conf', dest='output_confirmed_gff3', action='store_false', default=True)
    outputs.add_argument('--output-struc', dest='output_struc', action='store_true', help='outputs all events in structured splicing syntax similar to astalavista [off]', default=False)
    outputs.add_argument('--output-struc-conf', dest='output_confirmed_struc', action='store_true', help='outputs confirmed events in structured splicing syntax similar to astalavista [off]', default=False)
    outputs.add_argument('--output-bed', dest='output_bed', action='store_true', help='output all events in BED format [off]', default=False)
    outputs.add_argument('--output-conf-bed', dest='output_confirmed_bed', action='store_true', help='output confirmed events in BED format [off]', default=False)
    outputs.add_argument('--output-conf-tcga', dest='output_confirmed_tcga', action='store_true', help='output confirmed events in format used for TCGA [off]', default=False)
    outputs.add_argument('--output-conf-icgc', dest='output_confirmed_icgc', action='store_true', help='output confirmed events in format used for ICGC [off]', default=False)
    outputs.add_argument('--sparse-bam', dest='sparse_bam', action='store_true', help='store BAM content as sparse representation for later use [off]', default=False)
    outputs.add_argument('--compress-text', dest='compress_text', action='store_true', help='compress text output [on]', default=True)
    outputs.add_argument('--no-compress-text', dest='compress_text', action='store_false', default=True)
    outputs.add_argument('--tmp-dir', dest='tmpdir', metavar='DIR', help='directory to store temporary data [<outdir>/tmp]', default='')

    graph = parser_build.add_argument_group('GRAPH GENERATION')
    graph.add_argument('-c', '--confidence', dest='confidence', metavar='INT', type=int, help='confidence level (0 lowest to 3 highest) [3]', default=3)
    graph.add_argument('-I', '--iterations', dest='insert_intron_iterations', metavar='INT', type=int, help='number of iterations to insert new introns into the graph [5]', default=5)
    graph.add_argument('-M', '--merge-strat', dest='merge', metavar='<STRAT>', help='merge strategy, where <STRAT> is one of: single, merge_bams, merge_graphs, merge_all [merge_graphs]', default='merge_graphs')
    graph.add_argument('--chunked-merge', dest='chunked_merge', metavar="LEVEL MAX_LEVEL START END", nargs='+', action='append', help='provide info for external merge with START being 0-based and END non-inclusive', default=[])
    graph.add_argument('--chunksize', dest='chunksize', metavar='INT', type=int, help='chunksize for chunked merge [10]', default=10)
    graph.add_argument('--insert-ir', dest='insert_ir', action='store_true', help='insert intron retentions [on]', default=True)
    graph.add_argument('--no-insert-ir', dest='insert_ir', action='store_false', default=True)
    graph.add_argument('--insert-es', dest='insert_es', action='store_true', help='insert cassette exons [on]', default=True)
    graph.add_argument('--no-insert-es', dest='insert_es', action='store_false', default=True)
    graph.add_argument('--insert-ni', dest='insert_ni', action='store_true', help='insert new intron edges [on]', default=True)
    graph.add_argument('--no-insert-ni', dest='insert_ni', action='store_false', default=True)
    graph.add_argument('--remove-se', dest='remove_se', action='store_true', help='remove short exons [off]', default=False)
    graph.add_argument('--no-remove-se', dest='remove_se', action='store_false', default=False)
    graph.add_argument('--validate-sg',  dest='validate_sg', action='store_true', help='validate splice graph [off]', default=False)
    graph.add_argument('--no-validate-sg', dest='validate_sg', action='store_false', default=False)
    graph.add_argument('--re-infer-sg', dest='infer_sg', action='store_true', help='re-infer splice graph [off] (advanced)', default=False)
    graph.add_argument('--no-re-infer-sg', dest='infer_sg', action='store_false', default=False)
    graph.add_argument('--validate-sg-count', dest='sg_min_edge_count', metavar='INT', type=int, help='number of samples supporting an edge for it to be kept [min(10, #samples)]', default=10)

    splice = parser_build.add_argument_group('AS EVENT EXTRACTION')
    splice.add_argument('--event-types', dest='event_types', metavar='STRING', help='list of alternative splicing events to extract [exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons]', default='exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons')
    splice.add_argument('--extract-ase', dest='extract_as', action='store_true', help='extract alternative splicing events [on]', default=True)
    splice.add_argument('--no-extract-ase', dest='extract_as', action='store_false', default=None)
    splice.add_argument('--ase-edge-limit', dest='detect_edge_limit', metavar='INT', type=int, help='max number of edges in the graph to still extract events for a gene [500]', default=500)
    splice.add_argument('--curate-alt-prime', dest='curate_alt_prime', action='store_true', help='curate alt prime events [on]', default=True)
    splice.add_argument('--no-curate-alt-prime', dest='curate_alt_prime', action='store_false', default=None)
    splice.add_argument('--quantify-graph', dest='quantify_graph', action='store_true', help='quantify graph [on]', default=True)
    splice.add_argument('--no-quantify-graph', dest='quantify_graph', action='store_false', default=None)
    splice.add_argument('--use-anno-support', dest='use_anno_support', action='store_true', help='use annotation for validating event introns [off]', default=False)
    splice.add_argument('--no-use-anno-support', dest='use_anno_support', action='store_false', default=None)
    splice.add_argument('--psi-min-reads', dest='psi_min_reads', metavar='INT', help='minimum number of spliced reads covering either isoform to compute PSI [10]', default=10)
    splice.add_argument('--qmode', dest='qmode', metavar='STRING', help='quantification mode: single, collect, all [all]', default='all')
    parser_build.set_defaults(func=spladder)

    ### RUN MODE "TEST"
    parser_test = subparsers.add_parser('test', help='run mode to differentially test events')
    required_test = parser_test.add_argument_group('MANDATORY')
    required_test.add_argument('-o', '--outdir', dest='outdir', metavar='DIR', help='spladder output directory', default='-')
    required_test.add_argument('-a', '--conditionA', dest='conditionA', metavar='idA1,idA2,idA3,...', help='comma separated list of alignment files for condition A', default='-')
    required_test.add_argument('-b', '--conditionB', dest='conditionB', metavar='idB1,idB2,idB3,...', help='comma separated list of alignment files for condition B', default='-')
    inputs_test = parser_test.add_argument_group('INPUT OPTIONS')
    inputs_test.add_argument('-n', '--readlen', dest='readlen', metavar='INT', type=int, help='read length [50]', default=50)
    inputs_test.add_argument('-c', '--confidence', dest='confidence', metavar='INT', type=int, help='confidence level (0 lowest to 3 highest) [3]', default=3)
    inputs_test.add_argument('-M', '--merge-strat', dest='merge', metavar='<STRAT>', help='merge strategy, where <STRAT> is one of: merge_bams, merge_graphs, merge_all [merge_graphs]', default='merge_graphs')
    inputs_test.add_argument('--event-types', dest='event_types', metavar='STRING', help='list of alternative splicing events to be tested [exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons]', default='exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons')
    inputs_test.add_argument('--validate-sg',  dest='validate_sg', action='store_true', help='splice graph is validated [off]', default=False)
    inputs_test.add_argument('--no-validate-sg', dest='validate_sg', action='store_false', default=False)
    testing = parser_test.add_argument_group('TESTING OPTIONS')
    testing.add_argument('-C', '--correction', dest='correction', metavar='STR', help='method for multiple testing correction (BH, Bonferroni, Holm, Hochberg, BY, TSBH) [BH]', default='BH')
    testing.add_argument('-0', '--max-zero-frac', dest='max_0_frac', metavar='FLOAT', type=float, help='max fraction of 0 values per event isoform quantification over all tested samples [0.5]', default=0.5)
    testing.add_argument('--dpsi', dest='min_dpsi', metavar='FLOAT', type=float, help='Delta PSI cutoff between tested groups for events to be considered for testing [0.05]', default=0.05)
    testing.add_argument('-i', '--min-count', dest='min_count', metavar='INT', type=int, help='min read count sum over all samples for an event isoform to be tested [10]', default=10)
    testing.add_argument('--cap-outliers', dest='cap_outliers', action='store_true', help='replace splice outliers with a max value [off]', default=False)
    testing.add_argument('--no-cap-outliers', dest='cap_outliers', action='store_false', default=False)
    testing.add_argument('--cap-exp-outliers', dest='cap_exp_outliers', action='store_true', help='replace expression outliers with a max value [on]', default=True)
    testing.add_argument('--no-cap-exp-outliers', dest='cap_exp_outliers', action='store_false', default=True)
    outputs_test = parser_test.add_argument_group('OUTPUT OPTIONS')
    outputs_test.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='user verbose output mode [off]', default=False)
    outputs_test.add_argument('-d', '--debug', dest='debug', action='store_true', help='use debug mode [off]', default=False)
    outputs_test.add_argument('-D', '--diagnose-plots', dest='diagnose_plots', action='store_true', help='generate diagnose plots [off]', default=False)
    outputs_test.add_argument('--timestamp', dest='timestamp', action='store_true', help='add timestamp to output directory [off]', default=False)
    outputs_test.add_argument('--labelA', dest='labelA', metavar='STRING', help='label for condition A (used for output naming)', default='condA')
    outputs_test.add_argument('--labelB', dest='labelB', metavar='STRING', help='label for condition B (used for output naming)', default='condB')
    outputs_test.add_argument('--out-tag', dest='out_tag', metavar='STRING', help='additional tag to label out directory', default='-')
    outputs_test.add_argument('-f', '--plot-format', dest='plot_format', metavar='STRING', help='format for diagnose/output plots [png]', default='png')
    experimental_test = parser_test.add_argument_group('EXPERIMENTAL - BETA STATE')
    experimental_test.add_argument('--parallel', dest='parallel', metavar='<INT>', type=int, help='use multiple processors [1]', default=1)
    experimental_test.add_argument('--non-alt-norm', dest='non_alt_norm', action='store_true', help='only use non alternative exon segments for gene expression counting [off]', default=False)
    experimental_test.add_argument('--high-memory', dest='high_memory', action='store_true', help='use more memory to decrease running time (mem-map counts) [off]', default=False)
    parser_test.set_defaults(func=spladder_test)

    ### RUN MODE "VIZ"
    parser_viz = subparsers.add_parser('viz', help='run mode to visualize graphs and events')
    required_viz = parser_viz.add_argument_group('MANDATORY')
    required_viz.add_argument('-o', '--outdir', dest='outdir', metavar='DIR', help='spladder directory containing the spladder results', required=True)
    required_viz.add_argument('--range', dest='range', metavar="TYPE SPECS", nargs='+', action='append', help='defines which genomic range should be plotted', default=[])
    required_viz.add_argument('--track', dest='data_tracks', metavar="TYPE [SAMPLES [SAMPLES]]", nargs='+', action='append', help='defines which type of plot should be generated on which samples', default=[])

    output_viz = parser_viz.add_argument_group('OUTPUT')
    output_viz.add_argument('-O', '--outbase', dest='outbase', metavar='NAME', help='name of the plot output file, excluding extension [gene_overview]. Will be placed in <outdir>/plots/', default='gene_overview')
    output_viz.add_argument('-m', '--mincount', dest='mincount', metavar='INT', type=int, help='minimum count of introns to be displayed in coverage plot [0]', default=0)
    output_viz.add_argument('-f', '--format', dest='format', metavar='STR', help='plot file format [pdf, png, d3]', default='pdf')
    output_viz.add_argument('-l', '--log', dest='log', action='store_true', help='plot coverage information in log scale [off]', default=False)

    general_viz = parser_viz.add_argument_group('GENERAL')
    general_viz.add_argument('-c', '--confidence', dest='confidence', metavar='INT', type=int, help='confidence level (0 lowest to 3 highest) [3]', default=3)
    general_viz.add_argument('-V', '--validate-sg', dest='validate_sg', action='store_true', help='validate splice graph [off]', default=False)
    general_viz.add_argument('-M', '--merge-strat', dest='merge', metavar='<STRAT>', help='merge strategy, where <STRAT> is one of: merge_bams, merge_graphs, merge_all [merge_graphs]', default='merge_graphs')
    general_viz.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='user verbose output mode [off]', default=False)
    general_viz.add_argument('-d', '--debug', dest='debug', action='store_true', help='use debug mode [off]', default=False)

    optional_viz = parser_viz.add_argument_group('EXPERIMENTAL - BETA STATE')
    optional_viz.add_argument('--test', dest='test', metavar='[GROUP EVENT_TYPE TOP_K]', nargs='+', action='append', help='plot results for differential test (optionally provide test name, event_type(s) and top k cutoff')
    optional_viz.add_argument('--testdir', dest='testdir', metavar='DIR', help='specify here in case testing results were written to different directory', default='-')
    parser_viz.set_defaults(func=spladder_viz)

    ### RUN MODE "PREP"
    parser_prep = subparsers.add_parser('prep', help='run mode to perform input preparations')

    align = parser_prep.add_argument_group('ALIGNMENTS')
    align.add_argument('-b', '--bams', dest='bams', metavar='FILE1,FILE2,...', help='alignment files in BAM format (comma separated list)', default='-', required=False)
    align.add_argument('--sparse-bam', dest='sparse_bam', action='store_true', help='store BAM content as sparse representation for later use [off]', default=False)
    align.add_argument('-n', '--readlen', dest='readlen', metavar='INT', type=int, help='read length (used for automatic confidence levele settings) [36]', default=50)
    align.add_argument('-c', '--confidence', dest='confidence', metavar='INT', type=int, help='confidence level (0 lowest to 3 highest) [3]', default=3)
    align.add_argument('--primary-only', dest='primary_only', action='store_true', help='only use primary alignments [on]', default=True)
    align.add_argument('--no-primary-only', dest='primary_only', action='store_false', default=True)
    align.add_argument('--var-aware', dest='var_aware', action='store_true', help='alignment files are variation aware (presence of XM and XG tags) [off]', default=False)
    align.add_argument('--no-var-aware', dest='var_aware', action='store_false', default=False)
    align.add_argument('--set-mm-tag', dest='mm_tag', help='sets the sequence of the mismatch tag used in alignments [NM]', default='NM')
    align.add_argument('--ignore-mismatches', dest='ignore_mismatches', action='store_true', help='does not filter by edit operations - does not require NM in BAM [off]', default=False)
    align.add_argument('--reference', dest='ref_genome', metavar='FILE', help='reference genome (only needed for CRAM file de-compression or consensus filtering)', default='')

    anno = parser_prep.add_argument_group('ANNOTATION')
    anno.add_argument('-a', '--annotation', dest='annotation', metavar='FILE', help='file name for annotation in GTF/GFF3 or format', default='-', required=False)
    anno.add_argument('--filter-overlap-genes', dest='filter_overlap_genes', action='store_true', help='remove genes from annotation that overlap each other [off]', default=False)
    anno.add_argument('--filter-overlap-exons', dest='filter_overlap_exons', action='store_true', help='remove exons from annotation that overlap each other [off]', default=False)
    anno.add_argument('--filter-overlap-transcripts', dest='filter_overlap_transcripts', action='store_true', help='remove transcripts from annotation that overlap each other [off]', default=False)
    parser_prep.set_defaults(func=spladder_prep)

    general_prep = parser_prep.add_argument_group('GENERAL')
    general_prep.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='user verbose output mode [off]', default=False)
    general_prep.add_argument('--parallel', dest='parallel', metavar='<INT>', type=int, help='use INT processors [1]', default=1)
    general_prep.add_argument('--tmp-dir', dest='tmpdir', metavar='DIR', help='directory to store temporary data [./tmp]', default='tmp')

    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)
    if len(argv) < 3:
        if argv[1] == 'build':
            parser_build.print_help()
        elif argv[1] == 'test':
            parser_test.print_help()
        elif argv[1] == 'viz':
            parser_viz.print_help()
        elif argv[1] == 'prep':
            parser_prep.print_help()
        else:
            parser.print_help()
        sys.exit(2)

    return (argv[1], parser.parse_args(argv[1:]))


def check_options(options, mode):
    
    if mode in 'build':
        if len(options.filter_consensus) > 0:
            if len(options.ref_genome) == 0:
                sys.stderr.write('\nERROR: using --filter-consensus requires setting a reference via --reference\n')
                sys.exit(1)
            if not options.filter_consensus in ['strict', 'lenient']:
                sys.stderr.write('\nERROR: --filter-consensus only allows the following choices: strict, lenient\n')
                sys. exit(1)

   
def main(argv=sys.argv):

    ### get command line options
    (mode, options)  = parse_options(argv)
    check_options(options, mode)
    options.func(options)

if __name__ == "__main__":
    main(sys.argv)
