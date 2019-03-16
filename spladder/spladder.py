import sys

from .spladder_build import spladder
from .spladder_viz import spladder_viz
from .spladder_test import spladder_test
from .rproc import spladder_pyproc 

def parse_options(argv):

    """Parses options from the command line """

    from argparse import ArgumentParser

    parser = ArgumentParser(prog='spladder')
    subparsers = parser.add_subparsers(help='Running modes', metavar='{build, test, viz}')

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
    inputs.add_argument('-n', '--readlen', dest='readlen', metavar='INT', type=int, help='read length (used for automatic confidence levele settings) [36]', default=36)
    inputs.add_argument('--primary-only', dest='primary_only', action='store_true', help='only use primary alignments [on]', default=True)
    inputs.add_argument('--no-primary-only', dest='primary_only', action='store_false', default=True)
    inputs.add_argument('--var-aware', dest='var_aware', action='store_true', help='alignment files are variation aware (presence of XM and XG tags) [off]', default=False)
    inputs.add_argument('--no-var-aware', dest='var_aware', action='store_false', default=False)
    inputs.add_argument('--labels', dest='labels', metavar='STRING', help='use labels instead of bam file names (comma separated list) [-]', default='-')
    #inputs.add_argument('-S', '--ref-strain', dest='refstrain', metavar='STRING', help='reference strain [-]', default='-')
    #inputs.add_argument('-x', '--same-genome', dest='same_genome', metavar='y|n', help='input alignments share the same genome [y]', default='y')
    inputs.add_argument('--input-graph', dest='spladderfile', metavar='FILE', help='use existing SplAdder graph as input (advanced) [-]', default='-')
    inputs.add_argument('--filter-overlap-genes', dest='filter_overlap_genes', action='store_true', help='remove genes from annotation that overlap each other [off]', default=False)
    inputs.add_argument('--filter-overlap-exons', dest='filter_overlap_exons', action='store_true', help='remove exons from annotation that overlap each other [off]', default=False)
    inputs.add_argument('--filter-overlap-transcripts', dest='filter_overlap_transcripts', action='store_true', help='remove transcripts from annotation that overlap each other [off]', default=False)
    inputs.add_argument('--ignore-mismatches', dest='ignore_mismatches', action='store_true', help='does not filter by edit operations - does not require NM in BAM [off]', default=False)
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
    graph = parser_build.add_argument_group('GRAPH OPTIONS')
    graph.add_argument('-c', '--confidence', dest='confidence', metavar='INT', type=int, help='confidence level (0 lowest to 3 highest) [3]', default=3)
    graph.add_argument('-I', '--iterations', dest='insert_intron_iterations', metavar='INT', type=int, help='number of iterations to insert new introns into the graph [5]', default=5)
    graph.add_argument('-M', '--merge-strat', dest='merge', metavar='<STRAT>', help='merge strategy, where <STRAT> is one of: single, merge_bams, merge_graphs, merge_all [merge_graphs]', default='merge_graphs')
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
    splice = parser_build.add_argument_group('SPLICE OPTIONS')
    splice.add_argument('--event-types', dest='event_types', metavar='STRING', help='list of alternative splicing events to extract [exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mtex_exons]', default='exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons')
    splice.add_argument('--extract-ase', dest='extract_as', action='store_true', help='extract alternative splicing events [on]', default=True)
    splice.add_argument('--no-extract-ase', dest='extract_as', action='store_false', default=True)
    splice.add_argument('--ase-edge-limit', dest='detect_edge_limit', metavar='INT', help='max number of edges in the graph to still extract events for a gene [500]', default=500)
    splice.add_argument('--curate-alt-prime', dest='curate_alt_prime', action='store_true', help='curate alt prime events [on]', default=True)
    splice.add_argument('--no-curate-alt-prime', dest='curate_alt_prime', action='store_false', default=True)
    splice.add_argument('--quantify-graph', dest='quantify_graph', action='store_true', help='quantify graph (implicitly on when --extract-ase is used) [off]', default=False)
    splice.add_argument('--no-quantify-graph', dest='quantify_graph', action='store_false', default=False)
    experimental = parser_build.add_argument_group('EXPERIMENTAL - BETA STATE')
    experimental.add_argument('--pyproc', dest='pyproc', action='store_true', help='use parallel implementation [off]', default=False)
    experimental.add_argument('--environment', dest='environment', metavar='STRING', help='conda environment to by used for pyproc', default=None)
    #experimental.add_argument('-R', '--replicates', dest='replicates', metavar='1,1,2,2,...', help='replicate structure of files (same number as alignment files) [all 1 - no replicated]', default='-')
    experimental.add_argument('--intron-cov', dest='intron_cov', action='store_true', help='count intron coverage [off]', default=False)
    experimental.add_argument('--qmode', dest='qmode', metavar='STRING', help='quantification mode: single, collect, all [all]', default='all')
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
    inputs_test.add_argument('--subset-samples', dest='subset_samples', action='store_true', help='gene expression counting will be only done on the tested subset of samples [off]', default=False)
    inputs_test.add_argument('--no-subset-samples', dest='subset_samples', action='store_false', default=False)
    testing = parser_test.add_argument_group('TESTING OPTIONS')
    testing.add_argument('-C', '--correction', dest='correction', metavar='STR', help='method for multiple testing correction (BH, Bonferroni, Holm, Hochberg, BY, TSBH) [BH]', default='BH')
    testing.add_argument('-0', '--max-zero-frac', dest='max_0_frac', metavar='FLOAT', type=float, help='max fraction of 0 values per event isoform quantification over all tested samples [0.5]', default=0.5)
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
    experimental_test.add_argument('--low-memory', dest='low_memory', action='store_true', help='use less memory at the cost of longer running time [off]', default=False)
    parser_test.set_defaults(func=spladder_test)

    ### RUN MODE "VIZ"
    parser_viz = subparsers.add_parser('viz', help='run mode to visualize graphs and events')
    required_viz = parser_viz.add_argument_group('MANDATORY')
    required_viz.add_argument('-o', '--outdir', dest='outdir', metavar='DIR', help='spladder directory containing the spladder results', default='-', required=True)
    optional_viz = parser_viz.add_argument_group('OPTIONAL')
    optional_viz.add_argument('-b', '--bams', dest='bams', metavar='FILE1A,FILE2A:FILE1B,FILE2B,,...', help='alignment files in BAM format (comma separated list,colon separated groups)', default='-')
    optional_viz.add_argument('-L', '--labels', dest='labels', metavar='LABEL_A,LABEL_B,...', help='group labels for alignment files groups (comma separated list)', default='')
    optional_viz.add_argument('-g', '--gene-name', dest='gene_name', metavar='STR', help='gene_name to be plotted', default=None)
    optional_viz.add_argument('-e', '--event-id', dest='event_id', metavar='STR', help='event to be plotted', default=None)
    optional_viz.add_argument('--test-result', dest='test_result', metavar='INT', type=int, help='plot top k significant events from test', default=0)
    optional_viz.add_argument('--test-labels', dest='test_labels', metavar='STR', type=str, help='labels used for the groups in the test (order matters) [condA:condB]', default='condA:condB')
    optional_viz.add_argument('-t', '--event-types', dest='event_types', metavar='EVENT1,EVENT2,...', help='list of alternative splicing events to extract [exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons]', default='exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons')
    optional_viz.add_argument('--testdir', dest='testdir', metavar='DIR', help='directory to testing output, if different from spladder outdir', default='-')

    output_viz = parser_viz.add_argument_group('OUTPUT')
    output_viz.add_argument('-m', '--mincount', dest='mincount', metavar='INT', type=int, help='minimum count of introns to be displayed in coverage plot [0]', default=0)
    output_viz.add_argument('-f', '--format', dest='format', metavar='STR', help='plot file format [pdf, png, d3]', default='pdf')
    output_viz.add_argument('--zoom-x', dest='zoom_x', metavar='percent_left,percent_right', help='zoom x axis from percent_left to percent_right [0.0,1.0]', default='0.0,1.0')
    output_viz.add_argument('-l', '--log', dest='log', action='store_true', help='plot coverage information in log scale [off]', default=False)

    user_viz = parser_viz.add_argument_group('USER')
    user_viz.add_argument('-u', '--user', dest='user', action='store_true', help='apply user mode (experimental) [off]', default=False)
    user_viz.add_argument('-T', '--transcripts', dest='transcripts', action='store_true', help='plot annotated transcripts', default=False)
    user_viz.add_argument('-s', '--splicegraph', dest='splicegraph', action='store_true', help='plot splicegraph structure', default=False)

    general_viz = parser_viz.add_argument_group('GENERAL')
    general_viz.add_argument('-c', '--confidence', dest='confidence', metavar='INT', type=int, help='confidence level (0 lowest to 3 highest) [3]', default=3)
    general_viz.add_argument('-V', '--validate-sg', dest='validate_sg', action='store_true', help='validate splice graph [off]', default=False)
    general_viz.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='user verbose output mode [off]', default=False)
    general_viz.add_argument('-d', '--debug', dest='debug', action='store_true', help='use debug mode [off]', default=False)
    parser_viz.set_defaults(func=spladder_viz)

    ### RUN MODE "VIZ"
    parser_pyproc = subparsers.add_parser('pyproc')
    required_pyproc = parser_pyproc.add_argument_group('MANDATORY')
    required_pyproc.add_argument('proc', metavar='PROC_FILE', help='pyproc pickle file', default='-')
    required_pyproc.add_argument('data', metavar='DATA_FILE', help='pyproc data file', default='-')
    parser_pyproc.set_defaults(func=spladder_pyproc)

    options = parser.parse_args(argv[1:])

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
        else:
            parser.print_help()
        sys.exit(2)

    return options

def main(argv=sys.argv):

    ### get command line options
    options = parse_options(argv)
    options.func(options)

if __name__ == "__main__":
    main(sys.argv)
