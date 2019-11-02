.. _spladder_run_modes:

Run modes
=========

SplAdder has different run modes that reflect the different steps of a typical analysis pipeline:

``build`` mode
    for constructing splicing graphs from RNA-Seq data and extracting alternative events
``test`` mode
    for the differential analysis between samples
``viz`` mode
    for the visualization of splicing graphs and alternative events

In the following, we will give a short overview of the different modes and how to use them. Special
use cases, for instance the handling of large sample cohorts, will be discussed as a separate topic.

The ``build`` mode
------------------

The ``build`` mode is the basic run mode in SplAdder. It is used to construct splicing graphs and
to extract alternative splicing events.

To display all available options for ``build``, one can simply type::

    spladder build --help

This first step of any SplAdder pipeline consists of several main phases (some of which can be
omitted) :

:ref:`1 Graph construction <graph_construction>`
    This is the very initial phase. It parses the given annotation file and summarizes all
    transcripts of a gene into a splicing graph. This graph will be the basis for all further steps
    in the workflow.
:ref:`2 Graph augmentation <graph_augmentation>`
    Given at least one alignment file, the splicing graph of each gene is augmented with new introns
    and exon segments that were detected in the alignment file. There are different ways how a more
    than one input alignment files can be combined into final splicing graphs. At the end of this
    phase, each gene contains an augmented graph that carries not only annotated splice connections
    but also any novel connections found in the data. Depending on the chosen confidence level, this
    graph will have a higher or lower density.
:ref:`3 Graph quantification <graph_quantification>`
    Once a graph is constructed, all nodes and edges (exons and introns, respectively) in the graph
    can be quantified using at least one input alignment file. The quantification values can then be
    used subsequently to quantify splicing events and to compute percent spliced in (PSI) values.
:ref:`4 Event detection <event_detection>`
    Based on the splicing graph of each gene, SplAdder can detect different types of alternative
    splicing events: exon skipping, intron retention, alternative 3' splice sites, alternative 5'
    splice sites, mutual exclusive exons and multiple (coordinated) exon skips. Each event can be
    quantified using the graph quantifications from the previous step.

In the following, we will provide more in-depth information for each of the phases and describe how
the result can be influenced through the choice of command line parameters.

.. _graph_construction:

1 Graph construction
^^^^^^^^^^^^^^^^^^^^

This phase runs implicitly before any other phase. We just describe it here for completeness, but
in general there is no reason to run this phase only by itself. What it does in the background,
though, is to transform the given annotation file::

    spladder build .. --annotation annotation.gtf ...

into a SplAdder specific format, containing all transcript information and the initial splicing
graphs per gene. These information will be stored at the same location as ``annotation.gtf`` and is
identified by the suffix ``.pickle``. The resulting file in our example would be named
``annotation.gtf.pickle``. Depending on the settings, additional files might be created, for
instance to mask out certain regions from the annotation.
This step is only performed once per annotation file. The summary files will then be re-used by any
subsequent SplAdder run using the same annotation file.

The user can influence how SplAdder uses the annotation information in certain situations of
ambiguity. However, none of these options is set by default.

In cases of annotations overlapping on the same strand, one can remove the annotation on three
different levels.

If two exons of different genes overlap on the same strand, one can remove them with::

    spladder build ... --filter-overlap-exons ...

If two transcripts of different genes overlap on the same strand, one can remove them with::

    spladder build ... --filter-overlap-transcripts ...

If two genes overlap on the same strand::

    spladder build ... --filter-overlap-genes ...

.. _graph_augmentation:

2 Graph augmentation
^^^^^^^^^^^^^^^^^^^^

The augmentation phase brings together alignment file and splicing graphs. Let's assume that you are
given an alignment file ``alignment.bam`` (which should also have an index ``alignment.bam.bai``)
and an annotation file in GTF format ``annotation.gtf``. You can the simply invoke::

    spladder build --bams alignment.bam \
                   --annotation annotation.gtf \
                   --outdir spladder_out 

All three parameters are mandatory for a SplAdder run in ``build`` mode. Due to the default values
of other parameters, this will carry out a full run of all phases. We will describe in the
following, which parameters you can change to either only run this phase or to adapt how the
splicing graph will be augmented. 

Multiple alignment files can be provided using comma-separated notation::

    spladder build --bams alignment1.bam,alignment2.bam,...

Alternatively, a text file, e.g., ``alignment_list.txt``, can be provided. This should contain the
absolute path to one alignment file per line. The filename has to end in ``.txt``. SplAdder can then
be invoked with::
    
    spladder build --bams alignment_list.txt

In its latest version, SplAdder also supports (on an experimental level) CRAM compressed alignment
files as input. If you are using such files, in addition to the input filenames of the
alignment files, also the path to the indexed reference sequence used for compression is required::

    spladder build --bams alignment1.cram,alignment2.cram,... --cram-reference path/to/cram_ref.fa 

**Alignment**
    By default, SplAdder only uses primary alignments (in SAM/BAM the ones not carrying the 256
    bit-flag). This can be changed by also allowing for secondary alignments to be used::

        spladder build ... --no-primary-only ...

    The quality of an alignment is partially determined by the number of mismatches it carries. The
    default tag in SAM/BAM for this is the ``NM:i:`` tag. To let SplAdder use a different tag, such
    as ``Nm:i:``, one can use::
        
        spladder build ... --set-mm-tag Nm ...

    Alternatively, one can also force SplAdder not to use any mismatch information (this is not
    recommended)::

        spladder build ... --ignore-mismatches ...
    
**Augmentation**
    Different types of augmentations are possible. The majority of them is switched on by default.
    For instance the insertion of new intron retentions is always carried out. To switch this step
    off, one would add::

        spladder build ... --no-insert-ir ...
    
    Similarly, the addition of novel cassette exons is also on by default. To switch this step off,
    one would add::

        spladder build ... --no-insert-es ...

    Also the addition of novel intron edges is switched on by default. To switch it off, one would
    add::

        spladder build ... --no-insert-ni ...

    On the other hand, additional steps for graph cleaning are not switched on by default. For
    instance the removal of exons shorter than 9nt from the graph can be add with::

        spladder build ... --remove-se ...

    Lastly, as SplAdder is a heuristic framework, the addition of novel nodes and edges to the graph
    depends on the input order of new introns and on the current state of the graph (that is the
    nodes and edges already present). To increase sensitivity, the addition of new intron edges is
    iterated a certain number of times (per default 5 times). One can increase the number if
    iterations, for instance to 10, by::

        spladder build ... --iterations 10 ...

**Confidence**
    The confidence level of a SplAdder run determines how strongly input alignments are filtered
    before new nodes and edges are added to the splicing graphs. In general, there are four
    confidence levels, with confidence increasing from 0 to 3. The default level is 3 and applies
    the highest level of filtering. To adapt this choice, e.g., to confidence level 2, one can use::

        spladder build ... --confidence 2 ...

    The read filter criteria are dependent on the read length. Here a short overview of the criteria
    for each of the levels:

    +----------+------------------------------+---------------------------------+
    | Level    | Criteria                     | Value                           |
    +==========+==============================+=================================+
    |        3 | Maximum number of mismatches | 0                               |
    +----------+------------------------------+---------------------------------+
    |        3 | Minimum number of alignments | 2                               |
    +----------+------------------------------+---------------------------------+
    |        3 | Minimum anchor length        | ceil(0.25 * readlength)         |
    +----------+------------------------------+---------------------------------+
    |        3 | Maximum intron length        | 350000                          |
    +----------+------------------------------+---------------------------------+
    +----------+------------------------------+---------------------------------+
    |        2 | Maximum number of mismatches | max(1, floor(0.01 * readlength) |
    +----------+------------------------------+---------------------------------+
    |        2 | Minimum number of alignments | 2                               |
    +----------+------------------------------+---------------------------------+
    |        2 | Minimum anchor length        | ceil(0.20 * readlength)         |
    +----------+------------------------------+---------------------------------+
    |        2 | Maximum intron length        | 350000                          |
    +----------+------------------------------+---------------------------------+
    +----------+------------------------------+---------------------------------+
    |        1 | Maximum number of mismatches | max(1, floor(0.02 * readlength) |
    +----------+------------------------------+---------------------------------+
    |        1 | Minimum number of alignments | 2                               |
    +----------+------------------------------+---------------------------------+
    |        1 | Minimum anchor length        | ceil(0.15 * readlength)         |
    +----------+------------------------------+---------------------------------+
    |        1 | Maximum intron length        | 350000                          |
    +----------+------------------------------+---------------------------------+
    +----------+------------------------------+---------------------------------+
    |        0 | Maximum number of mismatches | max(2, floor(0.03 * readlength) |
    +----------+------------------------------+---------------------------------+
    |        0 | Minimum number of alignments | 1                               |
    +----------+------------------------------+---------------------------------+
    |        0 | Minimum anchor length        | ceil(0.10 * readlength)         |
    +----------+------------------------------+---------------------------------+
    |        0 | Maximum intron length        | 350000                          |
    +----------+------------------------------+---------------------------------+

    In the above table, the `maximum number of mismatches` is used to remove reads that have low
    quality alignments, the `minimum number of alignments` is the number of split/spliced alignments
    necessary to confirm a new intron edge for being taken into the graph, the `minimum achor
    length` is the shortest overlap to an exon segment that a split/spliced alignment needs to have
    to be counted towards confirming an intron, and the `maximum intron length` is the upper
    threshold for new introns to be counted.

**Merging**
    As SplAdder can be run with multiple alignment files as input, there are several ways on how
    these files can be combined into forming augmented splicing graphs. This behavior is controlled
    with the setting of the `merging strategy` using ``--merge-strat``.

    The first way of merging is to generate a separate augmented splicing graph per given input
    alignment file. This strategy is called `single` and can be invoked as follows::

        spladder build ... --merge-strat single ...

    The second (and default) way of merging is to create a single splicing graph per input file and
    then merge all graphs into a joint single graph. (This happens for every gene independently.)
    This strategy is called `merge graphs` and can be invoked as follows::

        spladder build ... --merge-strat merge_graphs ...

    A third way of merging is to treat all input alignment files as technical replicates and
    directly form a splicing graph using all reads. (This makes a difference especially for the
    count thresholds.) This strategy is called `merge bams` and can be invoked as follows::

        spladder build ... --merge-strat merge_bams ...

    The fourth way of merging is a combination of ``merge_bams`` and ``merge_graphs``. In this
    setting, both steps are performed and both resulting graphs are integrated into a joint graph.
    The idea behind this setting is to generate maximum sensitivity. However, the improvement is in
    general marginal and we would not advise to use this setting in general. If you would like to
    try it nevertheless, you can do so with::

        spladder build ... --merge_strat merge_all ...

**Validation**
    SplAdder has the option to validate edges in the graph. This is relevant when working on larger
    cohorts of samples. In this filtering step an edge is removed if it is not present in the
    initial annotation and is supported in less than a certain number of input samples. By default
    this threshold is 10 or the number of input samples in cases where less than 10 samples are
    given. The threshold can be adapted using ``--validate-sg-count``. If nodes get orphaned
    through the pruning process, they will be also removed from the graph. Following an example that
    removes all edges from the graph that are present in less than 5 input samples::

        spladder build ... --validate-sg --validate-sg-count 5 ...


.. _graph_quantification:

3 Graph quantification
^^^^^^^^^^^^^^^^^^^^^^

In the step of graph quantification, the augmented graph is evaluated again against all given input
alignment files, to determine edge and node weights based on the respective expression. If
alternative splicing events are to be extracted (next step), this step is carried out automatically.
If the user decided not to extract alternative splicing events (explained in the next section), but
the graph should be quantified anyways, this can be achieved with::

    spladder build ... --quantify-graph ...

Especially for larger cohorts, it can be challenging to process through all the alignment files for
quantification. (We will provide more detailed explanations for this scenario in `Working with large
cohorts`.) Here, we will just mention, that the quantification step can be invoked in different
modes, called `qmodes`. Let us assume, that two alignment files were provided to SplAdder,
``aligment1.bam`` and ``alignment2.bam``. Then the default is that all files processed sequentially.
This quantification mode is called ``all`` and (despite being used implicitly per default), can also
be explicitly set with::

    spladder build ... --bams alignment1.bam,alignment2.bam \
                       --qmode all ...

As an alternative, one can also provide a single alignment file at a time to SplAdder. This strategy
is called ``single`` and can be used to parallelize SplAdder processes across alignment files. It
can be invoked via::

    spladder build .. --bams alignment1.bam --qmode single ...
    spladder build .. --bams alignment2.bam --qmode single ...

The ``single`` command always needs to be accompanied by an additional run of SplAdder, that
integrates the quantification files for the single alignment files into a joint data structure. 
For this, all alignment files are provided as input and the quantification mode ``collect`` is
chosen::

    spladder build .. --bams alignment1.bam,alignment2.bam \
                      --qmode collect ...

.. _event_detection:

4 Event detection
^^^^^^^^^^^^^^^^^

In this last phase of the ``build`` mode, the graphs are used for the extraction of alternative
splicing events. Event extraction is performed per default. The user can choose to omit this step
entirely (for instance to carry it out at a later point in time). This is done via::

    spladder build ... --no-extract-ase ...

SplAdder can currently extract 6 different types of alternative splicing events:

- exon skips (`exon_skip`)
- intron retentions (`intron_retention`)
- alternative 3' splice sites (`alt_3prime`)
- alternative 5' splice sites (`alt_5prime`)
- mutually exclusive exons (`mutex_exons`)
- multiple (coordinated) exons skips (`mult_exon_skip`)

Per default all events of all types are extracted from the graph. To specify a single type or a
subset of types (e.g., exon skips and mutually exclusive exons only), the user can specify the short
names of the event types (as shown in parentheses above) as follows::

    spladder build ... --event-types exon_skip,mutex_exons ...

In some cases (for instance when integrating hundreds of alignment samples), the splicing graphs can
grow very complex. To limit the running time, an upper bound for the maximum number of edges in the
splicing graph of a gene to be used for event extraction is set. This threshold is 500 per default.
To adapt this threshold, e.g., to 250, the user can specify::
    
    spladder build ... --ase-edge-limit 250 ...

The ``test`` mode
-----------------

This SplAdder mode is for differentially testing the usage of alternative event between two groups
of samples. A prerequisite for this is that all samples that are involved in testing have been
subjected to a joint analysis in the ``build`` mode. However, not the full set of samples collected
in the ``build`` mode has to be subjected to testing, but subsets of samples can be used instead. 

It is recommended that for each sample condition to be tested (e.g., wild type and some mutant), the
number of available replicates is at least three. Further, the mean-variance relationship for intron
counts are estimated on the set of tested events. It the number of events to be tested becomes too
small, then this estimate becomes unstable and might result in an error.

For the invocation of the testing mode, three different input parameters are mandatory::
    
    spladder test --conditionA aligmmentA1.bam,alignmentA2.bam \
                  --conditionB alignmentB1.bam,alignmentB2.bam \
                  --outdir spladder_out

In detail, these are the two lists of alignment files representing the samples for conditions A and
B, respectively, as well as the SplAdder output directory. This is the same output directory, as
has been used for the ``build`` mode.
Analog to the way a list of alignments can be provided in ``build`` mode, also in ``test`` mode the
comma-separated file list can be substituted with a file containing the paths to the respective
files::

    spladder test --conditionA alignmentsA_list.txt \
                  --conditionB alignmentsB_list.txt \
                  --outdir spladder_out

By default all event types will be subjected to testing (if they were extracted from the graph prior
to testing). If only a specific event type or subset of types should be tested, e.g., exon skips and
mutual exclusive exons, the same syntax as in build mode can be applied::

    spladder test ... --event-types exon_skip,mutex_exons ...

If you have built the SplAdder graphs using non-default setting, for instance an adapted confidence
level of 2, these parameters also need to be passed in ``test`` mode, so the correct input files are
chosen from the project directory::

    spladder test ... --confidence 2 ...

By default expression outliers are removed in a preprocessing step. If you would like to keep genes
that show outlier expression, this behavior can be disabled with::

    spladder test ... --no-cap-exp-outliers

Similarly, you can also switch on the capping of splice outliers, which is not done by default::

    spladder test ... --cap-outliers ...

Sometimes it is useful to assign labels to the two groups being tested, especially is multiple
different groupings are analyzed. Groups A and B can be assigned arbitrary labels, such as `Mutant`
and `Wildtype`,  using::

    spladder test ... --labelA Mutant --labelB Wildtype

In addition, you can also provide a separate tag that will be appended to the output directory name.
This is useful, if several rounds of testing or different parameter choices are explored. To tag the
output directory with `Round1` you would use::

    spladder test ... --out-tag Round1 ...

The ``test`` mode is capable of generating several summary plots for diagnosing issues and getting a
better understanding of the data being tested. Per default, the plots are generated in `png` format,
but other formats such as `pdf` or `eps` can be chosen as well. Per default, the diagnose plots are
switched off. To generate them, for instance in `pdf` format, you would use::

    spladder test ... --diagnose-plots --plot-format pdf ...

If several compute cores are available, the computation of the testing can be accelerated by
allowing parallel access. If 4 cores should be used::

    spladder test ... --parallel 4 ...

The ``viz`` mode
----------------

The purpose of this mode is to generate visual overviews of splicing graphs and events and the
associated coverage available in the underlying RNA-Seq samples.

General organisation
^^^^^^^^^^^^^^^^^^^^

In general, the plots are organized as individual tracks, which can be stacked to visualize several
sources of information jointly. Thereby, the order, number and repetition of tracks can be defined
by the user. This allows for the generation of simple overview plots as well as for more complex
multi-track visualizations. If more than one track is present, all tracks share the same joint
coordinate system on the x axis.

To determine which genomic range is plotted, all elements provided in ``--tracks`` are considered
and a region including all of them is generated. This logic can be overruled using the ``--range``
parameter to specify a specific range. However, there are also data track elements that do not
necessarily carry any range information (such as a coverage track). In this case the ``--range``
argument would be required. In the following, we will first explain the definition of tracks in more
detail and will then provide some information on how to define a specific range.


Data tracks
^^^^^^^^^^^

This parameter is concerned with defining which data tracks should be visualized in the plot and in
which order. The general syntax for specific a data track is as follows::

    spladder viz --track TYPE [TYPE_INFO [TYPE_INFO ...] ]

Here, ``TYPE`` describes one of the following possibilities (where ``TYPE_INFO`` is specifically
defined for each type):

    - **splicegraph** shows the structure of the splicing graph for each of the given genes. If no
      ``TYPE_INFO`` is provided, the gene(s) from the ``--range`` argument are used. To plot the
      splicing graph for gene with ID `gene1`, one would use::
        
        spladder viz --track splicegraph gene1

    - **transcript** shows the structure of all annotated transcripts for each of the given genes.
      If not ``TYPE_INFO`` is provided, the gene(s) from the ``--range`` argument are used. To plot
      the splicing graph for gene with ID `gene1`, one would use::

        spladder viz --track transcript gene1

    - **event** shows the structure of the given events, where each event can be specified by its
      ID. For instance to show the structure of events `exon_skip_2` and `alt_3prime_5`, one can
      use::

        spladder viz --track event exon_skip_2 alt_3prime_5

      If not a specific event ID is given but only the event type, all events of that type
      for the genes given in ``--range`` are shown. So to show all `exon_skip` events of gene
      `gene1`, the correct call would be::

        spladder viz --range gene gene1 --track event exon_skip

      If all events of a given gene should be shown, then one can use the special keyword ``any`` to
      achieve this::
        
        spladder viz --range gene gene1 --track event any

    - **coverage** shows the coverage information in the given range for all samples provided in
      ``TYPE_INFO``. To show coverage for samples `alignment1.bam` and `alignment2.bam`, one would
      use::

        spladder viz --track coverage alignment1.bam alignment2.bam

      If the coverages of both files should be added up, one can also define them as a group::

        spladder viz --track coverage alignment1.bam,alignment2.bam
      
      Sometimes it is useful to assign descriptive labels to single or multiple samples. Given the
      samples `alignment1.bam` - `alignment4.bam`, which can be separated into groups `wildtype` and
      `mutant`, respectively, one can use these labels in the plot as follows::

        spladder viz --track coverage \
                             wildtype:alignment1.bam,alignment2.bam \
                             mutant:alignment3.bam,alignment4.bam

    - **segments** shows the coverage information in the given range as internally used by SplAdder
      in the splicing graph, quantifying each exonic segment. The usage is analog to ``--coverage``. 

Order of multiple tracks
^^^^^^^^^^^^^^^^^^^^^^^^

The order of the tracks is determined by the order they are provided in at the command line. This is
true for both the order of keywords within a single ``--track`` parameter, as well as for the order
of multiple ``--track`` parameters. 

Let us consider the following example::
    
    spladder viz --range gene gene1 \
                 --track coverage,segments alignment1.bam,alignment2.bam \
                 --track event any \
                 --track splicegraph \
                 --track event exon_skip

This plot will have **five** tracks: coverage, segments, events (any), splicing graph, events (only
exon skips). This means, even the same track can be plotted multiple times, if requested.

Plotting range
^^^^^^^^^^^^^^

Using the ``--range`` parameter, the user determines exactly which genomic range is to be considered
for plotting the data tracks. This information can be provided as coordinates or the ID information
of one or many elements. The usage of ``--range`` overrules any range determined based on the
elements given via ``--tracks``. The syntax thereby is as follows::

    spladder viz --range TYPE TYPE_INFO [TYPE_INFO ...]

Here, ``TYPE`` describes one of the following possibilities (where ``TYPE_INFO`` is specifically
defined for each type):

    - **gene** allows for providing at least one gene ID to be considered. If multiple genes should
      be used, just list them after the ``gene`` keyword::

        spladder viz --range gene geneID1 geneID2

    - **event** allows for providing at least one event ID to be considered. If multiple events
      should be used, just list them after the ``event`` keyword::

        spladder viz --range event eventID1 eventID2

    - **coordinate** allows for specifying a coordinate range to be used. Here, the type info
      contains the list of coordinates to be used. As all ranges will be combined into a joint range
      eventually, there is little use in providing several coordinate ranges, as the union would be
      taken. For specifying the genome range of positions 100000 to 101000 on chr1, one would
      specify::

        spladder viz --range coordinate chr1 100000 101000

.. note:: The ``--range`` parameter can be used multiple times to combine several ranges. Please note that all provided ranges will be combined into a joint range including all other ranges before plotting. Also note that plotting ranges on different chromosomes is currently not supported as well as plotting ranges exceeding a total length of 1 000 000 bases.

Output names and formats
^^^^^^^^^^^^^^^^^^^^^^^^

The user has to choose an output file name for each plot generated. This is specified as the
basename of the output file, not containing the output directory or the file ending (which is chosen
based on the format). The relevant parameter for this is ``--outbase`` (or short ``-O``). The
default format of the plots is `pdf`, but any format supported by Matplotlib can be used. The
following two calls for using the output basename `mytest` and the format `png` are equivalent::

    spladder viz ... --outbase mytest --format png ...
    spladder viz ... -O mytest --format png ...

Please note that when using the special plotting mode ``--test`` and providing a test directory with
``--testdir`` (see below), the plots are not placed in the SplAdder output directory but in the
given test directory.

Plotting test results
^^^^^^^^^^^^^^^^^^^^^

For visualizing events based on the outcome of the testing mode, there is a special track mode
available, which is called ``--test``. In principle it works as the other tracks but follows a
specific syntax of its elements. There is a default set of 2 tracks that is generated using this
option: an event track showing the event of interest and a segements track showing the
quantification of segments used for the test for each of the two groups. The general structure is::

    --test TESTCASE EVENT_TYPE TOP_K

While any of the elements is optional, the order is important and elements can only be omitted at
the end but not in the middle.

Depending on how the groups are named in testing mode, the output can be found in different
subdirectories. So if you groups were `WT` and `MUT`, your output name used by SplAdder would be
`testing_WT_vs_MUT`. However, if you did not use any group names, you can just use `default`.
Following are two examples for using the default and the specific group mode, respectively::

    spladder viz ... --test default ...
    spladder viz ... --test testing_WT_vs_MUT ...

The ``EVENT_TYPE`` specifies the test result of which event type should be considered. For each of the
`k` top events, a separate plot will be generated. You can comma-separate multiple event types or
write `any`, for all event types. Here two examples for plotting exon skips and intron retentions or
any event, respectively::

    spladder viz ... --test default exon_skip,intron_retention ...
    spladder viz ... --test default any

Lastly, the user can specify the number of top events (following the ranking in the testing result
file) that should be plotted. If the value is omitted, the default of 1 is used. To plot for
instance the top 5 exon skip events, one would use::

    spladder viz ... --test default exon_skip 5

This will create 5 separate plots. The output name will have descriptive suffixes, to tell them
apart.

It can happen that the output for testing with SplAdder was written into a user-defined directory
and not into the default SplAdder output directory. In this case, the directory can be specified
using ``--testdir``. For instance, if the test results can be found in `mytestingdir`, the SplAdder
call would need to be adapted as follows::

    spladder viz ... --test default exon_skip 5 --testdir mytestingdir

As already noted earlier, this will also influence where the plots for the test are placed. For the
above example, all plots will be written to ``mytestingdir/plots/``.

.. note:: If in addition to ``--test`` further tracks are also defined with ``--track``, then each of the tracks is added to **each** of the plots generated for the test results.

