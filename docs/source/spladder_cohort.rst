.. _spladder_cohorts:

Use on large cohorts
====================

While SplAdder is often run on a smaller set of samples, it can also be applied to larger cohorts
containing hundreds or thousands of samples. In this setting, it is often advisable to distribute
the computation over a high-performance compute cluster and use some workflow management framework
(such as `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_, `Nextflow
<https://www.nextflow.io/>`_, or even bash) to coordinate the individual steps.

In the following, we will provide a short overview on how the computation of splicing graphs and
their quantification can be split up into individual parts, so they can be run independently. Please
note, that the output directory for the whole computation is project specific and stays the same for
all (even parallel) runs. SplAdder will not re-compute existing files and when called on different
input files, name the outputs accordingly.

The process can be separated into four subsequent logical steps:

1. **Single graphs**: Creating an individual splicing graph per input sample
2. **Merged graph**: Merging all individual graphs into a joint graph representation
3. **Quantification**: Quantifying edges and nodes in the joint graph on each individual input sample
4. **Event Calling**: Calling events (optionally perform testing) on the joint, quantified graph.

For the description of the following steps, we will create a small hypothetical setup of samples,
that will be used throughout all commands. Assume that we have a cohort of 10 samples `S1` to `S10`.
The aligned (and indexed) bam files for the samples are available as ``S1.bam`` ... ``S10.bam``. We
are operating on a given annotation file ``annotation.gtf`` and our SplAdder project directory that
will hold all results will be ``spladder_out``. 

We will use simple ``bash`` commands to emulate the distribution of individual tasks. Please note
that the code as stated here, would just sequentially compute all single tasks and not generate any
parallelization benefit. For this to happen, you would need to submit the individual tasks to a
compute cluster or similar. (If you have a machine with many cores available, you could also
trivially parallelize by running several tasks in the background simultaneously.)

1. Single graphs
^^^^^^^^^^^^^^^^
In this first step, we will generate a splicing graph for each input sample. 

.. note:: In general a splicing graph is generated per gene. We will abstract from this here and
          only describe how the graphs will be treated across samples. This implicitly means that
          this is done per gene.

For each of the given input samples `Si`, we invoke the splicing graph generation separately::

    for i in $(seq 1 10)
    do
        spladder build -o spladder_out \
                       -a annotation.gtf \
                       -b S${i}.bam \
                       --merge-strat single \
                       --no-extract-ase \
                       --no-quantify-graph
    done

This will create an individual splicing graph in ``spladder_out/spladder`` for each sample. Please
note that we added the ``--no-extract-ase`` and ``--no-quantify-graph`` options here. This is to prevent 
SplAdder from quantifying the per-sample graphs and from automatically proceeding as is done by default 
in the smaller cohort analyses. With this option present, we gain a more fine-grained controlled over 
the graph building process. This option will also be present in the subsequent steps.

2. Merged graph
^^^^^^^^^^^^^^^

As a subsequent step, we now need to integrate the graphs across samples, to form one merged
splicing graph (per gene). We can just invoke SplAdder again on the same output directory now using
a different merging strategy::

    spladder build -o spladder_out \
                   -a annotation.gtf \
                   -b S1.bam,S2.bam,S3.bam,S4.bam,S5.bam,S6.bam,S7.bam,S8.bam,S9.bam,S10.bam \
                   --merge-strat merge_graphs \
                   --no-extract-ase \
                   --no-quantify-graph

Please note that now all alignment files of the cohort need to be provided. For larger cohorts it is
useful to collect all alignment files in a separate files, e.g. ``alignments.txt``. Then the merging
step could also be invoked as follows::

    spladder build -o spladder_out \
                   -a annotation.gtf \
                   -b alignments.txt \
                   --merge-strat merge_graphs \
                   --no-extract-ase \
                   --no-quantify-graph

Also not, that again we have provided the ``--no-extract-ase`` and ``--no-quantify-graph`` options,
as we are only interested in merging the sample-graphs, but not quantify per sample, yet, as this
can later be done in parallel.

As a last point, the merging of single sample-graphs is done progressively over samples and takes
some time. This is also hard to parallelize over genes, as there is the option to mere neighboring
genes, if they become too overlapping. If you really have many samples to merge (like hundreds or
thousands), SplAdder offers the possibility to parallelize the merge over chunks of samples, using
the option ``--chunked-merge LEVEL MAX_LEVEL START END``, where ``LEVEL`` is the current level of
the merging tree, ``MAX_LEVEL`` is the height of the merging tree, and ``START`` and ``END``
describe the node range at the current merge level. It is easiest to describe this using a small
example. Assume we have 100 single sample graphs to be merged. Instead of merging all 100 samples
progressively in one go, we could always merge 10 samples at a time in parallel. This will give us
10 ten-sample graphs. These 10 graphs can now again be merged, resulting in the 100-sample graph.
The total number of levels (``MAX_LEVEL``) can be computed as the ceiling of the chunksize-base
logarithm of the total number of samples. For our example this would be log_10 100 = 2.

Applying this to our 100 sample example::
    
    ### first level, merging samples in groups of 10
    for i in $(seq 0 10 90)
    do
        spladder build -o spladder_out \
                   -a annotation.gtf \
                   -b alignments.txt \
                   --merge-strat merge_graphs \
                   --no-extract-ase \
                   --no-quantify-graph \
                   --chunked-merge 0 1 ${i} $((${i} + 10))
    done

    ### second level, merging the 10 ten-sample graphs
    spladder build -o spladder_out \
               -a annotation.gtf \
               -b alignments.txt \
               --merge-strat merge_graphs \
               --no-extract-ase \
               --no-quantify-graph \
               --chunked-merge 1 1 0 10

Per default, the chunk size is 10, but can be adapted via the ``--chunksize`` option.

3. Quantification
^^^^^^^^^^^^^^^^^

Having the merged graph at hand, we can now proceed to quantifying nodes and edges of the graph
based on the alignment data. Each quantification will be done independently::

    for i in $(seq 1 10)
    do
        spladder build -o spladder_out -a annotation.gtf -b S${i}.bam \
                       --merge-strat merge_graphs \
                       --no-extract-ase \
                       --quantify-graph \
                       --qmode single
    done

Please note that now the merging strategy is still ``merge_graphs``, as we are quantifying the
merged graph and not the individual sample graphs. Also note that we have added the ``--qmode
single`` option.

As a second step to this phase, we need to collect the individual quantifications and aggregate them
in a joint database::

    spladder build -o spladder_out \
                   -a annotation.gtf \
                   -b S1.bam,S2.bam,S3.bam,S4.bam,S5.bam,S6.bam,S7.bam,S8.bam,S9.bam,S10.bam \
                   --merge-strat merge_graphs \
                   --no-extract-ase \
                   --quantify-graph \
                   --qmode collect

4. Event Calling
^^^^^^^^^^^^^^^^

Now one can proceed analog to the analysis of smaller cohorts. The joint graph is fully quantified
and we can move on to use it for downstream analyses, for instance to extract exon skipping events::

    spladder build -o spladder_out \
                   -a annotation.gtf \
                   -b S1.bam,S2.bam,S3.bam,S4.bam,S5.bam,S6.bam,S7.bam,S8.bam,S9.bam,S10.bam
                   --event-types exon_skip

In the above call we have omitted the ``--no-extract-ase`` option and SplAdder will automatically
proceed to this step. As all the intermediate quantification steps are already done, no step will be
carried out twice.

General Notes
^^^^^^^^^^^^^

When I/O is an issue, SplAdder has the option to generate a compressed summary for each input
alignment file. The information contained in that summary is comparable to a wiggle file but has
also information on the introns. Using this format will need some additional disk space, but allows
SplAdder to perform quantification and querying of intron coverage much more efficiently. You can
switch on the use of alignment summaries by::

    spladder build ... --sparse-bam ...

