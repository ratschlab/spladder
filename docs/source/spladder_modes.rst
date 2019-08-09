SplAdder - the run modes
========================

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

:ref:`Graph construction <graph_construction>`
    This is the very initial phase. It parses the given annotation file and summarizes all
    transcripts of a gene into a splicing graph. This graph will be the basis for all further steps
    in the workflow.
:ref:`Graph augmentation <graph_augmentation>`
    Given at least one alignment file, the splicing graph of each gene is augmented with new introns
    and exon segments that were detected in the alignment file. There are different ways how a more
    than one input alignment files can be combined into final splicing graphs. At the end of this
    phase, each gene contains an augmented graph that carries not only annotated splice connections
    but also any novel connections found in the data. Depending on the chosen confidence level, this
    graph will have a higher or lower density.
:ref:`Graph quantification <graph_quantification>`
    Once a graph is constructed, all nodes and edges (exons and introns, respectively) in the graph
    can be quantified using at least one input alignment file. The quantification values can then be
    used subsequently to quantify splicing events and to compute percent spliced in (PSI) values.
:ref:`Event detection <event_detection>`
    Based on the splicing graph of each gene, SplAdder can detect different types of alternative
    splicing events: exon skipping, intron retention, alternative 3' splice sites, alternative 5'
    splice sites, mutual exclusive exons and multiple (coordinated) exon skips. Each event can be
    quantified using the graph quantifications from the previous step.

In the following, we will provide more in-depth information for each of the phases and describe how
the result can be influenced through the choice of command line parameters.

.. _graph_construction:

Graph construction
^^^^^^^^^^^^^^^^^^

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

    spladder build .. --filter-overlap-transcripts ...

If two genes overlap on the same strand::

    spladder build .. --filter-overlap-genes

.. _graph_augmentation:

Graph augmentation
^^^^^^^^^^^^^^^^^^

The augmentation phase brings together alignment file and splicing graphs. Let's assume that you are
given an alignment file ``alignment.bam`` (which should also have an index ``alignment.bam.bai``)
and an annotation file in GTF format ``annotation.gtf``. You can the simply invoke::

    spladder build --bams alignment.bam --annotation annotation.gtf --outdir spladder_out 

All three parameters are mandatory for a SplAdder run in ``build`` mode. Due to the default values
of other parameters, this will carry out a full run of all phases. We will describe in the
following, which parameters you can change to either only run this phase or to adapt how the
splicing graph will be augmented. 

**Alignment**
    By default, SplAdder only uses primary alignments (in SAM/BAM the ones not carrying the 256
    bit-flag). This can be changed by also allowing for secondary alignments to be used::

        spladder build ... --no-primary-only ...

    The quality of an alignment is partially determined by the number of mismatches it carries. The
    default tag in SAM/BAM for this is the ``NM:i:`` tag. Te let SplAdder use a different tag, such
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


.. _graph_quantification:

Graph quantification
^^^^^^^^^^^^^^^^^^^^

Text

.. _event_detection:

Event detection
^^^^^^^^^^^^^^^

Text
