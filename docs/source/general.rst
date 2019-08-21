.. _general_info:

General information
===================

In the following, we provide some general information on how the setup of SplAdder works in
principle and mention several useful things to keep in mind when using the software. 

The SplAdder directory setup
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The general purpose of SplAdder is to build and quantify augmented splicing graphs from RNA-Sequencing
data and to utilize them for the detection for alternative splicing events. To achieve this,
SplAdder operates on a single `project` at a time, where a project is characterized by a single
shared output (or project) directory. All data of subsequent steps will be written to that
directory and a certain substructure will directly be created by SplAdder.

SplAdder has different run modes that reflect the different steps of a typical analysis pipeline:

``build`` mode
    for constructing splicing graphs from RNA-Seq data and extracting alternative events
``test`` mode
    for the differential analysis between samples
``viz`` mode
    for the visualization of splicing graphs and alternative events

All of these modes will operate on the same output directory. Please note, that the ``build`` mode
always has to precede the ``testing`` and ``viz`` modes, as this creates the splicing graph
structures the latter modes operate on.

Please have a look at the :ref:`SplAdder run modes <spladder_run_modes>` page for further information.

SplAdder is a heuristic
^^^^^^^^^^^^^^^^^^^^^^^

We would like to reiterate here that SplAdder is a heuristic approach, that employs a system of
empirical filter rules on RNA-Seq data to extend a splicing graph pre-defined by a given annotation.
Due to this heuristic nature, there is a list of things one should keep in mind when working with
the software:

- Although the algorithm is deterministic, it is sensitive to the order of input data. Especially
  when you are working with many input samples and integrate their information to generate a single
  splicing graph, the order of input files might influence the outcome. However, in most cases the
  generated splicing graphs are robust against changes of the order of the input.
- Depending on the annotation file that you are using, some annotated gene regions might be subject to
  filtering. SplAdder will automatically ignore all introns in the input data that could be assigned
  to more than one gene in the annotation. (This decision is strand-specific, that is both genes and
  the intron need to be on the same strand.) 
- As SplAdder will sample all possible events from the graph, there can be a certain redundancy in
  the events (although all events are unique as a whole). For instance, if there are 3 possible
  donors paired with an acceptor in a gene on the positive strand, SplAdder will output all three
  possible pairs of alternative 3 prime splice site events. 

Default parameters
^^^^^^^^^^^^^^^^^^

For most of the settings available in SplAdder default parameters are assumed. In a basic call in
build-mode (``spladder build``), SplAdder requires at least three parameters: the annotation file
(via ``-a``), a comma-separated list of alignment files (via ``-b``) and an output
directory where results files are stored (via ``-o``)::

    spladder build -o output_directory -b bam_file -a annotation_file

This will run SplAdder in its default configuration, which consists of the following steps:

- transform annotation into splicing graph representation
- generate an augmented splicing graph for each alignment file by inferring and
  adding the following elements:

    - insert intron retentions
    - insert cassette exons
    - insert new intron edges
- merge the augmented splicing graphs into a common splicing graph
- extract the following alternative splicing events:
    - exon skip
    - intron retention
    - alternative 3'/5' splice site
    - multiple exon skip
    - mutually exclusive exons
- quantify all alternative splicing events on each of the provided alignment
  files

Working with large datasets
^^^^^^^^^^^^^^^^^^^^^^^^^^^

SplAdder can be scaled to larger cohorts and has been run successfully on studies including as many
as 10000 samples. Working on such large datasets still needs some consideration and planning and
might require to call individual steps differently than how it is done for smaller sample sets.

We provide more detailed information in the section describing how to handle :ref:`Large cohorts
<spladder_cohorts>`.

