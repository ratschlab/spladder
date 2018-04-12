Welcome to the Python implementation of SplAdder
------------------------------------------------

This README describes the Matlab version of the software SplAdder. Briefly, the
software takes a given annotation and RNA-Seq read alignments in standardized
formats, transforms the annotation into a splicing graph representation,
augments the splicing graph with additional information extracted from the read data,
extracts alternative splicing events from the graph and quantifies the events
based on the alignment data. The quantified events can then be used for
differential analysis.

Dependencies and Installation
-----------------------------

The Python version of SplAdder requires only few standard packages that are part
of most Python package managers (e.g., [conda](http://conda.pydata.org/)):
* scipy (version >= 0.12 tested)
* pysam (version >= 0.7 required)
* h5py (version >= 2.2.0 tested)
* intervaltree (version >= 2.1.0 tested)

SplAdder will not run without these packages installed. 

Authors
-------

Information on how to contact the authors can be found in the AUTHORS file.

License and Disclaimer
----------------------

All licensing information can be found in the COPYRIGHT file.

Documentation
-------------

This README provides a quick walk-through of a basic SplAdder run. For further
reading, please consider the [SplAdder Wiki](https://github.com/ratschlab/spladder/wiki).

The main SplAdder script *spladder.py* can be found at top level of this directory.
Invoking the executable without any parameters will print a description of the
command line interface to the screen.

In a basic call, SplAdder is invoked with three parameters: the annotation file
(via -a), a comma separated list of alignment files (via -b) and an output
directory where results files are stored (via -o). This will run SplAdder in its
default configuration, which consists of the following steps:

* transform annotation into splicing graph representation
* generate an augmented splicing graph for each alignment file by inferring and
  adding the following elements:
    - insert intron retentions
    - insert cassette exons
    - insert new intron edges
* merge the augmented splicing graphs into a common splicing graph
* extract the following alternative splicing events:
    - exon skip
    - intron retention
    - alternative 3'/5' splice site
    - multiple exon skip
* quantify all alternative splicing events on each of the provided alignment
  files
