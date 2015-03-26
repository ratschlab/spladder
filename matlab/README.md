What is SplAdder?
-----------------

This README describes the software SplAdder, short for splicing adder, a toolbox
for alternative splicing analysis based on RNA-Seq alignment data. Briefly, the
software takes a given annotation and RNA-Seq read alignments in standardized
formats, transforms the annotation into a splicing graph representation,
augments the splicing graph with additional information extracted from the read data,
extracts alternative splicing events from the graph and quantifies the events
based on the alignment data, The quantified events can then be used for
differential analysis.

Dependencies
------------

SplAdder is written in Matlab code and requires Matlab to run. As it makes use
of the HDF5 storage format, currently difficulties arise from the usage of
Octave. In case the compatibility issues between Matlab and Octave regarding
HDF5 should get resolved, SplAdder might be able to run also under Octave.

Installation
------------

Please see the file INSTALL for details on how to install SplAdder on your
system.

Authors
-------

Information on how to contact the authors can be found in the AUTHORS file.

License and Disclaimer
----------------------

All licensing information can be found in the COPYRIGHT file.

Documentation
-------------

This README provides a quick walk-through of a basic SplAdder run. For further
reading, please consider the file DOCUMENTATION.

The SplAdder executable *spladder.sh* can be found in the bin directory.
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
