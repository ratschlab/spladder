What is SplAdder?
-----------------

This README describes the software SplAdder, short for Splicing Adder, a toolbox
for alternative splicing analysis based on RNA-Seq alignment data. Briefly, the
software takes a given annotation and RNA-Seq read alignments in standardized
formats, transforms the annotation into a splicing graph representation,
augments the splicing graph with additional information extracted from the read data,
extracts alternative splicing events from the graph and quantifies the events
based on the alignment data, The quantified events can then be used for
differential analysis.

Implementation and Dependencies
-------------------------------

Currently, there are two implementations available, which produce the same
results on the same input data. The development of SplAdder was started in
Matlab and has been moved to Python in the past two years, mostly for the reason
of having fewer dependencies. Both implementations are available in this
repository und the respective directories.

Note, that the current Matlab implementation makes extensive use of the HDF5
storage format and uses low level access to the files. Unfortunately, this
functionality is not available in Octave, therefor requiring a working Matlab
installation. In case the compatibility issues between Matlab and Octave
regarding low level HDF5 should get resolved in the near future, SplAdder might
be able to run also under Octave.

However, future development is planned to focus on the Python branch. The Python
version of the code has been developed in python 2.7 and requires the following
packages:
* scipy (version >= 0.12 tested)
* pysam (version >= 0.7 required)
* h5py (version >= 2.2.0 tested)
* matplotlib (for the visualization code only; version >= 1.4.0 tested)

Installation
------------

Please see the respective INSTALL files in the matlab or python directories for
details on how to install SplAdder on your system. 

Authors
-------

Information on how to contact the authors can be found in the AUTHORS file.

License and Disclaimer
----------------------

All licensing information can be found in the COPYRIGHT file.

Documentation
-------------

This README provides a quick walk-through of a basic SplAdder run. For further
reading, please consider the file DOCUMENTATION in the matlab or python
directory. The following description is generic for both implementations.

Invoking SplAdder without any parameters will print a description of the
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
    - mutually exclusive exons
* quantify all alternative splicing events on each of the provided alignment
  files
