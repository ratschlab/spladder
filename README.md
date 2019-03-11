What is SplAdder?
-----------------

This README describes the software SplAdder, short for Splicing Adder, a toolbox
for alternative splicing analysis based on RNA-Seq alignment data. Briefly, the
software takes a given annotation and RNA-Seq read alignments in standardized
formats, transforms the annotation into a splicing graph representation,
augments the splicing graph with additional information extracted from the read data,
extracts alternative splicing events from the graph and quantifies the events
based on the alignment data. The quantified events can then be used for
differential analysis.

Dependencies and Installation
-----------------------------
SplAdder relies on Python3 and requires only few standard packages that can be
installed using a Python package manager of your choice (e.g.,
[conda](http://conda.pydata.org/)):

* scipy (version >= 0.12 tested)
* pysam (version >= 0.7 required)
* h5py (version >= 2.2.0 tested)
* intervaltree (version >= 3.0.1 tested)
* matplotlib (for visualization mode; version >= 1.4.0 tested)
* statsmodels (for testing mode; version >= 0.9.0 tested)

SplAdder will not run without these packages installed. A complete list of dependencies can be also
found in [requirements.txt](requirements.txt)

Versions for Matlab and Python 2.7
----------------------------------
Previous versions of SplAdder were provided for both Matlab and Python. Since 2019, the
Matlab code is no longer provided as part of the SplAdder package. If you are
interested in the Matlab code, please download the initial release. 

If you are interested in previous versions of SplAdder capable of running under Python 2.7, please use
release 1.2.0. Please note, that the Python 2.7 code will be no longer maintained.

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
reading, please consider the [SplAdder Wiki](https://github.com/ratschlab/spladder/wiki).

After installation, the command `spladder` becomes available in your path.
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
