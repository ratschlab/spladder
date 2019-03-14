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
------------------------------------

SplAdder relies on Python3 and requires only few standard packages. A complete list of dependencies can be found in [requirements.txt](requirements.txt)

We recommend using a python package manager such as anaconda to organize dependencies.

* Installation via `pip`: `pip install spladder`
* Installation from source:
   1. clone this repository
   2. Within the SplAdder root directory, run `make install`


Running Spladder
--------------------

To start SplAdder and see a list of all available options, you can simply type:
```
spladder
```


Test with Example Data
---------------------------

If you installed SplAdder from source, you can test the installation by invoking ```make test``` in the SplAdder root directory.


Versions for Matlab and Python 2.7
----------------------------------
Previous versions of SplAdder were provided for both Matlab and Python. Since 2019, the
Matlab code is no longer provided as part of the SplAdder package. If you are
interested in the Matlab code, please download the initial release. 

If you are interested in previous versions of SplAdder capable of running under Python 2.7, please use
release 1.2.1. Please note, that the Python 2.7 code will be no longer maintained.


Authors
-------

Information on how to contact the authors can be found in the AUTHORS file.

License and Disclaimer
----------------------

All licensing information can be found in the COPYRIGHT file.

Documentation
-------------

This README provides only a high-level overview of a basic SplAdder run. For further
reading, please consider the [SplAdder Wiki](https://github.com/ratschlab/spladder/wiki).

After installation, the command `spladder` becomes available in your path.
Invoking SplAdder without any parameters will print a description of the
command line interface to the screen.

In a basic call in build-mode (`spladder build`), SplAdder requires at least three parameters: the annotation file
(via `-a`), a comma-separated list of alignment files (via `-b`) and an output
directory where results files are stored (via `-o`). This will run SplAdder in its
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
