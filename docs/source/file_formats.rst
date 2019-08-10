File formats
============

Input Formats
-------------

Annotation Files
^^^^^^^^^^^^^^^^
SplAdder accepts two different formats for annotation files: GTF and GFF. It will automatically
detect the format from the file name ending, so please make sure that your annotation ends with
either ``gtf`` or ``gff``.
Most sources for genome annotation provide their files in one of these two formats. If you would
like to generate your own annotation files, please follow the respective specifications:

GTF:
    `ensembl`_
GFF:
    `broad`_, `ucsc`_

Alignment Files
^^^^^^^^^^^^^^^

All alignment files are expected to be in BAM format, following the `SAM format specification`_. We
have successfully tested SplAdder with the following aligners:
- `STAR`_
- `PALMapper`_
- `TopHat`_

Output Formats
--------------
SplAdder produces a variety of different output files. Here we will mainly discuss files that are
aimed at the user and omit intermediate files that mainly necessary for internal processes of
SplAdder. Most of the latter will be stored in the ``spladder`` subdirectory in the output
directory.

After completing a SplAdder run, you will find several different output files in the output
directory. Following, we will describe each file type.

Annotation Files in GFF3 Format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These files have the general pattern
``merge_graphs_<event_type>_C<confidence_level>.confirmed.gff3`` and contain the events that have
been detected by SplAdder. Each event is shown as a mini gene consisting of two different isoforms.
If instance an exon skip would be described as::

    ##gff-version 3
    Chr1    exon_skip       gene    7616027 7616726 .       +       .       ID=exon_skip.1;GeneName="AT1G21690"
    Chr1    exon_skip       mRNA    7616027 7616726 .       +       .       ID=exon_skip.1_iso1;Parent=exon_skip.1;GeneName="AT1G21690"
    Chr1    exon_skip       exon    7616027 7616107 .       +       .       Parent=exon_skip.1_iso1
    Chr1    exon_skip       exon    7616603 7616726 .       +       .       Parent=exon_skip.1_iso1
    Chr1    exon_skip       mRNA    7616027 7616726 .       +       .       ID=exon_skip.1_iso2;Parent=exon_skip.1;GeneName="AT1G21690"
    Chr1    exon_skip       exon    7616027 7616107 .       +       .       Parent=exon_skip.1_iso2
    Chr1    exon_skip       exon    7616266 7616332 .       +       .       Parent=exon_skip.1_iso2
    Chr1    exon_skip       exon    7616603 7616726 .       +       .       Parent=exon_skip.1_iso2

For a definition of the different columns, please refer to one of the available GFF3 specifications
at `ensembl`_ or `ucsc`_. This file will allow you to display the events in a genome viewer such as
`UCSC Genome Browser`_, `IGV`_ or `GBrowse`_.

Event Files in HDF5 Format
^^^^^^^^^^^^^^^^^^^^^^^^^^

The event files contain all relevant event information and are stored in the hierarchical data
format `HDF5`_, allowing for efficient query, addition of data and interoperability between
different platforms and languages.
You can easily peek into the content of a hdf5 file::

    $> h5ls -r merge_graphs_exon_skip_C3.counts.hdf5

    /                        Group
    /conf_idx                Dataset {1}
    /event_counts            Dataset {4, 7, 2}
    /event_features          Group
    /event_features/alt_3prime Dataset {5}
    /event_features/alt_5prime Dataset {5}
    /event_features/exon_skip Dataset {7}
    /event_features/intron_retention Dataset {6}
    /event_features/mult_exon_skip Dataset {10}
    /event_features/mutex_exons Dataset {9}
    /event_pos               Dataset {2, 6}
    /gene_chr                Dataset {1}
    /gene_idx                Dataset {2}
    /gene_names              Dataset {1}
    /gene_pos                Dataset {1, 2}
    /gene_strand             Dataset {1}
    /strains                 Dataset {4}
    /verified                Dataset {2, 4}

This example is taken from the tutorial and lists the contents of the exon_skip event hdf5 file. The
tree that is shown looks a little bit like a file system tree and this is also the best analogy to
how the file is organized. Directories in the file system would correspond to groups in hdf5 and
files in file system to datasets in hdf5. Each group can contain more groups or datasets. 

The event hdf5 is structured as follows:

- **conf_idx**: 0-based index set, containing the index of the events that are confirmed in the provided samples
- **event_counts**: 3-dimensional matrix (S x F x E) containing counts for each of the E events, F features and S samples
- **event_features**: group that contains the description of the counted features per event type

    * **features alt3_prime / alt_5prime**: 
        + **valid**: contains a 1 if the event is valid and 0 otherwise
        + **exon_diff_cov**: mean coverage of the exonic segment that between the two alternative splice sites 
        + **exon_const_cov**: mean coverage of the remaining exonic segments in the event
        + **intron1_conf**: number of spliced alignments spanning the longer intron
        + **intron2_conf**: number of spliced alignments spanning the shorter intron
    * **features exon_skip**:
        + **valid**: contains a 1 if the event is valid and 0 otherwise
        + **exon_cov**: mean coverage of the cassette exon
        + **exon_pre_cov**: mean coverage of the left flanking exon (in genomic coordinates, ignoring strand)
        + **exon_aft_cov**: mean coverage of the right flanking exon (in genomic coordinates, ignoring strand)
        + **exon_pre_exon_conf**: number of spliced alignments spanning from left flanking to cassette exon
        + **exon_exon_aft_conf**: number of spliced alignments spanning from cassette to right flanking exon
        + **exon_pre_exon_aft_conf**: number of spliced alignments spanning from left flanking to right flanking exon
    * **features intron_retention**:
        + **valid**: contains a 1 if the event is valid and 0 otherwise
        + **intron_cov**: mean coverage of the retained intron
        + **exon1_cov**: mean coverage of the left flanking exon (in genomic coordinates, ignoring strand)
        + **exon2_cov**: mean coverage of the right flanking exon (in genomic coordinates, ignoring strand)
        + **intron_conf**: number of spliced alignments spanning the intron
        + **intron_cov_region**: fraction of positions in the intron that have a coverage > 0
    * **features mult_exon_skip**:
        + **valid**: contains a 1 if the event is valid and 0 otherwise
        + **exon_pre_cov**: mean coverage of the left flanking exon (in genomic coordinates, ignoring strand)
        + **exons_cov**: mean coverage over all skipped exons
        + **exon_aft_cov**: mean coverage of the right flanking exon (in genomic coordinates, ignoring strand)
        + **exon_pre_exon_conf**: number of spliced alignments spanning from left flanking to cassette exon
        + **exon_exon_aft_conf**: number of spliced alignments spanning from cassette to right flanking exon
        + **exon_pre_exon_aft_conf**: number of spliced alignments spanning from left flanking to right flanking exon
        + **sum_inner_exon_conf**: number of spliced alignments spanning any of the introns between neighboring skipped exons
        + **num_inner_exon**: number of skipped exons
        + **len_inner_exon**: cumulative length of skipped exons
    * **features mutex_exons**:
        + **valid**: contains a 1 if the event is valid and 0 otherwise
        + **exon_pre_cov**: mean coverage of the left flanking exon (in genomic coordinates, ignoring strand)
        + **exon1_cov**: mean coverage of the first skipped exon (first defined by genomic coordinates)
        + **exon2_cov**: mean coverage of the second skipped exon (second defined by genomic coordinates)
        + **exon_aft_cov**: mean coverage of the right flanking exon (in genomic coordinates, ignoring strand)
        + **exon_pre_exon1_conf**: number of spliced alignments spanning from left flanking to first exon
        + **exon_pre_exon2_conf**: number of spliced alignments spanning from left flanking to second exon
        + **exon1_exon_aft_conf**: number of spliced alignments spanning from first to right flanking exon
        + **exon2_exon_aft_conf**: number of spliced alignments spanning from second to right flanking exon
- **event_pos**: position of all event exons encoded as start,stop pairs for each event (events are rows, coordinates are columns)
- **gene_chr**: chromosome for each gene in the gene list
- **gene_idx**: index that maps each event to a gene in the gene list (0-based)
- **gene_names**: gene name for each gene in the gene list
- **gene_pos**: position of each gene in the gene list encoded as start,stop pair
- **gene_strand**: strand for each gene in the gene list
- **strains**: names of the samples counted
- **verified**: bool matrix over events X samples that is 1 if an event was verified in a sample and 0 otherwise

The naming of all these fields could be much more systematic but is currently kept the way it is to
not break compatibility with existing analysis pipelines. On a long term we plan to describe the
events and their counts in a more systematic way.

Event Files in TXT Format
^^^^^^^^^^^^^^^^^^^^^^^^^

Event files in txt format contain essentially the same information as the HDF5 files in a tab
delimited column format with one line per event and the following entries per line::

    1: chromosome of the event
    2: strand of the event
    3: unique event_id
    4: name of gene the event is located in
    5-5+n: start and stop coordinates of the event exons
    5+n and following: count values for each of the samples with the following layout (features are event type specific as defined above for HDF5 files:
        <sample1>:<feature1>
        <sample1>:<feature2>
        <sample1>:<feature3>
        ...
        <sample2>:<feature1>
        ...

Files in PICKLE Format
^^^^^^^^^^^^^^^^^^^^^^

These files are for internal usage only and can be ignored. 
        

.. _ensembl: http://www.ensembl.org/info/website/upload/gff.html
.. _broad: http://www.broadinstitute.org/annotation/argo/help/gff3.html
.. _ucsc: http://genome.ucsc.edu/FAQ/FAQformat.html#format3
.. _SAM format specification: https://samtools.github.io/hts-specs/SAMv1.pdf
.. _STAR: https://github.com/alexdobin/STAR
.. _PALMapper: http://www.raetschlab.org/suppl/palmapper/genomemapper-qpalma
.. _TopHat: https://ccb.jhu.edu/software/tophat/index.shtml
.. _UCSC Genome Browser: https://genome.ucsc.edu/cgi-bin/hgGateway
.. _IGV: http://www.broadinstitute.org/igv/
.. _GBrowse: http://gmod.org/wiki/GBrowse
.. _HDF5: https://www.hdfgroup.org/HDF5/

