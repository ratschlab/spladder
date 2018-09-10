#!/bin/bash

set -e

### set directories
basedir=examples
bams=${basedir}/NMD_WT1.tiny.bam,${basedir}/NMD_WT2.tiny.bam:${basedir}/NMD_DBL1.tiny.bam,${basedir}/NMD_DBL2.tiny.bam

### check that all data is there
if [ ! -f ${basedir}/result_tiny/merge_graphs_exon_skip_C3.counts.hdf5 ]
then
    echo "Example data not available - please run example_run.sh first"
    exit 1
fi

### run visualization
python spladder_viz.py -o ${basedir}/result_tiny -b $bams -c 3 -g AT1G21690 -T y -e exon_skip_1 -u y
python spladder_viz.py -o ${basedir}/result_tiny -b $bams -c 3 -T y -v y -u y

