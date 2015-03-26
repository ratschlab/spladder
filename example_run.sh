#!/bin/bash

set -e

### set dirs and variables
basedir=examples
anno=${basedir}/TAIR10_GFF3_genes.tiny.gff
bams=${basedir}/NMD_WT1.tiny.bam,${basedir}/NMD_WT2.tiny.bam,${basedir}/NMD_DBL1.tiny.bam,${basedir}/NMD_DBL2.tiny.bam
outdir=${basedir}/result_tiny

### check if all data is there
if [ ! -f ${basedir}/NMD_WT1.tiny.bam -o ! -f ${basedir}/NMD_WT2.tiny.bam  -o ! -f ${basedir}/NMD_DBL1.tiny.bam ]
then
    echo "Example data has not been downloaded yet. Do you want to download the data now? [yes|no]"
    read DNLD
    if [ "$DNLD" == yes ]
    then
        cd ${basedir}
        [[ ! -f NMD_WT1.tiny.bam ]] && wget http://ftp.ratschlab.org/user/akahles/data/spladder/examples_small/NMD_WT1.tiny.bam
        [[ ! -f NMD_WT1.tiny.bam.bai ]] && wget http://ftp.ratschlab.org/user/akahles/data/spladder/examples_small/NMD_WT1.tiny.bam.bai
        [[ ! -f NMD_WT2.tiny.bam ]] && wget http://ftp.ratschlab.org/user/akahles/data/spladder/examples_small/NMD_WT2.tiny.bam
        [[ ! -f NMD_WT2.tiny.bam.bai ]] && wget http://ftp.ratschlab.org/user/akahles/data/spladder/examples_small/NMD_WT2.tiny.bam.bai
        [[ ! -f NMD_DBL1.tiny.bam ]] && wget http://ftp.ratschlab.org/user/akahles/data/spladder/examples_small/NMD_DBL1.tiny.bam
        [[ ! -f NMD_DBL1.tiny.bam.bai ]] && wget http://ftp.ratschlab.org/user/akahles/data/spladder/examples_small/NMD_DBL1.tiny.bam.bai
        [[ ! -f NMD_DBL2.tiny.bam ]] && wget http://ftp.ratschlab.org/user/akahles/data/spladder/examples_small/NMD_DBL2.tiny.bam
        [[ ! -f NMD_DBL2.tiny.bam.bai ]] && wget http://ftp.ratschlab.org/user/akahles/data/spladder/examples_small/NMD_DBL2.tiny.bam.bai
        cd ..
    else
        echo "Nothing to process. Exiting."
        exit 1
    fi
fi

### create output directory
rm -rf $outdir
mkdir -p $outdir

### run SplAdder.
./bin/spladder.sh -b $bams -o $outdir -a $anno -v y -M merge_graphs -t exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons,mult_exon_skip -c 3 -p n
