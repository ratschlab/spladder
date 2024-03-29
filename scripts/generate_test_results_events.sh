#!/bin/bash

set -e

testname="events"

datadir=tests/testcase_${testname}/data

cd ..
### merging splice graphs
outdir=tests/testcase_${testname}/results_merged
mkdir -p $outdir
bams=""
for i in $(seq 1 20)
do
    bams="$bams,${datadir}/align/testcase_${testname}_1_sample${i}.bam"
done
python -m spladder.spladder build -v -o ${outdir} -a ${datadir}/testcase_${testname}_spladder.gtf -b ${bams#,} --event-types exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons,mult_exon_skip --readlen 50 --output-conf-icgc --output-txt --output-txt-conf --output-gff3 --output-struc --output-struc-conf --output-bed --output-conf-bed --output-conf-tcga #--use-anno-support

bamsA=bamlistA.txt
rm -f $bamsA
for i in $(seq 1 10)
do
    echo "${datadir}/align/testcase_${testname}_1_sample${i}.bam" >> $bamsA
done
bamsB=bamlistB.txt
rm -f $bamsB
for i in $(seq 11 20)
do
    echo "align/testcase_${testname}_1_sample${i}.bam" >> $bamsB
done
python -m spladder.spladder test -o ${outdir} -v --diagnose-plots -f pdf --readlen 50 --merge-strat merge_graphs --event-types exon_skip -a $bamsA -b $bamsB --dpsi 0
rm bamlistA.txt bamlistB.txt
