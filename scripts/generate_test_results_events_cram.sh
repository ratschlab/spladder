#!/bin/bash

set -e

testname="events"

datadir=tests/testcase_${testname}/data

cd ..
genome=$(pwd)/${datadir}/genome.fa
### merging splice graphs
outdir=tests/testcase_${testname}/results_merged_cram
mkdir -p $outdir
crams=""
bams=""
for i in $(seq 1 20)
do
    crams="$crams,${datadir}/align/testcase_${testname}_1_sample${i}.cram"
    bams="$bams,${datadir}/align/testcase_${testname}_1_sample${i}.bam"
done
#export REF_PATH=$genome
python -m spladder.spladder build -v -o ${outdir} -a ${datadir}/testcase_${testname}_spladder.gtf -b ${crams#,} --event-types exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons,mult_exon_skip --readlen 50 --output-conf-icgc --output-txt --output-txt-conf --output-gff3 --output-struc --output-struc-conf --output-bed --output-conf-bed --output-conf-tcga

cramsA=cramlistA.txt
for i in $(seq 1 10)
do
    echo "${datadir}/align/testcase_${testname}_1_sample${i}.cram" >> $cramsA
done
cramsB=cramlistB.txt
for i in $(seq 11 20)
do
    echo "align/testcase_${testname}_1_sample${i}.cram" >> $cramsB
done
python -m spladder.spladder test -o ${outdir} -v --diagnose-plots -f ps --readlen 50 --merge-strat merge_graphs --event-types exon_skip -a $cramsA -b $cramsB
rm cramlistA.txt cramlistB.txt
