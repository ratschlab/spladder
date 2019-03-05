#!/bin/bash

set -e

testname="events"

datadir=tests/testcase_${testname}/data

cd ..
### merging splice graphs
outdir=tests/testcase_${testname}/results_merged
mkdir -p $outdir
bams=$(ls -1 ${datadir}/align/*.bam | tr '\n' ',' | head -c-1)
python -m spladder.spladder build -v -o ${outdir} -a ${datadir}/testcase_${testname}_spladder.gtf -b ${bams} --event-types exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons,mult_exon_skip --readlen 50

bamsA=bamlistA.txt
for i in $(seq 1 10)
do
    echo "${datadir}/align/testcase_${testname}_1_sample${i}.bam" >> $bamsA
done
bamsB=bamlistB.txt
for i in $(seq 11 20)
do
    echo "${datadir}/align/testcase_${testname}_1_sample${i}.bam" >> $bamsB
done
python -m spladder.spladder test -o ${outdir} -v --diagnose_plots --readlen 50 --merge-strat merge_graphs --event-types exon_skip -a $bamsA -b $bamsB
