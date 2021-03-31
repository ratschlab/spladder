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
python -m spladder.spladder build -v -o ${outdir} -a ${datadir}/testcase_${testname}_spladder.gtf -b ${crams#,} --event-types exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons,mult_exon_skip --readlen 50 --output-conf-icgc --output-txt --output-txt-conf --output-gff3 --output-struc --output-struc-conf --output-bed --output-conf-bed --output-conf-tcga --reference ${genome}

cramsA=cramlistA.txt
rm -f $cramsA
for i in 10 8 2 1 7 6 5 3 9 4
do
    echo "${datadir}/align/testcase_${testname}_1_sample${i}.cram" >> $cramsA
done
cramsB=cramlistB.txt
rm -f $cramsB
for i in 20 13 17 11 12 19 15 14 16 18
do
    echo "align/testcase_${testname}_1_sample${i}.cram" >> $cramsB
done
python -m spladder.spladder test -o ${outdir} -v --diagnose-plots -f pdf --readlen 50 --merge-strat merge_graphs --event-types exon_skip -a $cramsA -b $cramsB --dpsi 0
rm cramlistA.txt cramlistB.txt
