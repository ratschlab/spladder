#!/bin/bash

set -e

testname="basic"

datadir=tests/testcase_${testname}/data

cd ..
### merging splice graphs
testcases="pos neg"
for testcase in $testcases
do
    outdir=tests/testcase_${testname}/results_merged_${testcase}
    mkdir -p $outdir
    python -m spladder.spladder -o ${outdir} -a ${datadir}/annotation_${testcase}.gtf -b ${datadir}/align/${testcase}_1.bam,${datadir}/align/${testcase}_2.bam,${datadir}/align/${testcase}_3.bam,${datadir}/align/${testcase}_4.bam,${datadir}/align/${testcase}_5.bam -T y -t exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons,mult_exon_skip -n 15
done

### merging splice graphs
testcases="pos neg"
for testcase in $testcases
do
    outdir=tests/testcase_${testname}/results_single_${testcase}
    mkdir -p $outdir
    python -m spladder.spladder -o ${outdir} -a ${datadir}/annotation_${testcase}.gtf -b ${datadir}/align/${testcase}_1.bam -T y -t exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons,mult_exon_skip -n 15 --merge_strat single
done

