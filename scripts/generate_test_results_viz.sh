#!/bin/bash

set -e

testname="basic"

datadir=tests/testcase_${testname}/data

cd ..
### visualize splicing graphs
testcases="pos neg"
for testcase in $testcases
do
    outdir=tests/testcase_${testname}/results_merged_${testcase}
    mkdir -p $outdir
    python -m spladder.spladder viz -o ${outdir} -b ${datadir}/align/${testcase}_1.bam,${datadir}/align/${testcase}_2.bam,${datadir}/align/${testcase}_3.bam:${datadir}/align/${testcase}_4.bam,${datadir}/align/${testcase}_5.bam -L group1,group2 -f png
done

