#!/bin/bash

set -e

testname='testcase_events'
basedir=$(dirname $(pwd) )

test_data_dir=${basedir}/tests/${testname}/data
mkdir -p $test_data_dir

genome=${test_data_dir}/genome.fa

# create cram files
for fname in ${test_data_dir}/align/*.bam
do
    outbase=${fname%.bam}
    echo processing $outbase
    echo $genome
    #samtools view -h -F4 $fname | cramtools cram -O ${outbase}.cram -n --capture-tags 'NM' -L '*8' -R $genome
    samtools view -T $genome -C -o ${outbase}.cram $fname
    samtools index ${outbase}.cram 
done
