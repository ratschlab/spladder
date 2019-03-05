#!/bin/bash

set -e

testname='testcase_events'
basedir=$(dirname $(pwd) )

test_data_dir=${basedir}/tests/${testname}/data
mkdir -p $test_data_dir

# Step 1and2: make changes on the test1.gtf file and reference sequence
# make some changes on the expression data and Create fq file
python create_event_testcases.py ${test_data_dir} ${testname} ${testname}_blocks.tsv ${testname}_abundances.tsv

### generate spladder GTF
grep -v t[2-9] ${test_data_dir}/${testname}.gtf > ${test_data_dir}/${testname}_spladder.gtf

# Step 3: Generating genome indexes for STAR
mkdir -p ${test_data_dir}/genome_idx
STAR --runThreadN 3 --runMode genomeGenerate --genomeDir ${test_data_dir}/genome_idx --genomeFastaFiles ${test_data_dir}/genome.fa --sjdbGTFfile ${test_data_dir}/${testname}.gtf --sjdbOverhang 25 --sjdbInsertSave Basic
rm Log.out

## Step 4: Running mapping jobs
mkdir ${test_data_dir}/align
for fname in ${test_data_dir}/${testname}*.fq
do
    fbase=$(basename $fname)
    fbase=${fbase%.fq}
    mkdir -p ${test_data_dir}/align/${fbase}

    STAR --runThreadN 3 --genomeDir ${test_data_dir}/genome_idx --readFilesIn ${fname} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${test_data_dir}/align/${fbase}/ --outSAMattributes NH HI NM MD AS XS && ln -s ${test_data_dir}/align/${fbase}/Aligned.sortedByCoord.out.bam ${test_data_dir}/align/${fbase}.bam
    samtools index ${test_data_dir}/align/${fbase}.bam
done
