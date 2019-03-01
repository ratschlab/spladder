#!/bin/bash

set -e

testname='testcase_basic'

test_data_dir=../tests/${testname}/data
mkdir -p $test_data_dir
cp ../tests/templates/annotation_pos.gtf ${test_data_dir}


# Step 1and2: make changes on the test1.gtf file and reference sequence
# make some changes on the expression data and Create fq file
python create_bam_file.py ${test_data_dir} ${testname} 5

# Step 3: Generating genome indexes for STAR
mkdir -p ${test_data_dir}/genome/pos
mkdir -p ${test_data_dir}/genome/neg
STAR --runThreadN 3 --runMode genomeGenerate --genomeDir ${test_data_dir}/genome/pos --genomeFastaFiles ${test_data_dir}/genome_pos.fa --sjdbGTFfile ${test_data_dir}/annotation_pos.gtf --sjdbOverhang 14 --sjdbInsertSave Basic
STAR --runThreadN 3 --runMode genomeGenerate --genomeDir ${test_data_dir}/genome/neg --genomeFastaFiles ${test_data_dir}/genome_neg.fa --sjdbGTFfile ${test_data_dir}/annotation_neg.gtf --sjdbOverhang 14 --sjdbInsertSave Basic

## Step 4: Running mapping jobs
mkdir ${test_data_dir}/align
for i in $(seq 1 5)
do
    mkdir -p ${test_data_dir}/align/pos_${i}
    mkdir -p ${test_data_dir}/align/neg_${i}

    STAR --runThreadN 3 --genomeDir ${test_data_dir}/genome/pos --readFilesIn ${test_data_dir}/${testname}_${i}.fq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${test_data_dir}/align/pos_${i}/ --outSAMattributes NH HI NM MD AS XS && ln -s pos_${i}/Aligned.sortedByCoord.out.bam ${test_data_dir}/align/pos_${i}.bam
    samtools index ${test_data_dir}/align/pos_${i}.bam

    STAR --runThreadN 3 --genomeDir ${test_data_dir}/genome/neg --readFilesIn ${test_data_dir}/${testname}_${i}.fq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${test_data_dir}/align/neg_${i}/ --outSAMattributes NH HI NM MD AS XS && ln -s neg_${i}/Aligned.sortedByCoord.out.bam ${test_data_dir}/align/neg_${i}.bam
    samtools index ${test_data_dir}/align/neg_${i}.bam
done
