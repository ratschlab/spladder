#/bin/bash
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2013 Andre Kahles, Jonas Behr, Gunnar Raetsch
# Copyright (C) 2009-2011 Max Planck Society
# Copyright (C) 2012-2013 Memorial Sloan-Kettering Cancer Center
#
# SplAdder wrapper script to start the interpreter with the correct list of arguments

set -e


function usage () {
    echo "
    
    Usage: SplAdder [-OPTION VALUE] 

    Options (default values in [...]):

    MANDATORY:

    -b  FILE1,FILE2,...     alignment files in BAM format (comma separated list)
    -o  DIR                 output directory
    -a  FILE                annotation file name (annotation in *.mat format)

    OPTIONAL:
    -l  FILE                log file name [stdout]
    -u  FILE                file with user settings [-]
    -c  INT                 confidence level (0 lowest to 3 highest) [3]
    -I  INT                 number of iterations to insert new introns
                            into the graph [5]
    -M  <STRAT>             merge strategy, where <STRAT> is one on:
                            merge_bams, merge_graphs, merge_all [merge_graphs]
    -R  R1,R2,...           replicate structure of files (same number as
                            alignment files) [all R1 - no replicated]
    -L  STRING              label for current experiment [-]
    -S  STRING              reference strain [-]
    -d  y|n                 use debug mode [n]
    -p  y|n                 use rproc [n]
    -V  y|n                 validate splice graph [n]
    -v  y|n                 use verbose output mode [n]
    -A  y|n                 curate alt prime events [y]
    -x  y|n                 input alignments share the same genome [y]
    -i  y|n                 insert intron retentions [y]
    -e  y|n                 insert cassette exons [y]
    -E  y|n                 insert new intron edges [y]
    -r  y|n                 remove short exons [n]
    -s  y|n                 re-infer splice graph [n]
    "

    exit 0
}

[[ -z "$1" ]] && usage

### mandatory parameters
S_BAM_FNAME=""
S_ANNO_FNAME=""
S_OUT_DIRNAME=""

### parse parameters from command lines
while getopts "b:o:l:a:u:c:I:M:R:L:S:d:p:V:v:A:x:i:e:E:r:s:h" opt
do
    case $opt in
    b ) S_BAM_FNAME="$OPTARG" ;;
    o ) S_OUT_DIRNAME="$OPTARG" ;;
    l ) S_LOG_FNAME="$OPTARG" ;;
    a ) S_ANNO_FNAME="$OPTARG" ;;
    u ) S_USER_FNAME="$OPTARG" ;;
    c ) I_CONFIDENCE="$OPTARG" ;;
    I ) I_INSERT_INTRON_ITER="$OPTARG" ;;
    M ) S_MERGE_STRATEGY="$OPTARG" ;;
    R ) S_REPLICATE_IDX="$OPTARG" ;;
    L ) S_EXPERIMENT_LABEL="$OPTARG" ;;
    S ) S_REFERENCE_STRAIN="$OPTARG" ;;
    d ) F_DEBUG="$OPTARG" ;;
    p ) F_RPROC="$OPTARG" ;;
    V ) F_VALIDATE_SG="$OPTARG" ;;
    v ) F_VERBOSE="$OPTARG" ;;
    A ) F_CURATE_ALTPRIME="$OPTARG" ;;
    x ) F_SHARE_GENESTRUCT="$OPTARG" ;;
    i ) F_INSERT_IR="$OPTARG" ;;
    e ) F_INSERT_CE="$OPTARG" ;;
    E ) F_INSERT_IE="$OPTARG" ;;
    r ) F_REMOVE_SE="$OPTARG" ;;
    s ) F_INFER_SG="$OPTARG" ;;
    h ) usage ;;
    \?) echo -e "UNKNOWN PARAMETER: $opt\n\n"; usage ;;
    esac
done

### check for mandatory parameters
[[ -z "$S_BAM_FNAME" ]] && echo -e "\nAlignment file name is mandatory!\n" && usage && exit 1
[[ -z "$S_ANNO_FNAME" ]] && echo -e "\nAnnotation file name is mandatory!\n" && usage && exit 1
[[ -z "$S_OUT_DIRNAME" ]] && echo -e "\nOutput directory name is mandatory!\n" && usage && exit 1

### assemble parameter string
PARAMS=""
for opt in S_BAM_FNAME S_OUT_DIRNAME S_LOG_FNAME S_ANNO_FNAME S_USER_FNAME I_CONFIDENCE I_INSERT_INTRON_ITER F_DEBUG F_VERBOSE F_INSERT_IR F_INSERT_CE F_INSERT_IE F_REMOVE_SE F_INFER_SG F_VALIDATE_SG S_MERGE_STRATEGY F_SHARE_GENESTRUCT S_REPLICATE_IDX S_EXPERIMENT_LABEL S_REFERENCE_STRAIN F_CURATE_ALTPRIME F_RPROC 
do
    if [ ! -z "$opt" ]
    then
        eval "PARAMS=\"$PARAMS${opt}:\${$opt};\""
    fi
done

### Index Bam files
echo "Indexing BAM files"
SAMPLE_LIST=()
OLD_IFS=$IFS
IFS=',' 
for BAM_FILE in "${S_BAM_FNAME}"
do
    CURR_BAMFILE=$BAM_FILE
    if [ ! -f ${CURR_BAMFILE}.bai ]
    then
        echo "Indexing $CURR_BAMFILE"
        ${SPLADDER_SAMTOOLS_BIN_DIR} index $BAM_FILE
    else
        echo "$CURR_BAMFILE already indexed"
    fi
done
IFS=$OLD_IFS
