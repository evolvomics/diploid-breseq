#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output/evidence/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"

TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    --bowtie2-scoring \"--ma 1 --mp 1 --np 1 --rdg 1,2 --rfg 1,2\" \
    --bowtie2-stage1 '--local --score-min L,1,0.4 -k 20' \
    --bowtie2-stage2 '' \
    --bowtie2-junction '--local --score-min L,1,0.30  -k 20' \
    ${REFERENCE_ARG} \
    ${DATADIR}/lambda/lambda_mixed_population.fastq \
    "

do_test $1 ${SELF}
