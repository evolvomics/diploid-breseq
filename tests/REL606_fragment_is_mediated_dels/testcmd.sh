#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output/evidence/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/REL606/REL606.fragment.gbk"

TESTCMD=" \
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -b 0 \
    -o ${SELF} \
    ${REFERENCE_ARG} \
    ${DATADIR}/REL606/REL606.fragment.3.fastq.gz
	"

do_test $1 ${SELF}
