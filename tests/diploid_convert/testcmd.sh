#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.fasta"
EXPECTED_OUTPUTS[0]="${SELF}/expected.fasta"
CURRENT_OUTPUTS[1]="${SELF}/output.mfasta"
EXPECTED_OUTPUTS[1]="${SELF}/expected.mfasta"
CURRENT_OUTPUTS[2]="${SELF}/output.vcf"
EXPECTED_OUTPUTS[2]="${SELF}/expected.vcf"

TESTCMD="\
    ${BRESEQ} \
        CONVERT-REFERENCE \
        -f MULTIPLOID \
        -o ${SELF}/output \
        ${DATADIR}/REL606/REL606.fragment.fna \
        ${SELF}/input.vcf \
    "

do_test $1 ${SELF}
