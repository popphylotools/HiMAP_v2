#!/usr/bin/env bash

IN_DIR=${1}
OUT_DIR=${2}
CSV_DIR=${3}

rm -r ${CSV_DIR}/*

orthofinder -f ${IN_DIR} -o ${OUT_DIR}

cp $(find . -name "Orthogroups.tsv") ${CSV_DIR}
