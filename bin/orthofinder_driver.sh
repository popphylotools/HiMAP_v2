#!/usr/bin/env bash
eval "$(conda shell.bash hook)"
conda activate orthofinder

IN_DIR=${1}
OUT_DIR=${2}
CSV_DIR=${3}

orthofinder -f ${IN_DIR} -o ${OUT_DIR}

cp $(find . -name "Orthogroups.tsv") ${CSV_DIR}