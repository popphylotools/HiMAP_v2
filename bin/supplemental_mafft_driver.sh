#!/usr/bin/env bash

eval "$(conda shell.bash hook)"
conda activate mafft

IN_FILE=${1}
OUT_FILE=${2}

mafft --localpair --maxiterate 1000 ${IN_FILE} > ${OUT_FILE}
