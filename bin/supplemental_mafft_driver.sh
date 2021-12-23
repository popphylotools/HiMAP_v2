#!/usr/bin/env bash

IN_FILE=${1}
OUT_FILE=${2}

mafft --localpair --maxiterate 1000 ${IN_FILE} > ${OUT_FILE}
