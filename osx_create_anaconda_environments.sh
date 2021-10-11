#!/usr/bin/env bash

# create main python environment as well as environments for 3rd party tools
conda env create -f environment.yml
conda env create -f mafft_env.yml
conda env create -f tapir_env.yml
conda env create -f orthofinder_env.yml

# the tapir environment still needs hyphy2
source activate tapir

curl "http://s3.faircloth-lab.org/packages/hyphy2.osx.gz" -o "hyphy2.osx.gz"
gunzip hyphy2.*.gz
chmod 0700 hyphy2.*
mv hyphy2.* $(dirname $(which python))/hyphy2

source deactivate