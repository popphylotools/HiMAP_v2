#!/usr/bin/env bash

# create main python environment as well as environments for 3rd party tools
conda env create -f environment.yml

conda deactivate
