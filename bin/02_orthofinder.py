#!/usr/bin/env python
# coding=utf-8
""""""
import os
import logging
import configparser
from context import himap, project_dir
from himap import external_software
log = logging.getLogger(os.path.basename(__file__))


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--configPath', help='configPath', default='./config.ini')
    args = parser.parse_args()

    # read config file
    config = configparser.ConfigParser()
    config.read(args.configPath)

    # input directory
    pep_fasta_dir = config["Paths"]["pep_fasta_dir"]
    
    #output directory
    tsv_dir = config["Paths"]["tsv_dir"]
    out_dir = os.path.join(tsv_dir, "Orthofinder")
    
    himap.external_software.orthofind_driver_dir(pep_fasta_dir, out_dir, tsv_dir)
    
    print("Step_02 is complete\nOrthofinder results are written to ", out_dir)


# run main
try:
    exit(main())
except Exception as e:
    log.exception("Exception in main(): {}".format(e))
    exit(1)
