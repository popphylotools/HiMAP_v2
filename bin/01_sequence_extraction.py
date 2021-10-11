#!/usr/bin/env python
# coding=utf-8
""""""
import os
import multiprocessing as mp
import configparser
import logging
from context import himap, project_dir
from himap import extract_sequences
log = logging.getLogger("01_sequence_extraction.py")

def main():
    #get annotation type dictionary from config file
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--configPath', help='configPath', default='./config.ini')
    args = parser.parse_args()

    # load parameters from config file
    config = configparser.ConfigParser()
    config.read(args.configPath)

    dtype_list = config["Annotation"]["source"].split()

    input_dict = dict()
    for dtype in dtype_list:
        path = config["Paths"]["working_dir"]
        for sp in [fn.rstrip(".gff") for fn in os.listdir(os.path.join(path, "gff"))]:
            gff_fn = os.path.join(path, "gff", sp + ".gff")
            fasta_fn = os.path.join(path, "fasta", sp + ".fasta")
            input_dict[sp] = (gff_fn, fasta_fn, path, sp)


    cpus = min(len(input_dict), mp.cpu_count())
    with mp.Pool(cpus) as p:
        p.starmap(himap.extract_sequences.driver, list(input_dict.values()))


# run main
try:
    exit(main())
except Exception as e:
    log.exception("Exception in main(): {}".format(e))
    exit(1)
