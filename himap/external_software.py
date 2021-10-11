#!/usr/bin/env python
# coding=utf-8
""""""

from multiprocessing.pool import ThreadPool
import multiprocessing as mp
import subprocess
import os


def core_mafft_driver_file(file):
    p = subprocess.Popen(["./bin/core_mafft_driver.sh", file, file + ".aln"],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return out, err


def core_mafft_driver_path(path):
    # remove old alignments
    rm_files = [os.path.join(path, file) for file in os.listdir(path) if ".aln" in file]
    for file in rm_files:
        os.remove(file)

    # call maft on each fasta
    files = [os.path.join(path, file) for file in os.listdir(path) if ".fasta" in file]
    pool = ThreadPool(mp.cpu_count())
    pool.map(core_mafft_driver_file, files)
    pool.close()
    pool.join()


def supplemental_mafft_driver_file(file):
    p = subprocess.Popen(["./bin/supplemental_mafft_driver.sh", file, file + ".aln"],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return out, err


def supplemental_mafft_driver_path(path):
    # remove old alignments
    rm_files = [os.path.join(path, file) for file in os.listdir(path) if ".aln" in file]
    for file in rm_files:
        os.remove(file)

    # call maft on each fasta
    files = [os.path.join(path, file) for file in os.listdir(path) if ".fasta" in file]
    pool = ThreadPool(mp.cpu_count())
    pool.map(supplemental_mafft_driver_file, files)
    pool.close()
    pool.join()


def orthofind_driver_dir(path_in, path_out, tsv_dir):
    p = subprocess.Popen(["./bin/orthofinder_driver.sh", path_in, path_out, tsv_dir],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return out, err