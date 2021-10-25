#!/usr/bin/env python
# coding=utf-8
""""""
import os
import logging
import configparser
import numpy as np
import pandas as pd

log = logging.getLogger(os.path.basename(__file__))

# Define function factory for building filters
def fuzzy_factory(col_set, target_copies_count, min_isos_at_target, min_copies_count, max_copies_count):
    """
    sample_set            set of samples this test applies to

    target_copies_count   target number of copies per proteome (generally one)
    min_isos_at_target       min count of proteomes at target_copies_count

    min_copies_count            lower bound for proteomes off target count
    max_copies_count            upper bound for proteomes off target count
    """

    def fuzzy(rec):
        """
        rec    one row from counts dataframe
        """
        min_max_test = np.logical_and.reduce([min_copies_count <= rec[col] <= max_copies_count for col in col_set])
        count_on_target_test = sum([rec[col] == target_copies_count for col in col_set]) >= min_isos_at_target
        return min_max_test & count_on_target_test

    return fuzzy


def ortho_selection(tsv_dir, core_fuzzy, outgroup_fuzzy, supplemental_fuzzy, core_sp_list, outgroup_sp_list, supplemental_sp_list):
    from IPython.core.interactiveshell import InteractiveShell
    InteractiveShell.ast_node_interactivity = "all"

    # Load Data from Orthofinder in tsv format
    orthologs = pd.read_csv(os.path.join(tsv_dir, "Orthogroups.tsv"), delimiter="\t", low_memory=False)
    orthologs.rename(columns={'Orthogroup': 'ortho'}, inplace=True)
    orthologs.set_index('ortho', inplace=True)

    # Process Data to Have Counts Instead of a List of Gene Identifiers
    counts = orthologs.copy(deep=True)
    for col in counts.columns:
        counts.loc[counts[col].notnull(), col] = counts.loc[counts[col].notnull(), col].apply(
            lambda x: len(x.strip(",").split(",")))
        counts[col].fillna(0, inplace=True)


    # Apply each filter and print count of passing orthologs
    putative_single_copy_mask = np.logical_and.reduce([counts[col] == 1 for col in counts.columns])
    logging.info("putative_single_copy_mask passing: {}".format(sum(putative_single_copy_mask)))

    outgroup_mask = counts.apply(outgroup_fuzzy, axis=1)
    logging.info("outgroup_mask passing: {}".format(sum(outgroup_mask)))

    core_mask = counts.apply(core_fuzzy, axis=1)
    logging.info("core_mask passing: {}".format(sum(core_mask)))

    supplemental_mask = counts.apply(supplemental_fuzzy, axis=1)
    logging.info("supplemental_mask passing: {}".format(sum(supplemental_mask)))

    keep_mask = outgroup_mask & core_mask & supplemental_mask
    logging.info("keep_mask passing: {}".format(sum(keep_mask)))

    # ToDo: this substitution in the lambda seems potentially dangerous
    # Output Dataset
    orthologs[keep_mask].stack()[counts[keep_mask].stack() == 1].apply(
        lambda x: x.replace("__", "::")).unstack().to_csv(os.path.join(tsv_dir, "keep_orthos.csv"))


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--configPath', help='configPath', default='./config.ini')
    args = parser.parse_args()

    # load parameters from config file
    config = configparser.ConfigParser()
    config.read(args.configPath)

    # settings
    core_target_copies_count = config["Settings"].getint("core_target_copies_count")
    core_min_isos_at_target = config["Settings"]["core_min_isos_at_target"]
    core_min_copies_count = config["Settings"].getint("core_min_copies_count")
    core_max_copies_count = config["Settings"].getint("core_max_copies_count")

    outgroup_target_copies_count = config["Settings"].getint("outgroup_target_copies_count")
    outgroup_min_isos_at_target = config["Settings"].getint("outgroup_min_isos_at_target")
    outgroup_min_copies_count = config["Settings"].getint("outgroup_min_copies_count")
    outgroup_max_copies_count = config["Settings"].getint("outgroup_max_copies_count")

    supplemental_target_copies_count = config["Settings"].getint("supplemental_target_copies_count")
    supplemental_min_isos_at_target = config["Settings"].getint("supplemental_min_isos_at_target")
    supplemental_min_copies_count = config["Settings"].getint("supplemental_min_copies_count")
    supplemental_max_copies_count = config["Settings"].getint("supplemental_max_copies_count")

    # species lists
    core_sp_list = config["Species Groups"]["core_sp_list"].split()
    outgroup_sp_list = config["Species Groups"]["outgroup_sp_list"].split()
    supplemental_sp_list = config["Species Groups"]["supplemental_sp_list"].split()

    # input paths
    tsv_dir = config["Paths"]["tsv_dir"]

    # parse special keyword
    if core_min_isos_at_target == "all":
        core_min_isos_at_target = len(core_sp_list)
    else:
        core_min_isos_at_target = int(core_target_copies_count)

    # Define filters
    core_fuzzy = fuzzy_factory(col_set=core_sp_list,
                                   target_copies_count=core_target_copies_count,
                                   min_isos_at_target=core_min_isos_at_target,
                                   min_copies_count=core_min_copies_count,
                                   max_copies_count=core_max_copies_count)

    outgroup_fuzzy = fuzzy_factory(col_set=outgroup_sp_list,
                                   target_copies_count=outgroup_target_copies_count,
                                   min_isos_at_target=outgroup_min_isos_at_target,
                                   min_copies_count=outgroup_min_copies_count,
                                   max_copies_count=outgroup_max_copies_count)

    supplemental_fuzzy = fuzzy_factory(col_set=supplemental_sp_list,
                                       target_copies_count=supplemental_target_copies_count,
                                       min_isos_at_target=supplemental_min_isos_at_target,
                                       min_copies_count=supplemental_min_copies_count,
                                       max_copies_count=supplemental_max_copies_count)

    ortho_selection(tsv_dir, core_fuzzy, outgroup_fuzzy, supplemental_fuzzy, core_sp_list, outgroup_sp_list, supplemental_sp_list)
    print("Step_03 is complete\nResults are written to ", os.path.join(tsv_dir, "keep_orthos.csv"))
# run main
try:
    exit(main())
except Exception as e:
    log.exception("Exception in main(): {}".format(e))
    exit(1)
