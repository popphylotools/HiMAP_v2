#!/usr/bin/env python
# coding=utf-8
""""""
import os
import pandas as pd
import re
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import gffutils
import json
from Bio.Seq import Seq
from pyfaidx import Fasta
import logging
import configparser
from context import himap, project_dir, n_gap

from himap import external_software
from himap import consensus
from himap import parse_msa

from multiprocessing.pool import ThreadPool
import multiprocessing as mp


# set up logging
# logging configuration
# noinspection PyMissingOrEmptyDocstring
class OneLineExceptionFormatter(logging.Formatter):
    def formatException(self, exc_info):
        result = super().formatException(exc_info)
        return repr(result)

    def format(self, record):
        """Format error message to fit on a single line.
        """
        result = super().format(record)
        if record.exc_text:
            result = result.replace(r"\n", r"/n")
        return result


handler = logging.StreamHandler()
formatter = OneLineExceptionFormatter(logging.BASIC_FORMAT)
handler.setFormatter(formatter)
root = logging.getLogger()
root.setLevel(os.environ.get("LOGLEVEL", "DEBUG"))
root.addHandler(handler)
log = logging.getLogger(os.path.basename(__file__))


# define functions to parse coordinates of cds's from concatinated aligned fasta w/ n's and -'s
def findBreakpoints(seq, n_gap):
    n_count = n_gap
    breakpoints = []
    loc = 0
    regex = re.compile(r"n+[-+n]*")
    while True:
        # loc = seq.find(nnn, loc)
        match = regex.search(seq, loc)
        if not match:
            break
        if len(match.group().replace('-', '')) >= n_count:
            breakpoints.append(match.span())
        loc = match.end()
    return breakpoints


def findExonCoords(seq, n_gap):
    breakpoints = findBreakpoints(seq, n_gap)
    length = len(seq)

    if len(breakpoints) == 0:
        return [(0, length)]

    if len(breakpoints) == 1:
        bp = breakpoints[0]
        return [(0, bp[0]), (bp[1], length)]

    elif len(breakpoints) > 0:
        exonCoords = [(0, breakpoints[0][0])]

        for i in range(len(breakpoints) + 1)[1:-1]:  # all intermediate exons
            ex_start = breakpoints[i - 1][1]
            ex_end = breakpoints[i][0]
            exonCoords.append((ex_start, ex_end))

        exonCoords.append((breakpoints[-1][1], length))  # last exon
        return exonCoords

def gapPercent(seq):
    seq = str(seq)
    gappedLen = len(seq)
    gapCount = seq.count('-')
    return (100.0 * gapCount) / gappedLen


def longestGap(seq):
    seq = str(seq)
    gap_regex = re.compile(r"-+")
    gap_list = gap_regex.findall(seq)
    if gap_list:
        return sorted([len(gap) for gap in gap_list], reverse=True)[0]
    else:
        return 0


def get_frame(cds):
    if cds.frame == ".":
        return 0
    else:
        return int(cds.frame)


def prep_output(sp, iso_id, db, fasta, n_gap):
    iso = db[iso_id]

    cds_list = list(db.children(iso, order_by="start", featuretype="CDS"))
    n_seq_list = [Seq(c.sequence(fasta, use_strand=False)) for c in cds_list]

    if iso.strand == "-":
        n_seq_list = [n_seq.reverse_complement() for n_seq in n_seq_list[::-1]]
        frame = get_frame(cds_list[-1])
    else:
        frame = get_frame(cds_list[0])
    sep = "N" * n_gap
    n_seq = Seq(sep.join([str(c) for c in n_seq_list]))

    return SeqRecord(n_seq[frame:], id=sp, description="")

#melt ortho db to assign identifiers 

def maft_core_prep(gff_db_dict, fasta_dict, working_dir, tsv_dir, core_align_dir, add_outgroup_in_core_align,
                   core_sp_set, outgroup_sp_set, n_gap):

    orthoMelt = pd.read_csv(os.path.join(tsv_dir, "keep_orthos.csv")).melt(id_vars="ortho", var_name="species", value_name="iso_id")
    orthoMelt = orthoMelt.loc[orthoMelt["iso_id"].notnull()]

    orthoMelt["core"] = orthoMelt.species.apply(lambda sp: sp in core_sp_set)
    orthoMelt["outgroup"] = orthoMelt.species.apply(lambda sp: sp in outgroup_sp_set)

    orthoMelt["cds_cat"] = orthoMelt.apply(
        lambda rec: prep_output(rec["species"], rec["iso_id"], gff_db_dict[rec["species"]],
                                fasta_dict[rec["species"]], n_gap).format("fasta"), axis=1)

    for ortho in sorted(orthoMelt["ortho"].unique()):
        with open(os.path.join(core_align_dir, ortho + ".fasta"), "w") as f:
            if add_outgroup_in_core_align:
                ortho_df = orthoMelt[(orthoMelt["ortho"] == ortho) & (orthoMelt["core"] | orthoMelt["outgroup"])]
            else:
                ortho_df = orthoMelt[(orthoMelt["ortho"] == ortho) & (orthoMelt["core"])]
            ortho_df.index = ortho_df["species"]
            for sp in sorted(ortho_df["species"].unique()):
                f.write(ortho_df.loc[sp, "cds_cat"])

    # create parent_groups
    groups = dict()
    for ortho in sorted(orthoMelt["ortho"].unique()):
        groups[ortho] = dict()
        ortho_df = orthoMelt[(orthoMelt["ortho"] == ortho)]
        ortho_df.index = ortho_df["species"]
        for sp in sorted(ortho_df["species"].unique()):
            groups[ortho][sp] = ortho_df.loc[sp, "iso_id"]

    # output parent_groups to groups.json
    filename = os.path.join(working_dir, "groups.json")
    with open(filename, 'w') as f:
        json.dump(groups, f)


#######################################
# core post and supplemental Prep     #
#######################################


def maft_core_post_and_supp_prep(gff_db_dict, fasta_dict, working_dir, core_align_dir, sup_align_dir,
                                 add_outgroup_in_sup_align, core_sp_list, outgroup_sp_list, n_gap, min_exon_length, max_gap_percent, max_gap_length):

    # import ortholog groups
    with open(os.path.join(working_dir, "groups.json"), 'r') as f:
        parent_groups = json.load(f)

    # create handles for all .fasta files in aligned core fasta directory
    aligned_fasta_fn = {name.split('.fasta')[0]: os.path.join(core_align_dir, name) for name in
                        os.listdir(core_align_dir) if
                        ((".fasta.aln" in name) and (".fasta.aln.fai" not in name))}

    # read and parse fasta files for each species
    aligned_fasta = {}
    for ortho in aligned_fasta_fn.keys():
        aligned_fasta[ortho] = {seq_record.id: seq_record
                                for seq_record in SeqIO.parse(aligned_fasta_fn[ortho],
                                                              "fasta")}

    filtered_universal_core_exon_spans = {}
    sp_core_exon_spans = {}
    for ortho, seq_dict in aligned_fasta.items():
        sp_exon_list, filtered_universal_exon_spans = himap.parse_msa.universal_exon_spans(seq_dict, n_gap, min_exon_length, max_gap_percent, max_gap_length)
        if len(filtered_universal_exon_spans) > 0:
            sp_core_exon_spans[ortho] = sp_exon_list
            filtered_universal_core_exon_spans[ortho] = filtered_universal_exon_spans
#        else:
#            logging.debug("ortho:{} - no universal filtered core exons, excluding ortho".format(ortho))

    fasta_prep = {}
    for ortho in filtered_universal_core_exon_spans:
        fasta_prep[ortho] = {}
        cds_list = {}
        for coord in filtered_universal_core_exon_spans[ortho]:
            fasta_prep[ortho][coord] = {}
            for sp in core_sp_list:
                if sp not in cds_list.keys():
                    parent = gff_db_dict[sp][parent_groups[ortho][sp]]
                    strand = parent.strand
                    cds_list[sp] = [cds for cds in
                                    gff_db_dict[sp].children(parent, featuretype="CDS", order_by="start")]
                
                if parent.strand == "-":
                    cds_list[sp] = cds_list[sp][::-1]

                index = sp_core_exon_spans[ortho][sp].index(coord)
                cds = cds_list[sp][index]
                fasta_prep[ortho][coord][sp] = cds

    shutil.rmtree(sup_align_dir, ignore_errors=True)
    os.makedirs(sup_align_dir, exist_ok=True)
    driver_input = {}
    nnn = Seq('n' * n_gap)
    for ortho in fasta_prep:
        for coord in fasta_prep[ortho]:
            filename = os.path.join(sup_align_dir, ortho + "_" + str(coord[0]) + "-" + str(coord[1]) + ".fasta")
            driver_input[filename] = ""
            for sp in core_sp_list:
                cds = fasta_prep[ortho][coord][sp]
                start, end = coord
                seq = nnn + aligned_fasta[ortho][sp].seq[start:end] + nnn
                seqRec = SeqRecord(seq, id=sp, description=cds.id)
                driver_input[filename] += seqRec.format("fasta")

            for sp in sorted(parent_groups[ortho]):

                if sp in core_sp_list:
                    continue

                if (not add_outgroup_in_sup_align) and (sp in outgroup_sp_list):
                    continue

                parent = gff_db_dict[sp][parent_groups[ortho][sp]]
                strand = parent.strand
                cds_list = [cds for cds in gff_db_dict[sp].children(parent, featuretype="CDS", order_by="start")]
                cat_seq = Seq("")

                for cds in cds_list:
                    try:
                        cat_seq += Seq(str(cds.sequence(fasta=fasta_dict[sp], use_strand=False)))
                    except ValueError as e:
                        if "imply a diffent length than sequence" in str(e):
                            cat_seq += Seq(str(fasta_dict[sp][cds.chrom]))
                            logging.debug(
                                "coordinates from gff don't fall within with scaffold, grabbing entire scaffold")
                        else:
                            raise

                if strand == '-':
                    cat_seq = cat_seq.reverse_complement()
                seqRec = SeqRecord(cat_seq, id=sp, description=parent.id)
                driver_input[filename] += seqRec.format("fasta")

    def write_2nd_mafft_input(filename, content):
        with open(filename, "w") as f:
            f.write(content)

    pool = ThreadPool(mp.cpu_count())
    pool.starmap(write_2nd_mafft_input, driver_input.items())
    pool.close()
    pool.join()


#####################
# supplemental post #
#####################


def maft_supp_post(enhanced_alignment_path, final_exon_dir, core_species_list, min_exon_length, max_gap_percent, max_gap_length, min_num_sp, n_gap):
    # create handles for all .fasta files in aligned_full_fasta directory
    aligned_fasta_fn = {name.split('.fasta')[0]: os.path.join(enhanced_alignment_path, name) for name in
                        os.listdir(enhanced_alignment_path) if
                        ((".fasta.aln" in name) and (".fasta.aln.fai" not in name))}

    # read and parse fasta files for each species
    aligned_fasta = {}
    for ortho in aligned_fasta_fn.keys():
        aligned_fasta[ortho] = {seq_record.id: seq_record
                                for seq_record in SeqIO.parse(aligned_fasta_fn[ortho],
                                                              "fasta")}

    # parse coords from core species in aligned fasta's and trash entries w/ all gaps
    coords = {}  # coords[ortho][sp] = [coord, ]
    for ortho in aligned_fasta:
        coords[ortho] = {}
        for sp in core_species_list:
            seq = str(aligned_fasta[ortho][sp].seq)
            core_coords = findExonCoords(str(aligned_fasta[ortho][sp].seq), n_gap)
            for start, end in core_coords:
                cds = seq[start:end]
                if len(cds) != cds.count('-'):
                    if sp not in coords[ortho]:
                        coords[ortho][sp] = (start, end)
                    elif type(coords[ortho][sp]) is list:
                        coords[ortho][sp].append((start, end))
                    else:
                        core = coords[ortho][sp]
                        coords[ortho][sp] = [core, (start, end)]

    # sanity check for multiple non gap core cds's per ortho,sp
    for ortho in coords:
        for sp in coords[ortho]:
            if type(coords[ortho][sp]) is list:
                 del coords[ortho][sp]
                 logging.error("problem parsing coords. multiple non-gap core cds's for {},{}: {}".format(ortho, sp, coords[ortho][sp]))

    # Filter aligned exons
    ortho_coords = {}
    for ortho in coords:
        ortho_coords[ortho] = {}
        for sp in coords[ortho]:
            coord = coords[ortho][sp]

            # filter for length
            start, end = coord
            length = end - start
            # if not min_exon_length <= length
            if not min_exon_length <= length:
                continue

            # filter for gap percent
            seq = str(aligned_fasta[ortho][sp].seq[start:end])
            if gapPercent(seq) > max_gap_percent:
                continue

            # filter for gap length
            if longestGap(seq) > max_gap_length:
                continue

            # prep to filter for species membership of ortho
            if coord not in ortho_coords[ortho].keys():
                ortho_coords[ortho][coord] = set()
            ortho_coords[ortho][coord].add(sp)

    # set of coords per ortho which were represented in all species
    universal_ortho_coords = {}
    for ortho in ortho_coords:
        for coord in ortho_coords[ortho]:
            sp_set = ortho_coords[ortho][coord]
            if len(sp_set) == len(core_species_list):
                if ortho not in universal_ortho_coords.keys():
                    universal_ortho_coords[ortho] = set()
                universal_ortho_coords[ortho].add(coord)
#            else:
#                logging.info("2nd alignment broke up {} {}, has only {}. excluding".format(ortho, coord, sp_set))

    # fasta prep
    fasta_prep = {}
    for ortho in universal_ortho_coords:
        fasta_prep[ortho] = []
        for coord in universal_ortho_coords[ortho]:
            core_sp_list = []
            for sp in sorted(aligned_fasta[ortho]):
                start, end = coord
                seq = aligned_fasta[ortho][sp].seq[start:end]
                des = aligned_fasta[ortho][sp].description
                seqReq = SeqRecord(seq, id=sp, description=des)
                if sp in core_species_list:
                    fasta_prep[ortho].append(seqReq)
                else:
                    core_sp_list.append(seqReq)

            fasta_prep[ortho].extend(core_sp_list)

    for ortho in fasta_prep:
        fasta_prep[ortho] = [seqReq for seqReq in fasta_prep[ortho] if
                             (gapPercent(seqReq.seq) <= max_gap_percent) and (
                                 longestGap(seqReq.seq) <= max_gap_length)]

    fasta_prep = {ortho: seq_list for ortho, seq_list in fasta_prep.items() if len(seq_list) >= min_num_sp}

    # fasta output
    shutil.rmtree(final_exon_dir, ignore_errors=True)
    os.makedirs(final_exon_dir, exist_ok=True)
    for ortho in fasta_prep:
        filename = os.path.join(final_exon_dir, ortho + ".fasta")
        with open(filename, "w") as f:
            for seq_req in fasta_prep[ortho]:
                f.write(seq_req.format("fasta"))


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--configPath', help='configPath', default='./config.ini')
    args = parser.parse_args()

    # load parameters from config file
    config = configparser.ConfigParser()
    config.read(args.configPath)

    # input paths
    gff_db_dir = config["Paths"]["gff_db_dir"]
    fasta_dir = config["Paths"]["fasta_dir"]
    tsv_dir = config["Paths"]["tsv_dir"]

    # output paths
    working_dir = config["Paths"]["working_dir"]
    os.makedirs(working_dir, exist_ok=True)

    core_align_dir = config["Paths"]["core_align_dir"]
    os.makedirs(core_align_dir, exist_ok=True)

    sup_align_dir = config["Paths"]["sup_align_dir"]
    os.makedirs(sup_align_dir, exist_ok=True)

    final_exon_dir = config["Paths"]["final_exon_dir"]
    os.makedirs(final_exon_dir, exist_ok=True)

    # remove old files
    out_dirs = [core_align_dir, sup_align_dir, final_exon_dir]
    for dir in out_dirs:
        rm_files = [os.path.join(dir, file) for file in os.listdir(dir)]
        for file in rm_files:
            os.remove(file)

    # settings
    min_exon_length = config["Settings"].getint("min_exon_length")
    max_gap_percent = config["Settings"].getint("max_gap_percent")
    max_gap_length = config["Settings"].getint("max_gap_length")
    min_num_sp = config["Settings"].getint("min_num_sp")
    add_outgroup_in_core_align = config["Settings"].getboolean("add_outgroup_in_core_align")
    add_outgroup_in_sup_align = config["Settings"].getboolean("add_outgroup_in_sup_align")

    # species lists
    core_sp_list = config["Species Groups"]["core_sp_list"].split()
    outgroup_sp_list = config["Species Groups"]["outgroup_sp_list"].split()


    gff_db_dict = dict()
    for db_fn in os.listdir(gff_db_dir):
        sp = os.path.splitext(db_fn)[0]
        db = gffutils.FeatureDB(os.path.join(gff_db_dir, db_fn), keep_order=True)
        gff_db_dict[sp] = db

    fasta_dict = dict()
    for fasta_fn in os.listdir(fasta_dir):
        if ".fai" in fasta_fn or ".fasta" not in fasta_fn:
            continue
        sp = os.path.splitext(fasta_fn)[0]
        fasta = Fasta(os.path.join(fasta_dir, fasta_fn))
        fasta_dict[sp] = fasta


    # run everything
    maft_core_prep(gff_db_dict, fasta_dict, working_dir, tsv_dir, core_align_dir, add_outgroup_in_core_align,
                   set(core_sp_list), set(outgroup_sp_list), n_gap)

    himap.external_software.core_mafft_driver_path(core_align_dir)

    maft_core_post_and_supp_prep(gff_db_dict, fasta_dict, working_dir, core_align_dir, sup_align_dir,
                                 add_outgroup_in_sup_align, core_sp_list, outgroup_sp_list, n_gap, min_exon_length, max_gap_percent, max_gap_length)

    himap.external_software.supplemental_mafft_driver_path(sup_align_dir)

    maft_supp_post(sup_align_dir, final_exon_dir, core_sp_list, min_exon_length, max_gap_percent, max_gap_length, min_num_sp, n_gap)
    
    print("Step_04 is complete\nPadded exons are written to: ", core_align_dir, "\n", len(os.listdir(sup_align_dir)), "Raw exons are written to: ", sup_align_dir, "\n", len(os.listdir(final_exon_dir)), "Final filtered exons are written to: ", final_exon_dir)

# run main
try:
    exit(main())
except Exception as e:
    log.exception("Exception in main(): {}".format(e))
    exit(1)
