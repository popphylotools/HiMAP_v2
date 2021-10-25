#!/usr/bin/env python
# coding=utf-8
""""""
import logging
import os

import Bio
import gffutils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta

log = logging.getLogger("extract_sequences.py")


def group_gemoma(iso):
    key = "_".join(iso.attributes["ID"][0].split("_")[:2])
    return key


def group_ncbi(iso):
    key = iso.attributes["Parent"][0]
    return key


def group_new_rnaseq(iso):
    key = iso.attributes["Parent"][0]
    return key


def group_old_rnaseq(iso):
    key = iso.attributes["ID"][0].split("::")[0]
    return key


group_by_funcs = {
    "GeMoMa": group_gemoma,
    "NCBI": group_ncbi,
    "TrinityOldRNAseq": group_old_rnaseq,
    "TrinityNewRNAseq": group_new_rnaseq,
}


def rank_isos_by_len_and_id(iso, db):
    iso_len = db.children_bp(iso, child_featuretype='CDS', merge=False, ignore_strand=False)
    return -iso_len, iso.attributes["ID"]


def get_frame(iso, dtype):
    if iso.frame == ".":
        return 0
    else:
        return int(iso.frame)


def driver(gff_fn, fasta_fn, path, sp):
    dtype = sp.split("_")[0]
    group_by = group_by_funcs[dtype]
    db = create_db(gff_fn, sp)
    selected_isos = select_isoforms(sp, db, group_by, rank_isos_by_len_and_id)
    key_list, gff_groups, n_seqRecs, p_seqRecs = prep_output(sp, selected_isos, db, fasta_fn, dtype)

    write_gff(os.path.join(path, "01_simple_gff", sp + ".gff"), key_list, gff_groups)
    write_nuc_fasta(os.path.join(path, "01_nuc_fasta", sp + ".fasta"), key_list, n_seqRecs)
    write_pep_fasta(os.path.join(path, "01_pep_fasta", sp + ".fasta"), key_list, p_seqRecs)

def create_db(gff_fn, sp):
    fn = os.path.join("data", "01_gff_db", sp + ".db")
    if os.path.isfile(fn):
        db = gffutils.FeatureDB(fn, keep_order=True)
    else:
        db = gffutils.create_db(data=gff_fn,
                                dbfn=fn,
                                force=True,
                                merge_strategy='merge',
                                id_spec=['ID', 'Name'])
    return db


def get_cds_parents(db, level=1):
    cds_list = list(db.features_of_type(featuretype='CDS'))
    parents = set()
    for cds in cds_list:
        parents.update(db.parents(cds, level=level))
    parents = sorted(list(parents), key=lambda p: (p.seqid, p.start, p.id))
    return parents


def select_isoforms(sp, db, group_by, rank_by):
    iso_groups = dict()
    selected_isos = dict()

    # group with group_by function
    for iso in get_cds_parents(db, level=1):
        try:
            group_key = group_by(iso)
        except KeyError as e:
            message = "KeyError - sp:{} - key:{} - error:{}".format(sp, iso.attributes["ID"], e)
            log.debug(message)
            continue
        if group_key not in iso_groups:
            iso_groups[group_key] = []
        iso_groups[group_key].append(iso)

    # sort groups with rank_by function and select first from each group
    for iso_list in iso_groups.values():
        iso = sorted(iso_list, key=lambda iso: rank_by(iso, db))[0]
        selected_isos[iso.attributes["ID"][0]] = iso

    return selected_isos


def prep_output(sp, selected_isos, db, fasta_path, dtype):
    fasta = Fasta(fasta_path)

    key_list = []
    gff_groups = dict()
    n_seqRecs = dict()
    p_seqRecs = dict()

    key_sort = dict()

    for key, iso in selected_isos.items():
        cds_list = list(db.children(iso, order_by="start", featuretype="CDS"))
        try:
            n_seq = Seq("".join([c.sequence(fasta, use_strand=False) for c in cds_list]))
        except KeyError as e:
            message = "KeyError - sp:{} - key:{} - error:{}".format(sp, key, e)
            log.debug(message)
            continue
        except ValueError as e:
            message = "ValueError - sp:{} - key:{} - error:{}".format(sp, key, e)
            log.debug(message)
            continue
        if n_seq is "":
            message = "Excluding - sp:{} - key:{} - n_seq is empty".format(sp, key)
            log.debug(message)
            continue
        if iso.strand == "-":
            n_seq = n_seq.reverse_complement()
            frame = get_frame(cds_list[-1], dtype)
        else:
            frame = get_frame(cds_list[0], dtype)

        trans_n_seq = n_seq[frame:]  # the portion of n_seq to translate
        pad = (3 - (len(trans_n_seq) % 3)) % 3
        if pad != 0:
            trans_n_seq += "N" * pad
        try:
            p_seq = trans_n_seq.translate()
        except Bio.Data.CodonTable.TranslationError as e:
            message = "TranslationError - sp:{} - key:{} - error:{}".format(sp, key, e)
            log.debug(message)
            continue
        if "*" in p_seq[:-1]:
            message = "Excluding - sp:{} - key:{} - internal stop codon".format(sp, key)
            log.debug(message)
            continue
        if len(p_seq) == 0:
            message = "Excluding - sp:{} - key:{} - p_seq is empty".format(sp, key)
            log.debug(message)
            continue

        key_sort[key] = (iso.seqid, iso.start)

        key_list.append(key)
        n_seqRecs[key] = SeqRecord(n_seq, id=key, description="")
        p_seqRecs[key] = SeqRecord(p_seq, id=key, description="")
        gff_groups[key] = list(db.parents(iso, level=1)) + [iso] + cds_list

    key_list = sorted(key_list, key=lambda k: key_sort[k])

    message = "sp: {} - isoforms kept: {} - isoforms discarded: {} - percent discarded: {}"
    log.info(message.format(sp, len(key_list), len(selected_isos) - len(key_list),
                            (len(selected_isos) - len(key_list)) / len(selected_isos) * 100))

    return key_list, gff_groups, n_seqRecs, p_seqRecs


def write_gff(path, key_list, gff_groups):
    with open(path, "w") as f:
        for key in key_list:
            for gff_line in gff_groups[key]:
                f.write(str(gff_line) + "\n")


def write_nuc_fasta(path, key_list, n_seqRecs):
    with open(path, "w") as f:
        for key in key_list:
            f.write(n_seqRecs[key].format("fasta"))


def write_pep_fasta(path, key_list, p_seqRecs):
    with open(path, "w") as f:
        for key in key_list:
            f.write(p_seqRecs[key].format("fasta"))
