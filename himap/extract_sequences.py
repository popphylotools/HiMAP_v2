#!/usr/bin/env python
# coding=utf-8
""""""

import re

from context import himap, project_dir
from himap import consensus


def gap_percent(seq):
    seq = str(seq)
    gappedLen = len(seq)
    gapCount = seq.count('-')
    return (100.0 * gapCount) / gappedLen


def longest_gap(seq):
    seq = str(seq)
    gap_regex = re.compile(r"-+")
    gap_list = gap_regex.findall(seq)
    if gap_list:
        return sorted([len(gap) for gap in gap_list], reverse=True)[0]
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

        # filter for gap percent
        if gap_percent(consensus) > max_gap_percent:
            continue

        # filter for gap length
        if longest_gap(consensus) > max_gap_length:
            continue

        filtered_universal_exon_spans.append((start,end))
    return sp_exon_spans, filtered_universal_exon_spans
