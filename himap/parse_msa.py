#!/usr/bin/env python
# coding=utf-8
""""""

import re

from context import himap, project_dir, n_gap
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
        return 0


def find_exon_spans(seq, intron_spans):
    length = len(seq)

    if len(intron_spans) == 0:
        return [(0, length)]

    elif len(intron_spans) == 1:
        bp = intron_spans[0]
        return [(0, bp[0]), (bp[1], length)]

    elif len(intron_spans) > 1:
        exonCoords = [(0, intron_spans[0][0])]  # first exon

        for i in list(range(len(intron_spans) + 1))[1:-1]:  # all intermediate exons
            ex_start = intron_spans[i - 1][1]
            ex_end = intron_spans[i][0]
            exonCoords.append((ex_start, ex_end))

        exonCoords.append((intron_spans[-1][1], length))  # last exon
        return exonCoords


def universal_exon_spans(seq_dict, n_gap, min_exon_length, max_gap_percent, max_gap_length):
    re_no_end_gaps = re.compile(r"n[-n]*n[-n]n*")
    re_end_gaps = re.compile(r"[-n]*n[-n]*")

    # parse the set of canidate exon endpoints which are conserved
    sp_canidate_endpoints = {}
    for sp, seq in seq_dict.items():
        sp_canidate_endpoints[sp] = set()
        loc = 0
        while True:
            match_no_end_gaps = re_no_end_gaps.search(str(seq.seq), loc)
            match_end_gaps = re_end_gaps.search(str(seq.seq), loc)
            if not (match_no_end_gaps and match_end_gaps):
                break
            if (len(match_no_end_gaps.group().replace('-', '')) >= n_gap and
                   len(match_end_gaps.group().replace('-', '')) >= n_gap):
                sp_canidate_endpoints[sp].add(match_no_end_gaps.start())
                sp_canidate_endpoints[sp].add(match_no_end_gaps.end())
                sp_canidate_endpoints[sp].add(match_end_gaps.start())
                sp_canidate_endpoints[sp].add(match_end_gaps.end())
            loc = match_no_end_gaps.end()
    universal_canidate_endpoints = set.intersection(*sp_canidate_endpoints.values())
    # universal_canidate_endpoints

    sp_intron_spans = {}
    for sp, seq in seq_dict.items():
        sp_intron_spans[sp] = []
        loc = 0
        while True:
            match_no_end_gaps = re_no_end_gaps.search(str(seq.seq), loc)
            match_end_gaps = re_end_gaps.search(str(seq.seq), loc)
            if not (match_no_end_gaps and match_end_gaps):
                break
            if (len(match_no_end_gaps.group().replace('-', '')) >= n_gap and
                   len(match_end_gaps.group().replace('-', '')) >= n_gap):

                start_set = {match_no_end_gaps.start(), match_end_gaps.start()} & universal_canidate_endpoints
                end_set = {match_no_end_gaps.end(), match_end_gaps.end()} & universal_canidate_endpoints

                if len(start_set) == 0:
                    start = match_no_end_gaps.start()  # no start in universal_canidate_endpoints, be conservative
                else:
                    start = sorted(start_set)[0]  # first start in universal_canidate_endpoints

                if len(end_set) == 0:
                    end = match_no_end_gaps.end()  # no end in universal_canidate_endpoints, be conservative
                else:
                    end = sorted(end_set)[-1]  # last end in universal_canidate_endpoints
                                    
                sp_intron_spans[sp].append((start, end))

            loc = match_no_end_gaps.end()

    #sp_intron_spans

    sp_exon_spans = {}
    for sp,seq in seq_dict.items():
        sp_exon_spans[sp] = find_exon_spans(seq, sp_intron_spans[sp])


    universal_exon_spans = sorted(set.intersection(*[set(lis) for lis in sp_exon_spans.values()]), key=lambda x: x[0])
    #universal_exon_spans

    # filter
    filtered_universal_exon_spans = []
    for start,end in universal_exon_spans:
        seq_list = [str(seq_dict[sp].seq[start:end]) for sp in seq_dict]
        consensus = himap.consensus.degenerate_consensus(seq_list).strip("-")

        # filter for length
        if len(consensus) < min_exon_length:
            continue

        # filter for gap percent
        if gap_percent(consensus) > max_gap_percent and longest_gap(consensus) > max_gap_length:
            continue

        filtered_universal_exon_spans.append((start,end))
    return sp_exon_spans, filtered_universal_exon_spans
