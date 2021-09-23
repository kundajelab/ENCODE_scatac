"""
Modified by Austin Wang

Original description:

This tool identifies barcodes in SingleCell ATAC 10x Genomics data likely to have been
caused by barcode multiplets or gel bead multiplets.  Detection is done by identifying
fragment adjacency similarities between the contaminant and primary barcodes.

Requires a path to a completed output folder from the CellRanger ATAC data.
In particular, we require the deduplicated and position-corrected fragments file along
with its index.

Required memory: ~16GB.

Version 1.0, released November 4, 2019

Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""

# from __future__ import absolute_import, division, print_function

# import argparse
import pysam
import os
import sys
import time
import gzip
import pandas as pd
import numpy as np
from collections import Counter


MAXIMUM_STOP_DICT_SIZE = 2500
MAXIMUM_FRAGMENT_SIZE = 800
MAXIMUM_PILEUP = 20

MINIMUM_COUNTS_FOR_GB_MULTIPLETS = 250
MAXIMUM_BARCODES_FOR_GB_MULTIPLETS = 15000

MINIMUM_COUNTS_FOR_BC_MULTIPLETS = 100

SELF_SIGNAL_THRESHOLD_MULTIPLIER = 0.4


# def query_path_for_file(path, filename):
#     if os.path.exists(os.path.join(path, filename)):
#         return os.path.join(path, filename)
#     elif os.path.exists(os.path.join(path, "outs", filename)):
#         return os.path.join(path, "outs", filename)
#     else:
#         print("Could not find required file {}.".format(filename))
#         print(
#             "Check that the input path is correct and that this tool is being "
#             "run on a completed CellRanger ATAC analysis."
#         )
#         sys.exit(1)


# def check_output_path(path, prefix, filename):
#     full_path = os.path.abspath(os.path.join(path, "{}{}".format(prefix, filename)))
#     if not os.access(path, os.W_OK):
#         print(
#             "Putative output file {} is not writeable.  Please check file permissions.".format(
#                 full_path
#             )
#         )
#         sys.exit(1)
#     return full_path


# def parse_args():
#     parser = argparse.ArgumentParser(description=__doc__)
#     parser.add_argument(
#         "input_path",
#         type=str,
#         help="Path to the output directory of the CellRanger ATAC analysis to be corrected.",
#     )
#     parser.add_argument(
#         "--output_path",
#         type=str,
#         default="./",
#         required=False,
#         help="Path where results will be written.",
#     )
#     parser.add_argument(
#         "--prefix",
#         type=str,
#         default="",
#         required=False,
#         help="Prefix to be added to output files.",
#     )
#     args = parser.parse_args()

#     # Check that the path has the correct outputs needed
#     args.fragments = query_path_for_file(args.input_path, "fragments.tsv.gz")
#     args.fragments_index = query_path_for_file(args.input_path, "fragments.tsv.gz.tbi")
#     args.singlecell = query_path_for_file(args.input_path, "singlecell.csv")

#     # Hook up the outputs and make sure we can write them out
#     args.cell_barcodes = check_output_path(
#         args.output_path, args.prefix, "cell_barcodes.csv"
#     )
#     args.summary = check_output_path(args.output_path, args.prefix, "summary.txt")
#     args.excluded_barcodes = check_output_path(
#         args.output_path, args.prefix, "excluded_barcodes.csv"
#     )

#     return args


# def get_barcode_gem_group(barcode):
#     """10x Genomics barcodes have a key after the barcode sequence that identifies
#     the original library the barcode was derived from.  This returns that key.
#     Only necessary for aggregated data (single samples should all have the same key).
#     """
#     split = barcode.split("-")
#     assert len(split) == 2
#     barcode, gem_group = split
#     return gem_group


def parsed_fragments_from_file(fragments_filename, fragments_index=None):
    """Iterates through all fragments in an index input file."""
    fragments = pysam.TabixFile(fragments_filename, index=fragments_index)
    for fragment in fragments.fetch():
        contig, start, stop, barcode, dup_count = fragment.split("\t")
        yield contig, int(start), int(stop), barcode


def group_fragments_by_start(filename, index=None):
    """Iterates through fragments, grouping up sets of fragments that share the same
    start position.
    """
    current = None
    fragment_list = []
    for fragment in parsed_fragments_from_file(filename, index):
        contig, start = fragment[0], fragment[1]
        if current is None:
            current = (contig, start)
        if current == (contig, start):
            fragment_list.append(fragment)
        else:
            yield current[0], fragment_list
            current = (contig, start)
            fragment_list = [fragment]
    if fragment_list:
        yield current[0], fragment_list


def iterate_through_fragment_adjacencies(
    filename, index=None, logoutfile=None, starttime=0
):
    """Iterates through the full fragments file, identifying all pairs of fragment
    linkages between barcodes and yielding them iteratively.
    """
    stops = {}
    current_contig = None
    for contig, fragment_list in group_fragments_by_start(filename, index):

        if current_contig != contig:
            # We've rolled over to a new contig, reset the stops dictionary
            current_contig = contig
            stops = {}
            if logoutfile is not None:
                print_and_log(
                    "Beginning processing of contig {}".format(contig),
                    logoutfile,
                    starttime,
                )

        for _, start, stop, barcode in fragment_list:

            # Periodically clean up old keys to reduce the memory footprint
            if len(stops) > MAXIMUM_STOP_DICT_SIZE:
                for key in stops.keys():
                    if key < (start - MAXIMUM_FRAGMENT_SIZE):
                        stops.pop(key)

            # Add each fragment into the stops dictionary
            if stop not in stops:
                stops[stop] = Counter()
            stops[stop][barcode] += 1

            # Check each fragment for any adjacency links
            if start not in stops:
                continue

            # Avoid linking fragments at very high pileup regions due to likely
            # erroneous connections
            if len(stops[stop]) + len(fragment_list) > MAXIMUM_PILEUP:
                continue

            # gem_group = get_barcode_gem_group(barcode)
            for linked_barcode, count in stops[start].iteritems():
                # # Avoid linking barcodes from different libraries
                # if gem_group != get_barcode_gem_group(linked_barcode):
                #     continue

                # Yield the linkage, along with strength of it
                yield barcode, linked_barcode, count


def split_barcode(barcode, return_gg=False):
    """Splits a barcode sequence into part A, part C, and part B sequences."""
    split = barcode.split("-")
    if len(split) == 2:
        # Remove the gem group from the barcode
        barcode, gem_group = split
    else:
        gem_group = None
    assert len(barcode) == 16
    part_a, part_c, part_b = barcode[:7], barcode[7:9], barcode[9 : 9 + 7]
    if return_gg:
        return part_a, part_c, part_b, gem_group
    return part_a, part_c, part_b


def merge_barcode(a, c, b, gem_group=None):
    barcode = "{}{}{}".format(a, c, b)
    if gem_group is not None:
        return "{}-{}".format(barcode, gem_group)
    return barcode


def query_barcode_subsequences(valid_barcodes):
    """Breaks down a list of valid barcode sequences into unique part A, B, and C
    subsequences."""
    part_a_seqs = {}
    part_b_seqs = {}
    part_c_seqs = set()
    gem_group_seqs = set()

    for barcode in valid_barcodes:
        a, c, b, gem_group = split_barcode(barcode, return_gg=True)
        part_c_seqs.add(c)
        gem_group_seqs.add(gem_group)
        if c not in part_a_seqs:
            part_a_seqs[c] = set()
        if c not in part_b_seqs:
            part_b_seqs[c] = set()
        part_a_seqs[c].add(a)
        part_b_seqs[c].add(b)

    part_c_seqs = sorted(part_c_seqs)
    for c in part_c_seqs:
        part_a_seqs[c] = sorted(part_a_seqs[c])
        part_b_seqs[c] = sorted(part_b_seqs[c])

    return part_a_seqs, part_b_seqs, part_c_seqs, gem_group_seqs


def nearest_neighbor(linkage_matrix, index):
    """Finds the non-self index in a linkage matrix with the highest linkage."""
    row = linkage_matrix[index]
    indices = np.arange(len(linkage_matrix))
    mask = indices != index
    return indices[mask][np.argmax(row[mask])]


def query_linkage_for_gb_multiplets(
    linkage_matrix, valid_barcodes, fragment_counts_by_barcode
):
    putative_multiplets = []
    for index in range(len(valid_barcodes)):
        neighbor = nearest_neighbor(linkage_matrix, index)
        if index >= neighbor:
            # Prevents reporting the same pair twice
            continue

        if nearest_neighbor(linkage_matrix, neighbor) == index:
            putative_multiplets.append(
                (valid_barcodes[index], valid_barcodes[neighbor])
            )

    exclusions = {}
    for pair in putative_multiplets:
        if fragment_counts_by_barcode[pair[0]] < fragment_counts_by_barcode[pair[1]]:
            excluded_bc, major_bc = pair
        else:
            major_bc, excluded_bc = pair
        exclusions[excluded_bc] = (major_bc, "gel_bead_multiplet")
    return exclusions


def print_and_log(text, outfile, starttime=0):
    logtime = time.process_time() - starttime
    if logtime < 60:
        logtime = "{:,}s".format(logtime)
    else:
        logtime = "{:,}m {:,}s".format(logtime // 60, logtime % 60)
    outfile.write("{} - {}\n".format(logtime, text))
    print("{} - {}".format(logtime, text))


def get_barcodes(fragments):
    bc_counter = Counter()
    with gzip.open(fragments, "rt") as f:
        for line in f:
            entries = line.strip("\n").split("\t")
            bc = entries[3]
            count = int(entries[4])

            bc_counter[bc] += count
    return pd.DataFrame.from_records([(k, v) for k, v in bc_counter.items()], columns=['barcode', 'passed_filters'])


def filter_fragments(frag_in, frag_out, exclusions):
    blacklist = set(exclusions)
    with gzip.open(frag_in, "rt") as fi, open(frag_out, "wt") as fo:
        for line in fi:
            entries = line.strip("\n").split("\t")
            bc = entries[3]
            if bc in blacklist:
                print(f"exc {bc}") ####
                continue
            fo.write(line)


def main(fragments, fragments_index, fragments_out, excluded_barcodes, summary):
    logout = open(summary, "w")
    starttime = time.process_time() 

    print_and_log(
        "Identifying candidate barcode for gel bead multiplet detection...",
        logout,
        starttime,
    )
    singlecell_df = get_barcodes(fragments)
    # Query the singlecell file to identify barcodes to examine for possible gel bead
    # multiplets
    mincounts = MINIMUM_COUNTS_FOR_GB_MULTIPLETS
    all_barcodes = singlecell_df["barcode"]
    valid_mask = all_barcodes != "NO_BARCODE"
    fragment_counts = singlecell_df["passed_filters"]
    fragment_counts_by_barcode = {
        bc: count
        for bc, count in zip(all_barcodes[valid_mask], fragment_counts[valid_mask])
    }

    while True:
        count_mask = fragment_counts > mincounts
        count_mask &= valid_mask
        if sum(count_mask) < MAXIMUM_BARCODES_FOR_GB_MULTIPLETS:
            break
        mincounts += 50
    valid_barcodes_for_gbmult = all_barcodes[count_mask].values
    print_and_log(
        "Identified {:,} total barcodes for gel bead multiplet detection.\n"
        "Minimum fragment count for included barcodes of {:,}.".format(
            sum(count_mask), mincounts
        ),
        logout,
        starttime,
    )

    # Initialize the linkage matrices
    # One matrix for the GB multiplet detection of only high-fragment barcodes
    # Two for the barcode multiplet detection of only barcodes that share a common
    # prefix or suffix.
    print_and_log("Initializing linkage matrices...", logout, starttime)
    gb_index_by_barcode = {bc: i for i, bc in enumerate(valid_barcodes_for_gbmult)}
    gb_linkage_matrix = np.zeros(
        (len(valid_barcodes_for_gbmult), len(valid_barcodes_for_gbmult)),
        dtype=np.uint32,
    )
    print_and_log(
        "Gel bead linkage matrix initialized, memory used: {:,} bytes.".format(
            gb_linkage_matrix.nbytes
        ),
        logout,
        starttime,
    )

    part_a_seqs, part_b_seqs, part_c_seqs, gem_group_seqs = query_barcode_subsequences(
        all_barcodes[valid_mask]
    )
    part_a_count = max([len(part_a_seqs[c]) for c in part_c_seqs])
    part_b_count = max([len(part_b_seqs[c]) for c in part_c_seqs])
    part_c_count = len(part_c_seqs)
    gem_group_count = len(gem_group_seqs)

    index_by_part_a = {
        c: {a: i for i, a in enumerate(part_a_seqs[c])} for c in part_c_seqs
    }
    index_by_part_b = {
        c: {b: i for i, b in enumerate(part_b_seqs[c])} for c in part_c_seqs
    }
    index_by_part_c = {c: i for i, c in enumerate(part_c_seqs)}
    index_by_gem_group = {gg: i for i, gg in enumerate(gem_group_seqs)}

    part_a_linkage_matrix = np.zeros(
        (gem_group_count, part_c_count, part_b_count, part_a_count, part_a_count),
        dtype=np.uint32,
    )

    print_and_log(
        "Part A linkage matrix initialized:\nshape {}\nmemory used: {:,} bytes".format(
            part_a_linkage_matrix.shape, part_a_linkage_matrix.nbytes
        ),
        logout,
        starttime,
    )

    part_b_linkage_matrix = np.zeros(
        (gem_group_count, part_c_count, part_a_count, part_b_count, part_b_count),
        dtype=np.uint32,
    )

    print_and_log(
        "Part B linkage matrix initialized:\nshape {}\nmemory used: {:,} bytes".format(
            part_b_linkage_matrix.shape, part_b_linkage_matrix.nbytes
        ),
        logout,
        starttime,
    )

    # Loop through all fragment adjacencies and use them to build up the linkage
    # matrices.

    print_and_log("Beginning processing of fragment file...", logout, starttime)

    for barcode, linked_barcode, count in iterate_through_fragment_adjacencies(
        fragments, fragments_index, logout, starttime
    ):
        # First increment the GB linkage matrix
        if barcode in gb_index_by_barcode and linked_barcode in gb_index_by_barcode:
            index1 = gb_index_by_barcode[barcode]
            index2 = gb_index_by_barcode[linked_barcode]
            gb_linkage_matrix[index1, index2] += count
            if barcode != linked_barcode:
                # Don't double count self-linkages
                gb_linkage_matrix[index2, index1] += count

        # Next increment the partwise linkage matrices
        a1, c1, b1, gem_group1 = split_barcode(barcode, return_gg=True)
        a2, c2, b2, gem_group2 = split_barcode(linked_barcode, return_gg=True)

        if c1 != c2 or gem_group1 != gem_group2:
            continue

        shared_c_index = index_by_part_c[c1]
        shared_gemgroup_index = index_by_gem_group[gem_group1]

        if a1 == a2:
            shared_a_index = index_by_part_a[c1][a1]
            b1_index = index_by_part_b[c1][b1]
            b2_index = index_by_part_b[c1][b2]
            part_b_linkage_matrix[
                shared_gemgroup_index,
                shared_c_index,
                shared_a_index,
                b1_index,
                b2_index,
            ] += count
            if b1 != b2:
                # Don't double count self-linkages
                part_b_linkage_matrix[
                    shared_gemgroup_index,
                    shared_c_index,
                    shared_a_index,
                    b2_index,
                    b1_index,
                ] += count

        if b1 == b2:
            shared_b_index = index_by_part_b[c1][b1]
            a1_index = index_by_part_a[c1][a1]
            a2_index = index_by_part_a[c1][a2]
            part_a_linkage_matrix[
                shared_gemgroup_index,
                shared_c_index,
                shared_b_index,
                a1_index,
                a2_index,
            ] += count
            if a1 != a2:
                # Don't double count self-linkages
                part_a_linkage_matrix[
                    shared_gemgroup_index,
                    shared_c_index,
                    shared_b_index,
                    a2_index,
                    a1_index,
                ] += count

    # Call out to get gel bead multiplet exclusions
    print_and_log("Identifying gel bead multiplets", logout, starttime)
    exclusions = query_linkage_for_gb_multiplets(
        gb_linkage_matrix, valid_barcodes_for_gbmult, fragment_counts_by_barcode
    )

    # Because we need all the indices we generated, find barcode multiplets here in
    # main()
    print_and_log("Identifying barcode multiplets", logout, starttime)
    putative_barcode_multiplets = {}
    for major_barcode, count in fragment_counts_by_barcode.iteritems():
        if count < MINIMUM_COUNTS_FOR_BC_MULTIPLETS:
            # Too small to examine as a potential major barcode
            continue
        part_a, part_c, part_b, gem_group = split_barcode(major_barcode, return_gg=True)
        a_index = index_by_part_a[part_c][part_a]
        b_index = index_by_part_b[part_c][part_b]
        c_index = index_by_part_c[part_c]
        gem_group_index = index_by_gem_group[gem_group]

        for other_part_a in part_a_seqs[part_c]:
            if other_part_a == part_a:
                continue
            minor_barcode = merge_barcode(other_part_a, part_c, part_b, gem_group)
            other_a_index = index_by_part_a[part_c][other_part_a]
            self_signal = part_a_linkage_matrix[
                gem_group_index, c_index, b_index, other_a_index, other_a_index
            ]
            major_signal = part_a_linkage_matrix[
                gem_group_index, c_index, b_index, other_a_index, a_index
            ]
            if major_signal > self_signal * SELF_SIGNAL_THRESHOLD_MULTIPLIER:
                if minor_barcode not in putative_barcode_multiplets:
                    putative_barcode_multiplets[minor_barcode] = major_barcode
                else:
                    old_major = putative_barcode_multiplets[minor_barcode]
                    old_a, _, _ = split_barcode(old_major)
                    old_a_index = index_by_part_a[part_c][old_a]
                    old_signal = part_a_linkage_matrix[gem_group_index, c_index, b_index, other_a_index, old_a_index]
                    if major_signal > old_signal:
                        putative_barcode_multiplets[minor_barcode] = major_barcode

        for other_part_b in part_b_seqs[part_c]:
            if other_part_b == part_b:
                continue
            minor_barcode = merge_barcode(part_a, part_c, other_part_b, gem_group)
            other_b_index = index_by_part_b[part_c][other_part_b]
            self_signal = part_b_linkage_matrix[
                gem_group_index, c_index, a_index, other_b_index, other_b_index
            ]
            major_signal = part_b_linkage_matrix[
                gem_group_index, c_index, a_index, other_b_index, b_index
            ]
            if major_signal > self_signal * SELF_SIGNAL_THRESHOLD_MULTIPLIER:
                if minor_barcode not in putative_barcode_multiplets:
                    putative_barcode_multiplets[minor_barcode] = major_barcode
                else:
                    old_major = putative_barcode_multiplets[minor_barcode]
                    _, _, old_b = split_barcode(old_major)
                    old_b_index = index_by_part_b[part_c][old_b]
                    old_signal = part_b_linkage_matrix[gem_group_index, c_index, a_index, other_b_index, old_b_index]
                    if major_signal > old_signal:
                        putative_barcode_multiplets[minor_barcode] = major_barcode

    # Post-screen barcode multiplets for mutually linked pairs, and remove the pair
    # where we've excluded the larger barcode
    for minor_barcode in putative_barcode_multiplets.keys():
        if minor_barcode not in putative_barcode_multiplets:
            # Because we've popped it off before we got here
            continue
        major_barcode = putative_barcode_multiplets[minor_barcode]
        if (
            major_barcode in putative_barcode_multiplets
            and putative_barcode_multiplets[major_barcode] == minor_barcode
        ):
            if (
                fragment_counts_by_barcode[major_barcode]
                > fragment_counts_by_barcode[minor_barcode]
            ):
                putative_barcode_multiplets.pop(major_barcode)
            else:
                putative_barcode_multiplets.pop(minor_barcode)

    # Post-screen barcode multiplets for those where the major barcode is itself
    # linked to another barcode
    for minor_barcode, major_barcode in putative_barcode_multiplets.iteritems():
        while major_barcode in putative_barcode_multiplets:
            major_barcode = putative_barcode_multiplets[major_barcode]
            putative_barcode_multiplets[minor_barcode] = major_barcode

    # Merge these exclusions with the previous gel bead multiplet exclusions
    # Note that this will overwrite older exclusions
    for excluded_barcode, major_barcode in putative_barcode_multiplets.iteritems():
        exclusions[excluded_barcode] = (major_barcode, "barcode_multiplet")

    # Write out the excluded barcodes to a CSV
    with open(excluded_barcodes, "w") as excluded_out:
        excluded_out.write(
            "{}\t{}\t{}\n".format(
                "Excluded Barcode", "Linked Barcode", "Exclusion Reason"
            )
        )
        for excluded_barcode in sorted(exclusions.keys()):
            linked_barcode, reason = exclusions[excluded_barcode]
            excluded_out.write(
                "{}\t{}\t{}\n".format(excluded_barcode, linked_barcode, reason)
            )

    orig_cell_barcodes = all_barcodes.values

    print_and_log(
        "Original run had {:,} total cell barcodes".format(
            len(orig_cell_barcodes)
        ),
        logout,
        starttime,
    )

    new_cell_mask = np.array([
        bc in orig_cell_barcodes and bc not in exclusions
        for bc in all_barcodes
    ], dtype=int)

    print_and_log(
        "After multiplet exclusions, have {:,} total cell barcodes".format(
            sum(new_cell_mask)
        ),
        logout,
        starttime,
    )


    # # Write out the input singlecell.csv file but with modified cell barcodes
    # # Hacky way to get the species list without requiring a reference
    # species_list = [
    #     col[3:-13]
    #     for col in singlecell_df.columns
    #     if col.startswith("is_") and col.endswith("_cell_barcode")
    # ]

    # for i, species in enumerate(species_list):
    #     # is_cell_data = singlecell_df["is_{}_cell_barcode".format(species)].values

    #     # print_and_log(
    #     #     "Excluding barcodes for species {}".format(species), logout, starttime
    #     # )
    #     # cell_mask = is_cell_data == 1
    #     orig_cell_barcodes = all_barcodes[cell_mask].values
    #     print_and_log(
    #         "Original run had {:,} total cell barcodes".format(
    #             len(orig_cell_barcodes)
    #         ),
    #         logout,
    #         starttime,
    #     )

    #     new_cell_mask = np.array([
    #         bc in orig_cell_barcodes and bc not in exclusions
    #         for bc in all_barcodes
    #     ], dtype=int)

    #     print_and_log(
    #         "After multiplet exclusions, have {:,} total cell barcodes".format(
    #             sum(new_cell_mask)
    #         ),
    #         logout,
    #         starttime,
    #     )

    #     singlecell_df["is_{}_cell_barcode".format(species)] = new_cell_mask

    # with open(args.cell_barcodes, "w") as cell_out:
    #     singlecell_df.to_csv(cell_out, columns=singlecell_df.columns, header=True, index=False)

    filter_fragments(fragments, fragments_out, exclusions)
    
    logout.close()

try:
    fragments = snakemake.input['frag']
    fragments_index = snakemake.input['frag_ind']
    
    fragments_out = snakemake.output['frag']
    excluded_barcodes = snakemake.output['barcodes']
    summary = snakemake.output['qc']

    main(fragments, fragments_index, fragments_out, excluded_barcodes, summary)

except NameError:
    pass 