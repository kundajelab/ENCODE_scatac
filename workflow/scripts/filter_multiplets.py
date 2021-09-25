import time
import itertools
import gzip
from collections import Counter

def print_and_log(text, outfile, starttime=0):
    logtime = time.process_time() - starttime
    if logtime < 60:
        logtime = "{:,}s".format(logtime)
    else:
        logtime = "{:,}m {:,}s".format(logtime // 60, logtime % 60)
    outfile.write("{} - {}\n".format(logtime, text))
    print("{} - {}".format(logtime, text))

def main(fragments, fragments_out, excluded_barcodes, summary, min_common=2, min_counts=500):
    logout = open(summary, "w")
    starttime = time.process_time() 

    print_and_log("Identifying candidate barcodes", logout, starttime)

    barcode_counts = Counter()
    with gzip.open(fragments, 'rt') as f:
        for line in f:
            line = line.split('\t')
            barcode = line[3]
            barcode_counts[barcode] += 1


    pair_counts = Counter()
    
    cur_clique = set()
    cur_coord = None

    i = 0

    print_and_log("Reading fragments", logout, starttime)

    with gzip.open(fragments, 'rt') as f:
        for line in f:
            line = line.split('\t')
            chr, start, _, barcode = line[:3]

            if barcode_counts[barcode] < min_counts:
                continue

            this_coord = (chr, start) # Use left insertion
            if this_coord != cur_coord:
                for x, y in itertools.combinations(cur_clique, 2):
                    x, y = x, y if x < y else y, x
                    pair_counts[(x, y)] += 1

                cur_clique = set([barcode])
                cur_coord = this_coord

            else:
                cur_clique.add(barcode)

            i += 1
            if i%1e7==0:
                print(i)

    print_and_log("Identifying barcode multiplets", logout, starttime)

    multiplet_data = {}
    bc_set_map = {}
    for x, y in pair_counts.items():
        if y >= min_common:
            a, b = x
            bca = barcode_counts[a]
            bcb = barcode_counts[b]
            multiplet_data[x] = (a, b, bca, bcb, y, y/(bca + bcb - y), None)

            a_set = bc_set_map.setdefault(a, set()).add(a) # initialize starting sets if needed
            b_set = bc_set_map.setdefault(b, set()).add(a)
            
            if a_set is not b_set: # check if they point to same object
                a_set |= b_set # merge sets if they aren't already
                bc_set_map[b] = a_set

    primary_barcodes = {}
    for s in bc_set_map.values():
        primary_bc = max(barcode_counts[b] for b in s)
        for i in s:
            primary_barcodes[s] = primary_bc

    blacklist = set()
    with open(excluded_barcodes, 'w') as f:
        f.write("Barcode1\tBarcode2\tBarcode1Counts\tBarcode2Counts\tCommon\tJaccardIndex\tPrimaryBarcode\n")
        for x, data in multiplet_data.items():
            a, b = x
            data[-1] = primary_barcodes[a]
            f.write("{}\t{}\t{}\t{}\t{}\t{:.4f}\t{}\n".format(*data))

    print_and_log(
        "Original run had {:,} total cell barcodes".format(
            len(barcode_counts)
        ),
        logout,
        starttime,
    )

    print_and_log(
        "After multiplet exclusions, have {:,} total cell barcodes".format(
            len(barcode_counts) - len(blacklist)
        ),
        logout,
        starttime,
    )

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