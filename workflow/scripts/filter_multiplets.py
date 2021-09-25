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

def main(fragments, fragments_out, multiplet_barcodes, summary, min_common=2, min_counts=500):
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
    primary_barcodes = {}
    bc_sets = {}
    for x, y in pair_counts.items():
        if y >= min_common:
            a, b = x
            bca = barcode_counts[a]
            bcb = barcode_counts[b]
            multiplet_data[x] = (a, b, bca, bcb, y, y/(bca + bcb - y), None)

            a_primary = primary_barcodes.setdefault(a, a)
            b_primary = primary_barcodes.setdefault(b, b)
            a_set = bc_sets.setdefault(a_primary, set(a)) # initialize starting sets if needed
            b_set = bc_sets.setdefault(b_primary, set(b))
            set_info = {a_primary: a_set, b_primary: b_set}
            

            if a_primary != b_primary:
                remaining_primary = max([a_primary, b_primary], key=barcode_counts.get) 
                other_primary = a_primary if remaining_primary == b_primary else b_primary   
                for k in set_info[other_primary]:
                    primary_barcodes[k] = remaining_primary

                a_set |= b_set
                bc_sets[remaining_primary] = a_set
                bc_sets.pop(other_primary) 

    blacklist = set()
    with open(multiplet_barcodes, 'w') as f:
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
        "Identified {:,} barcodes belonging to multiplets".format(
            len(barcode_counts)
        ),
        logout,
        starttime,
    )

    print_and_log(
        "Identified {:,} unique multiplets".format(
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
    multiplet_barcodes = snakemake.output['barcodes']
    summary = snakemake.output['qc']

    main(fragments, fragments_index, fragments_out, multiplet_barcodes, summary)

except NameError:
    pass 