import argparse
import gzip
from pathlib import Path

import pandas as pd
import numpy as np

import matcha

# Processes 10x ATAC-seq fastqs
# Inputs: 10x fastq files R1,R2,R3, and barcode whitelist
# Outputs: 
#   - R1.fastq.gz: Input R1 filtered for matching barcodes, with corrected barcode
#     appended to the read name
#   - R2.fastq.gz: Input R3 (NOT R2) filtered for matching barcodes, with corrected barcode
#     appended to the read name
#
# Limitations:
#  - Current QC output is just printing a few lines to stdout at the end of the run
#  - Currently doesn't have arguments for the output filenames
#  - Only outputs corrected barcode, skipping the uncorrected barcode.

def main():
    parser = argparse.ArgumentParser(description='Process raw fastq reads to append cell barcodes')
    parser.add_argument('fastq_I1', type=str, help='I1 fastq file')
    parser.add_argument('fastq_I2', type=str, help='I2 fastq file')
    parser.add_argument('fastq_R1', type=str, help='R1 fastq file')
    parser.add_argument('fastq_R2', type=str, help='R2 fastq file')
    parser.add_argument('i5_barcodes', type=str, help='tsv file with i5 barcodes')
    parser.add_argument('T7_barcodes', type=str, help='tsv file with T7 barcodes')
    parser.add_argument('output_dir', type=str, help="Directory for saving outputs (R1/R2/stats)")

    parser.add_argument("--reverse-complement-I1", action="store_true", help="Whether to reverse complement the I1 sequence")
    parser.add_argument("--reverse-complement-I2", action="store_true", help="Whether to reverse complement the I2 sequence")
    parser.add_argument("--gzip", action="store_true", help="If set, output fastq.gz files rather than fastq")
    parser.add_argument("--max-barcode-dist", type=int, default=2, help="Maximum edit distance allowed between a barcode and a whitelist entry")
    parser.add_argument("--ncores", default=4, help="Number of cores for parallel processing")
    
    args = parser.parse_args()

    output_path = Path(args.output_dir)

    f = matcha.FastqReader(threads = args.ncores)
    extension = "fastq.gz" if args.gzip else "fastq"
    f.add_sequence("R1", args.fastq_R1, output_path=output_path / f"R1.{extension}")
    f.add_sequence("R2", args.fastq_R2, output_path=output_path / f"R2.{extension}")
    f.add_sequence("I1", args.fastq_I1)
    f.add_sequence("I2", args.fastq_I2)

    valid_i5_barcodes = pd.read_csv(args.i5_barcodes, sep="\t")
    valid_T7_barcodes = pd.read_csv(args.T7_barcodes, sep="\t")

    if args.reverse_complement_I2:
        i5_sequences = [reverse_complement(b.encode()) for b in valid_i5_barcodes["sequence"]]
    else:
        i5_sequences = valid_i5_barcodes["sequence"]
    
    i5_barcode = matcha.HashMatcher(
        sequences = i5_sequences,
        labels = valid_i5_barcodes["sequence"],
        max_mismatches=args.max_barcode_dist,
        subsequence_count=2
    )

    if args.reverse_complement_I1:
        T7_sequences = [reverse_complement(b.encode()) for b in valid_T7_barcodes["sequence"]]
    else:
        T7_sequences = valid_T7_barcodes["sequence"]

    T7_barcode = matcha.HashMatcher(
        sequences = T7_sequences,
        labels = valid_T7_barcodes["sequence"],
        max_mismatches=args.max_barcode_dist,
        subsequence_count=2
    )

    f.add_barcode("i5", i5_barcode, "I2")
    f.add_barcode("T7", T7_barcode, "I1")

    f.set_output_names("{read_name} CB:Z:{i5}{T7}")

    barcode_counts = np.zeros((args.max_barcode_dist + 2, args.max_barcode_dist + 2), int)
    
    total_reads = 0
    total_pass = 0

    chunk_size = 10000

    dists = np.zeros((2, chunk_size), int)
    second_dists = np.zeros((2, chunk_size), int)
    while f.read_chunk(chunk_size):
        dists[0,:] = f.get_match_result("i5", "dist")
        second_dists[0,:] = f.get_match_result("i5", "second_best_dist")
        dists[1,:] = f.get_match_result("T7", "dist")
        second_dists[1,:] = f.get_match_result("T7", "second_best_dist")

        pass_filter = (dists[0] < args.max_barcode_dist) & \
            (dists[1] < args.max_barcode_dist) & \
            (dists[0] + dists[1] < second_dists[0] + second_dists[1])
            
        #import pdb; pdb.set_trace()

        total_reads += len(pass_filter)
        total_pass += pass_filter.sum()

        values, counts = np.unique(dists, axis = 1, return_counts=True)
        indices = np.minimum(values, args.max_barcode_dist+1)
        barcode_counts[(indices[0], indices[1])] += counts
        
        f.write_chunk(pass_filter)
    
    stats_output = open(output_path / "matching_stats.tsv", "w")
    print(f"{total_pass}/{total_reads} reads passing, ({total_pass/total_reads*100:.2f}%)\n", file=stats_output)
    print("mismatches_i5\tmismatches_T7\treads", file=stats_output)
    for i5_dist in range(args.max_barcode_dist + 2):
        for T7_dist in range(args.max_barcode_dist + 2):
            print(
                i5_dist if i5_dist <= args.max_barcode_dist else f">{args.max_barcode_dist}",
                T7_dist if T7_dist <= args.max_barcode_dist else f">{args.max_barcode_dist}",
                barcode_counts[i5_dist, T7_dist],
                sep = "\t",
                file=stats_output
            )

rev_comp = bytes.maketrans(b"ATGC",b"TACG")
def reverse_complement(seq):
    return bytes.translate(seq, rev_comp)[::-1]


if __name__ == "__main__":
    main()