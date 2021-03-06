import gzip
from pathlib import Path

import numpy as np

import matcha

# Processes 10x ATAC-seq fastqs
# Inputs: 10x fastq files R1,R2,R3, and barcode whitelist
# Outputs: 
#   - R1.fastq: Input R1 filtered for matching barcodes, with corrected barcode
#     appended to the read name
#   - R2.fastq: Input R3 (NOT R2) filtered for matching barcodes, with corrected barcode
#     appended to the read name
#
# Limitations:
#  - Current QC output is just printing a few lines to stdout at the end of the run
#  - Currently doesn't have arguments for the output filenames
#  - Only outputs corrected barcode, skipping the uncorrected barcode.

def main():
    parser = argparse.ArgumentParser(description='Process raw fastq reads to append cell barcodes')
    parser.add_argument('fastq1', type=str, help='R1 fastq file')
    parser.add_argument('fastq2', type=str, help='R2 fastq file')
    parser.add_argument('fastq3', type=str, help='R3 fastq file')
    parser.add_argument('barcodes', type=str, help='Text file with 10x barcode whitelist (optionally gzipped)')
    parser.add_argument('output_dir', type=str, help="Directory for saving outputs (R1/R2/stats)")

    parser.add_argument("--reverse-complement", action="store_true", help="Whether to reverse complement the R2 sequence (NextSeq support)")
    parser.add_argument("--gzip", action="store_true", help="If set, output fastq.gz files rather than fastq")
    parser.add_argument("--max-barcode-dist", type=int, default=2, help="Maximum edit distance allowed between a barcode and a whitelist entry")
    parser.add_argument("--ncores", default=4, help="Number of cores for parallel processing")
    
    args = parser.parse_args()

    output_path = Path(args.output_dir)

    f = matcha.FastqReader(threads = args.ncores)
    extension = "fastq.gz" if args.gzip else "fastq"
    f.add_sequence("R1", args.fastq1, output_path=output_path / f"R1.{extension}")
    f.add_sequence("R2", args.fastq2)
    f.add_sequence("R3", args.fastq3, output_path=output_path / f"R2.{extension}")

    valid_10x_barcodes = [b.strip() for b in read_gzip(args.barcodes)]
    if args.reverse_complement:
        barcode_sequences = [reverse_complement(b) for b in valid_10x_barcodes]
    else:
        barcode_sequences = valid_10x_barcodes
    cell_barcode = matcha.HashMatcher(
        sequences = barcode_sequences,
        labels = valid_10x_barcodes,
        max_mismatches=args.max_barcode_dist,
        subsequence_count=2
    )

    f.add_barcode("cell", cell_barcode, "R2")

    f.set_output_names("{read_name} CB:Z:{cell}")

    barcode_counts = np.zeros(args.max_barcode_dist + 2, int)

    total_reads = 0
    total_pass = 0

    chunk_size = 10000
    while f.read_chunk(chunk_size):
        pass_filter = (f.get_match_result("cell", "dist") <= args.max_barcode_dist) & \
            (f.get_match_result("cell", "second_best_dist") > f.get_match_result("cell", "dist"))

        total_reads += len(pass_filter)
        total_pass += pass_filter.sum()
        values, counts = np.unique(f.get_match_result("cell", "dist"), return_counts=True)
        barcode_counts[np.minimum(values, args.max_barcode_dist + 1)] += counts
        
        f.write_chunk(pass_filter)

    stats_output = open(output_path / "matching_stats.tsv", "w")
    print(f"{total_pass}/{total_reads} reads passing, ({total_pass/total_reads*100:.2f}%)\n", file=stats_output)
    print("mismatches\treads", file=stats_output)
    for dist in range(args.max_barcode_dist + 2):
        print(
            dist if dist <= args.max_barcode_dist else f">{args.max_barcode_dist}",
            barcode_counts[dist],
            sep = "\t",
            file=stats_output
        )

def read_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    else:
        return open(path, "rb")

rev_comp = bytes.maketrans(b"ATGC",b"TACG")
def reverse_complement(seq):
    return bytes.translate(seq, rev_comp)[::-1]

if __name__ == "__main__":
    main()
