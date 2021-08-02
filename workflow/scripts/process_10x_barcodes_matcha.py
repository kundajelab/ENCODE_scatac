import argparse
import gzip

import numpy as np
import tqdm

import matcha

def main():
    parser = argparse.ArgumentParser(description='Process raw fastq reads to append cell barcodes')
    parser.add_argument('fastq1', type=str, help='R1 fastq file')
    parser.add_argument('fastq2', type=str, help='R2 fastq file')
    parser.add_argument('fastq3', type=str, help='R3 fastq file')
    parser.add_argument('barcodes', type=str, help='Gzip file with 10x barcode whitelist')

    parser.add_argument("--max-barcode-dist", type=int, default=2, help="Maximum edit distance allowed between a barcode and a whitelist entry")
    parser.add_argument("--ncores", default=4, help="Number of cores for parallel processing")
    
    args = parser.parse_args()
    f = matcha.FastqReader(threads = args.ncores)
    f.add_sequence("R1", args.fastq1, output_path="R1.fastq.gz")
    f.add_sequence("R2", args.fastq2)
    f.add_sequence("R3", args.fastq3, output_path="R2.fastq.gz")

    valid_10x_barcodes = [b.strip() for b in gzip.open(args.barcodes)]

    cell_barcode = matcha.HashMatcher(
        valid_10x_barcodes,
        max_mismatches=args.max_barcode_dist,
        subsequence_count=2
    )

    f.add_barcode("cell", cell_barcode, "R2")


    f.set_output_names("{read_name} CB:Z:{cell}")

    barcode_counts = [0] * (args.max_barcode_dist + 1)
    total_reads = 0

    chunk_size = 10000
    progress = tqdm.tqdm(disable=None, unit="reads")
    while f.read_chunk(chunk_size):
        pass_filter = (f.get_match_result("cell", "dist") <= args.max_barcode_dist) & \
            (f.get_match_result("cell", "second_best_dist") > f.get_match_result("cell", "dist"))

        total_reads += len(pass_filter)
        values, counts = np.unique(f.get_match_result("cell", "dist"), return_counts=True)
        for v, c in zip(values, counts):
            if v <= args.max_barcode_dist:
                barcode_counts[v] += c
        
        f.write_chunk(pass_filter)
        progress.update(len(pass_filter))

    for val, i in enumerate(barcode_counts):
        print(val, i)

if __name__ == "__main__":
    main()