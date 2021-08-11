import gzip

import numpy as np
import tqdm

import matcha

fastq1 = snakemake.input['fastq1']
fastq2 = snakemake.input['fastq_bc']
fastq3 = snakemake.input['fastq2']
bc_whitelist = snakemake.input['bc_whitelist']

max_barcode_dist = snakemake.params['max_barcode_dist']

fastq1_out_path = snakemake.output['fastq1_bc']
fastq2_out_path = snakemake.output['fastq2_bc']

qc_path = snakemake.output['qc_barcode_corr']

threads = snakemake.threads

def get_open_fn(path):
    with open(path, "rb") as f:
        is_gzipped = (f.read(2) == b'\x1f\x8b')
    return gzip.open if is_gzipped else open

def main(fastq1, fastq2, fastq3, bc_whitelist, max_barcode_dist, fastq1_out_path, fastq2_out_path, qc_path, threads):
    f = matcha.FastqReader(threads = threads)
    f.add_sequence("R1", fastq1, output_path=fastq1_out_path)
    f.add_sequence("R2", fastq2)
    f.add_sequence("R3", fastq3, output_path=fastq2_out_path)

    open_fn = get_open_fn(bc_whitelist)
    with open_fn(bc_whitelist) as f:
        valid_10x_barcodes = [b.strip() for b in f]

    cell_barcode = matcha.HashMatcher(
        valid_10x_barcodes,
        max_mismatches=max_barcode_dist,
        subsequence_count=2
    )

    f.add_barcode("cell", cell_barcode, "R2")

    f.set_output_names("{read_name} CB:Z:{cell}")

    barcode_counts = [0] * (max_barcode_dist + 1)
    total_reads = 0
    total_reads_in = 0

    chunk_size = 10000
    progress = tqdm.tqdm(disable=None, unit="reads")
    while reads_in := f.read_chunk(chunk_size):
        pass_filter = (f.get_match_result("cell", "dist") <= max_barcode_dist) & \
            (f.get_match_result("cell", "second_best_dist") > f.get_match_result("cell", "dist"))

        total_reads += len(pass_filter)
        total_reads_in += reads_in
        values, counts = np.unique(f.get_match_result("cell", "dist"), return_counts=True)
        for v, c in zip(values, counts):
            if v <= max_barcode_dist:
                barcode_counts[v] += c
        
        f.write_chunk(pass_filter)
        progress.update(len(pass_filter))

    with open(qc_path, 'w') as f:
        f.write(f"Total barcodes processed\t{total_reads_in}\n")
        f.write(f"Barcodes accepted\t{total_reads}\t{total_reads/total_reads_in*100:.2f}\n")
        for val, i in enumerate(barcode_counts):
            f.write(f"Barcodes with edit distance {i}\t{val}\t{val/total_reads_in*100:.2f}\n")

main(fastq1, fastq2, fastq3, bc_whitelist, max_barcode_dist, fastq1_out_path, fastq2_out_path, qc_path, threads)
