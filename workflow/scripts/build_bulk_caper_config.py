import json
import os
import sys

GENOME_PATHS = {
    "GRCh38": "https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v3/hg38.tsv",
    "mm10": "https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v3/mm10.tsv"
}

TEMPLATE = {
    "atac.pipeline_type" : "atac",
    "atac.align_only" : False,
    "atac.paired_end" : True,
    "atac.multimapping" : 4
}

def build_caper_config(sample, genome, bucket, read_length, out_path):
    config = TEMPLATE.copy()
    config["atac.title"] = sample
    config["atac.nodup_bams"] = [f"gs://{bucket}/results/{sample}/filtering/filtered.bam"]
    config["atac.genome_tsv"] = GENOME_PATHS[genome]
    config["atac.read_len"] = [read_length]
    
    with open(out_path, "w") as f:
        json.dump(config, f, indent=4)

sample = snakemake.params["sample"]
genome = snakemake.params["genome"]
bucket = snakemake.params["bucket"]
read_length = snakemake.params["read_length"]
sample_data, = snakemake.input
out_path, = snakemake.output
