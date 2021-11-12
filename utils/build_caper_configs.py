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
    "atac.true_rep_only" : True,
    "atac.paired_end" : True,
    "atac.multimapping" : 4
}

def build_caper_config(sample, genome, bucket, out_dir):
    config = TEMPLATE.copy()
    config["atac.title"] = sample
    config["atac.nodup_bams"] = [f"gs://{bucket}/results/{sample}/filtering/filtered.bam"]
    config["atac.genome_tsv"] = GENOME_PATHS[genome]
    
    with open(os.path.join(out_dir, f"{sample}.json"), "w") as f:
        json.dump(config, f, indent=4)

def load_samples(sample_path):
    with open(sample_path) as sample_file:
        h = sample_file.readline().rstrip('\n').split('\t')
        exp_ind = h.index("Experiment")
        rep_ind = h.index("Replicate")
        mod_ind = h.index("Modality")
        gen_ind = h.index("Genome")
        samples = []
        for line in sample_file:
            if line.startswith("#"):
                continue
            if line.startswith("@"):
                continue
            if line.startswith("$"):
                continue
            entries = line.rstrip('\n').split('\t')
            exp = entries[exp_ind]
            rep = int(entries[rep_ind])
            mod = entries[mod_ind]
            gen = entries[gen_ind]
            data = [
                exp,
                rep,
                mod,
                gen
            ]
            samples.append(data)
    return samples

def read_completed(completed_path):
    completed = set()
    with open(completed_path) as f:
        for line in f:
            completed.add(line)
    return completed

def build_configs(sample_path, bucket, out_dir, completed_path):
    os.makedirs(out_dir, exist_ok=True)
    completed = read_completed(completed_path)
    print(completed)
    samples = load_samples(sample_path)
    for s in samples:
        exp, rep, _, gen = s
        name = f"{exp}-{rep}"
        if name in completed:
            continue
        build_caper_config(name, gen, bucket, out_dir)

if __name__ == '__main__':
    sample_path = sys.argv[1]
    bucket = sys.argv[2]
    out_dir = sys.argv[3]
    completed = sys.argv[4]
    build_configs(sample_path, bucket, out_dir, completed)