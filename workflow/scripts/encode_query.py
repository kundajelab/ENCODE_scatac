import os
import json
import encode_utils as eu
from encode_utils.connection import Connection

sample_data_file, = snakemake.output
sample = snakemake.params["sample"]
experiment = snakemake.params["experiment"]
replicate_num = snakemake.params["replicate"]
modality = snakemake.params["modality"]
log_dir, = snakemake.log

os.environ["DCC_API_KEY"] = snakemake.params["dcc_api_key"]
os.environ["DCC_SECRET_KEY"] = snakemake.params["dcc_secret_key"]

eu.connection.LOG_DIR = log_dir

conn = Connection("prod")
data = conn.get(sample)

r1 = {}
r2 = {}
bc = {}

replicate_id = None
for rep in data["replicates"]:
    if rep["biological_replicate_number"] == replicate_num:
        replicate_id = rep["uuid"]

platform = None
read_length = None
files = data["files"]
for f in files:
    id = f["@id"]
    if f["file_format"] != "fastq":
        continue
    if f["replicate"]["biological_replicate_number"] != replicate_num:
        continue

    p = f["platform"]["uuid"]
    if platform is not None and p != platform:
        raise ValueError("Multiple sequencing platforms detected in input")

    if f["output_type"] == "index reads":
        bc[id] = f
        continue
    
    l = f["read_length"]
    if read_length is not None and l != read_length:
        raise ValueError("Multiple read lengths detected in input")

    if f["paired_end"] == "1":
        r1[id] = f
    elif f["paired_end"] == "2":
        r2[id] = f

out_data = {
    "experiment": experiment,
    "replicate_num": replicate_num,
    "replicate_id": replicate_id,
    "platform": platform,
    "read_length": read_length,
    "fastq": {"R1": [], "R2": [], "BC": []},
    "accessions": {"R1": [], "R2": [], "BC": []}
}

for f in bc.values():
    acc_bc = f["accession"]
    m0, m1 = f["index_of"]

    if m0 in r1 and m1 in r2:
        r1_fq = r1[m0]["s3_uri"]
        r2_fq = r2[m1]["s3_uri"]
        r1_acc = r1[m0]["accession"]
        r2_acc = r2[m1]["accession"]

        out_data["fastq"]["R1"].append(r1_fq)
        out_data["fastq"]["R2"].append(r2_fq)
        out_data["accessions"]["R1"].append(r1_acc)
        out_data["accessions"]["R2"].append(r2_acc)

    elif m1 in r1 and m0 in r2:
        r1_fq = r1[m1]["s3_uri"]
        r2_fq = r2[m0]["s3_uri"]
        r1_acc = r1[m1]["accession"]
        r2_acc = r2[m0]["accession"]

        out_data["fastq"]["R1"].append(r1_fq)
        out_data["fastq"]["R2"].append(r2_fq)
        out_data["accessions"]["R1"].append(r1_acc)
        out_data["accessions"]["R2"].append(r2_acc)

    else:
        raise ValueError("Index FASTQ does not properly match with reads")

    bc_fq = bc["s3_uri"]
    bc_acc = bc["accession"]
    out_data["fastq"]["BC"].append(bc_fq)
    out_data["accessions"]["R2"].append(r2_acc)

with open(sample_data_file) as f:
    metadata = json.dump(out_data, f, indent=4)


