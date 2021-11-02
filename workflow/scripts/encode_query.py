import os
import json
import statistics
from urllib.parse import urljoin
import encode_utils as eu
from encode_utils.connection import Connection

sample_data_file, = snakemake.output
dcc_mode = snakemake.config["dcc_mode"]
experiment = snakemake.params["experiment"]
replicate_num = snakemake.params["replicate"]
modality = snakemake.params["modality"]
assembly = snakemake.params["assembly"]
log_dir, = snakemake.log

os.environ["DCC_API_KEY"] = snakemake.params["dcc_api_key"]
os.environ["DCC_SECRET_KEY"] = snakemake.params["dcc_secret_key"]

eu.connection.LOG_DIR = log_dir

conn = Connection(dcc_mode)
server = conn.dcc_url
data = conn.get(experiment)

r1 = {}
r2 = {}
bc = {}

replicate_id = None
for rep in data["replicates"]:
    if rep["biological_replicate_number"] == replicate_num:
        replicate_id = rep["uuid"]

platform = None
read_lengths = []
files = data["files"]
for f in files:
    id = f["@id"]
    if f["file_format"] != "fastq":
        continue
    if f["replicate"]["biological_replicate_number"] != replicate_num:
        continue
    if "derived_from" in f:
        continue

    p = f["platform"]["uuid"]
    if platform is not None and p != platform:
        raise ValueError("Multiple sequencing platforms detected in input")
    platform = p

    if f["output_type"] == "index reads":
        bc[id] = f
        continue
    
    l = f["read_length"]
    read_lengths.append(l)

    if f["paired_end"] == "1":
        r1[id] = f
    elif f["paired_end"] == "2":
        r2[id] = f

if max(read_lengths) - min(read_lengths) > 4:
    raise ValueError("Inconsistent read lengths in input FASTQs")

read_length = statistics.median_low(read_lengths)

out_data = {
    "experiment": experiment,
    "replicate_num": replicate_num,
    "replicate_id": replicate_id,
    "modality": modality,
    "platform": platform,
    "read_length": read_length,
    "assembly": assembly
}

if modality == "ren":
    out_data |= {
        "fastq": {"R1": [], "R2": []},
        "accessions": {"R1": [], "R2": []}
    }
    
    for k, v in r1.items():
        r1_fq = urljoin(server, v["href"])
        r1_acc = v["accession"]

        p2 = v["paired_with"]
        r2_fq = urljoin(server, r2[p2]["href"])
        r2_acc = r2[p2]["accession"]

        out_data["fastq"]["R1"].append(r1_fq)
        out_data["fastq"]["R2"].append(r2_fq)
        out_data["accessions"]["R1"].append(r1_acc)
        out_data["accessions"]["R2"].append(r2_acc)

else:
    out_data |= {
        "fastq": {"R1": [], "R2": [], "BC": []},
        "accessions": {"R1": [], "R2": [], "BC": []}
    }

    for f in bc.values():
        m0, m1 = f["index_of"]

        if m0 in r1 and m1 in r2:
            r1_fq = urljoin(server, r1[m0]["href"])
            r2_fq = urljoin(server, r2[m1]["href"])
            r1_acc = r1[m0]["accession"]
            r2_acc = r2[m1]["accession"]

            out_data["fastq"]["R1"].append(r1_fq)
            out_data["fastq"]["R2"].append(r2_fq)
            out_data["accessions"]["R1"].append(r1_acc)
            out_data["accessions"]["R2"].append(r2_acc)

        elif m1 in r1 and m0 in r2:
            r1_fq = urljoin(server, r1[m1]["href"])
            r2_fq = urljoin(server, r2[m0]["href"])
            r1_acc = r1[m1]["accession"]
            r2_acc = r2[m0]["accession"]

            out_data["fastq"]["R1"].append(r1_fq)
            out_data["fastq"]["R2"].append(r2_fq)
            out_data["accessions"]["R1"].append(r1_acc)
            out_data["accessions"]["R2"].append(r2_acc)

        else:
            raise ValueError("Index FASTQ does not properly match with reads")
        
        bc_fq = urljoin(server, f["href"])
        bc_acc = f["accession"]
        out_data["fastq"]["BC"].append(bc_fq)
        out_data["accessions"]["BC"].append(bc_acc)

with open(sample_data_file, 'w') as f:
    metadata = json.dump(out_data, f, indent=4)


