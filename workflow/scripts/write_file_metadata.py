from collections import OrderedDict
import json
import os


def file_header(sample_data, config, out_path, preds, step_run, parse_preds=True):
    lab = config["dcc_lab"]
    experiment = sample_data["experiment"]
    replicate = sample_data["replicate_num"]
    alias = f"{lab}:{experiment}${replicate}${os.path.basename(out_path)}"
    if parse_preds:
        pred_ids = [f"{lab}:{experiment}${replicate}${os.path.basename(p)}" for p in preds]
    else:
        pred_ids = preds
    h = OrderedDict({
        "lab": lab,
        "award": config["dcc_award"],
        "dataset": experiment,
        "aliases": [alias],
        "submitted_file_name": os.path.abspath(out_path),
        "derived_from": pred_ids,
        "step_run": step_run
    })
    return h

def fastq_metadata(sample_data, pair, other):
    d = OrderedDict({
        "file_format": "fastq",
        "run_type": "paired-ended",
        "output_type": "reads",
        "platform": sample_data["platform"],
        "read_length": sample_data["read_length"],
        "replicate": sample_data["replicate_id"],
        "paired_end": pair,
    })
    if pair == "2":
        d["paired_with"] = other
    return d

def bam_metadata(sample_data):
    d = OrderedDict({
        "file_format": "bam",
        "assembly": "GRCh38",
        "mapped_run_type": "paired-ended",
        "output_type": "alignments",
        "mapped_read_length": sample_data["read_length"],
    })
    return d

def fragments_metadata(sample_data):
    d = OrderedDict({
        "file_format": "tar",
        "assembly": "GRCh38",
        "output_type": "fragments",
    })
    return d

def analyses_metadata(sample_data): #TODO
    d = OrderedDict({
        "file_format": "tar",
        "assembly": "GRCh38",
    })
    return d

def write_json(data, out_path):
    with open(out_path, "w") as f:
        json.dump(data, f, indent=4)

try:
    out_group = snakemake.params['output_group']
    sample_data = snakemake.params['sample_data']
    config = snakemake.config

    if out_group == "fastqs":
        step_run = "7f3f3341-e03f-40ce-b962-44851b80aa88" #TODO Replace with final value

        r1 = snakemake.input['r1']
        r2 = snakemake.input['r2']
        out1 = snakemake.output['r1']
        out2 = snakemake.output['r2']

        preds = list(sample_data["accessions"].values())
        
        h1 = file_header(sample_data, config, r1, preds, step_run, parse_preds=False)
        h2 = file_header(sample_data, config, r2, preds, step_run, parse_preds=False)

        d1 = fastq_metadata(sample_data, "1", h2["aliases"][0])
        d2 = fastq_metadata(sample_data, "2", h1["aliases"][0])
        
        s1 = h1 | d1
        s2 = h2 | d2

        write_json(s1, out1)
        write_json(s2, out2)

    elif out_group == "mapping":
        step_run = "7f3f3341-e03f-40ce-b962-44851b80aa88" #TODO Replace with final value

        bam = snakemake.input['bam']
        r1_pred = snakemake.input['fq_R1']
        r2_pred = snakemake.input['fq_R2']
        out, = snakemake.output

        h = file_header(sample_data, config, bam, [r1_pred, r2_pred], step_run)
        d = bam_metadata(sample_data)
        s = h | d

        write_json(s, out)

    elif out_group == "filtering":
        step_run = "7f3f3341-e03f-40ce-b962-44851b80aa88" #TODO Replace with final value

        bam = snakemake.input['bam']
        r1_pred = snakemake.input['fq_R1']
        r2_pred = snakemake.input['fq_R2']
        out, = snakemake.output

        h = file_header(sample_data, config, bam, [r1_pred, r2_pred], step_run)
        d = bam_metadata(sample_data)
        s = h | d

        write_json(s, out)

    elif out_group == "fragments":
        step_run = "7f3f3341-e03f-40ce-b962-44851b80aa88" #TODO Replace with final value

        fragments = snakemake.input['fragments']
        pred = snakemake.input['bam']
        out, = snakemake.output

        h = file_header(sample_data, config, fragments, [pred], step_run)
        d = fragments_metadata(sample_data)
        s = h | d

        write_json(s, out)

    elif out_group == "analyses":
        step_run = "7f3f3341-e03f-40ce-b962-44851b80aa88" #TODO Replace with final value

        tarball = snakemake.input['archr']
        pred = snakemake.input['fragments']
        out, = snakemake.output

        h = file_header(sample_data, config, tarball, [pred], step_run)
        s = h 

        write_json(s, out)


except NameError:
    pass