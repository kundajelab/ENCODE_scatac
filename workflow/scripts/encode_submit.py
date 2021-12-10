import os
import json
from collections import Mapping
import encode_utils as eu
from encode_utils.connection import Connection

def retype_href(href, type):
    data = href.split(";", 1)[1]
    return f"data:{type};{data}"


def apply_patches(metadata, schema):
    if schema == "sc_atac_alignment_quality_metric":
        if "mito_stats" in metadata:
            type = "text/plain"
            href = metadata["mito_stats"]["href"]
            href_new = retype_href(href, type)
            entry = {
                "download": "frac_mito.tsv",
                "href": href_new,
                "type": type
            }
            metadata["mito_stats"] = entry

        alias_entries = metadata["aliases"][0].split("$")
        if alias_entries[-1] == "alignments_lib_comp_qc_metadata.json":
            alias_new = "$".join(alias_entries[:-1] + ["alignments_filtered_qc_metadata.json"])
            metadata["aliases"] = [alias_new]

    elif schema == "sc_atac_library_complexity_quality_metric":
        metadata.pop("positions_with_two_reads")

    elif schema == "sc_atac_multiplet_quality_metric":
        if "barcode_pairs_expanded" in metadata:
            type = "application/gzip"
            href = metadata["barcode_pairs_expanded"]["href"]
            href_new = retype_href(href, type)
            entry = {
                "download": "barcode_pairs_expanded.tsv.gz",
                "href": href_new,
                "type": type
            }
            metadata["barcode_pairs_expanded"] = entry

    elif schema == "file": # ALT
        alias_entries = metadata["aliases"][0].split("$")
        if alias_entries[-1] == "filtered.bam":
            deriv_from_new = "$".join(alias_entries[:-1] + ["raw.bam"])
            metadata["derived_from"] = [deriv_from_new]

    elif schema in ["sc_atac_analysis_quality_metric", "sc_atac_counts_summary_quality_metric"]: # ALT
        metric_of_entries = metadata["quality_metric_of"][0].split("$")
        metric_of_new = "$".join(metric_of_entries[:-1] + ["archr_project.tar.gz"])
        metadata["quality_metric_of"] = [metric_of_new]


def set_attachments(conn, payload):
    attachments = []
    for key, val in payload.items():
        if isinstance(val, Mapping) and ("path" in val):
            attachments.append(key)

    for key in attachments:
        val = payload[key]
        attachment = conn.set_attachment(document=val["path"])
        payload[key] = attachment

metadata_file = snakemake.input["json"]
dcc_mode = snakemake.config["dcc_mode"]
schema = snakemake.params["schema"]
step_run = snakemake.params["step_run"]
log_dir, = snakemake.log

os.environ["DCC_API_KEY"] = snakemake.params["dcc_api_key"]
os.environ["DCC_SECRET_KEY"] = snakemake.params["dcc_secret_key"]

with open(metadata_file) as f:
    metadata_raw = json.load(f)
metadata = {k: v for k, v in metadata_raw.items() if not k.startswith("_")}
metadata["step_run"] = step_run

eu.connection.LOG_DIR = log_dir

metadata[Connection.PROFILE_KEY] = schema
conn = Connection(dcc_mode, submission=True)
set_attachments(conn, metadata)
apply_patches(metadata, schema)
conn.post(metadata)





