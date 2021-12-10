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
            if "mito_stats" in metadata:
                metadata["barcode_pairs_expanded"] = entry


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





