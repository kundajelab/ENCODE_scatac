import os
import json
import encode_utils as eu
from encode_utils.connection import Connection

metadata_file = snakemake.input["json"]
dcc_mode = snakemake.config["dcc_mode"]
schema = snakemake.params["schema"]
log_dir, = snakemake.log

os.environ["DCC_API_KEY"] = snakemake.params["dcc_api_key"]
os.environ["DCC_SECRET_KEY"] = snakemake.params["dcc_secret_key"]

with open(metadata_file) as f:
    metadata = json.load(f)

eu.connection.LOG_DIR = log_dir

metadata[Connection.PROFILE_KEY] = schema
conn = Connection(dcc_mode, submission=True)
conn.post(metadata)





