import os

import encode_utils as eu
from encode_utils.connection import Connection

resubmit_files = snakemake.input
dcc_mode = snakemake.config["dcc_mode"]
log_dir, = snakemake.log

os.environ["DCC_API_KEY"] = snakemake.params["dcc_api_key"]
os.environ["DCC_SECRET_KEY"] = snakemake.params["dcc_secret_key"]

eu.connection.LOG_DIR = log_dir

conn = Connection(dcc_mode, submission=True)

for k, v in resubmit_files.items():
    conn.upload_file(file_id=k ,file_path=v)





