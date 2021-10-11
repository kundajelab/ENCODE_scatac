
"""
ENCODE portal submission
"""
# Note: these rules are not run by default

rule submit_fastq_1:
    """
    Submit FASTQ pair 1
    """
    input: 
        json = "metadata/{sample}/R1_trim_metadata.json"
    output: 
        touch("submit/{sample}/R1_trim_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/R1_trim_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_fastq_2:
    """
    Submit FASTQ pair 2
    """
    input: 
        json = "metadata/{sample}/R2_trim_metadata.json",
        fastq1 = "submit/{sample}/R1_trim_submit.done"
    output: 
        touch("submit/{sample}/R2_trim_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/R2_trim_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_reads_qc: 
    """
    Submit reads QC
    """
    input: 
        json = "metadata/{sample}/reads_qc_metadata.json",
        fastq1 = "submit/{sample}/R1_trim_submit.done"
    output: 
        touch("submit/{sample}/reads_qc_metadata_submit.done")
    params:
        schema = "file", #TODO
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/reads_qc_metadata_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_bam_raw:
    """
    Submit raw BAM
    """
    input: 
        json = "metadata/{sample}/raw_bam_metadata.json",
        fastq1 = "submit/{sample}/R1_trim_submit.done",
        fastq2 = "submit/{sample}/R2_trim_submit.done"
    output: 
        touch("submit/{sample}/raw_bam_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/raw_bam_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_alignments_raw_qc: 
    """
    Submit raw alignments QC
    """
    input: 
        json = "metadata/{sample}/alignments_raw_qc_metadata.json",
        raw_bam = "submit/{sample}/raw_bam_submit.done"
    output: 
        touch("submit/{sample}/alignments_raw_qc_metadata_submit.done")
    params:
        schema = "file", #TODO
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/alignments_raw_qc_metadata_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_bam_filtered:
    """
    Submit filtered BAM
    """
    input: 
        json = "metadata/{sample}/filtering/filtered_bam_metadata.json",
        fastq1 = "submit/{sample}/R1_trim_submit.done",
        fastq2 = "submit/{sample}/R2_trim_submit.done"
    output: 
        touch("submit/{sample}/filtered_bam_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/filtered_bam_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_alignments_filtered_qc: 
    """
    Submit filtered alignments QC
    """
    input: 
        json = "metadata/{sample}/filtering/alignments_filtered_qc_metadata.json",
        filtered_bam = "submit/{sample}/filtered_bam_submit.done"
    output: 
        touch("submit/{sample}/alignments_filtered_qc_metadata_submit.done")
    params:
        schema = "file", #TODO
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/alignments_filtered_qc_metadata_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_lib_comp_qc: 
    """
    Submit library complexity QC
    """
    input: 
        json = "metadata/{sample}/filtering/alignments_lib_comp_qc_metadata.json",
        filtered_bam = "submit/{sample}/filtered_bam_submit.done"
    output: 
        touch("submit/{sample}/alignments_lib_comp_qc_metadata_submit.done")
    params:
        schema = "file", #TODO
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/alignments_lib_comp_qc_metadata_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_fragments:
    """
    Submit fragment file
    """
    input: 
        json = "metadata/{sample}/fragments/fragments_metadata.json",
        bam = "submit/{sample}/filtered_bam_submit.done"
    output: 
        touch("submit/{sample}/fragments_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/filtered_bam_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_done:
    """
    Touch flag file upon group completion
    """
    input: 
        "submit/{sample}/R1_trim_submit.done",
        "submit/{sample}/R2_trim_submit.done",
        "submit/{sample}/reads_qc_metadata_submit.done",
        "submit/{sample}/raw_bam_submit.done",
        "submit/{sample}/filtered_bam_submit.done",
        "submit/{sample}/fragments_submit.done"
    output:
        touch("submit/{sample}/submit.done")
    group: 
        "submit"
