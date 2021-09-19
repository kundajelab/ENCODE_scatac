
"""
ENCODE portal submission
"""
# Note: these rules are not run by default

rule submit_fastq_1:
    """
    Submit FASTQ pair 1
    """
    input: 
        json = "results/{sample}/fastqs/R1_trim_metadata.json"
    output: 
        touch("results/{sample}/submit/R1_trim_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/R1_trim_submit")
    conda:
        "envs/portal.yaml"
    group: 
        "submit"
    script: 
        "scripts/encode_submit.py"

rule submit_fastq_2:
    """
    Submit FASTQ pair 2
    """
    input: 
        json = "results/{sample}/fastqs/R2_trim_metadata.json",
        fastq1 = "results/{sample}/submit/R1_trim_submit.done"
    output: 
        touch("results/{sample}/submit/R2_trim_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/R2_trim_submit")
    conda:
        "envs/portal.yaml"
    group: 
        "submit"
    script: 
        "scripts/encode_submit.py"

rule submit_bam_raw:
    """
    Submit raw BAM
    """
    input: 
        json = "results/{sample}/mapping/raw_bam_metadata.json",
        fastq1 = "results/{sample}/submit/R1_trim_submit.done",
        fastq2 = "results/{sample}/submit/R2_trim_submit.done"
    output: 
        touch("results/{sample}/submit/raw_bam_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/raw_bam_submit")
    conda:
        "envs/portal.yaml"
    group: 
        "submit"
    script: 
        "scripts/encode_submit.py"

rule submit_bam_filtered:
    """
    Submit filtered BAM
    """
    input: 
        json = "results/{sample}/filtering/filtered_bam_metadata.json",
        fastq1 = "results/{sample}/submit/R1_trim_submit.done",
        fastq2 = "results/{sample}/submit/R2_trim_submit.done"
    output: 
        touch("results/{sample}/submit/filtered_bam_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/filtered_bam_submit")
    conda:
        "envs/portal.yaml"
    group: 
        "submit"
    script: 
        "scripts/encode_submit.py"

rule submit_fragments:
    """
    Submit fragment file
    """
    input: 
        json = "results/{sample}/fragments/fragments_metadata.json",
        bam = "results/{sample}/submit/filtered_bam_submit.done"
    output: 
        touch("results/{sample}/submit/fragments_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/filtered_bam_submit")
    conda:
        "envs/portal.yaml"
    group: 
        "submit"
    script: 
        "scripts/encode_submit.py"

rule submit_done:
    """
    Touch flag file upon group completion
    """
    input: 
        "results/{sample}/submit/R1_trim_submit.done",
        "results/{sample}/submit/R2_trim_submit.done",
        "results/{sample}/submit/raw_bam_submit.done",
        "results/{sample}/submit/filtered_bam_submit.done",
        "results/{sample}/submit/fragments_submit.done"
    output:
        touch("results/{sample}/submit/submit.done")
    group: 
        "submit"
