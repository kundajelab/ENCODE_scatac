"""
FASTQ processing
"""

rule strip_fastq:
    """
    Strip FASTQ read descriptions
    """
    input:
        lambda w: S3.remote(sample_data[w.sample]["fastq"][w.read], keep_local=config["keep_inputs"]) 
    output:
        pipe("temp/{sample}/fastqs/stripped_{read}.fastq")
    conda:
        "envs/fastqs.yaml"
    group: 
        "fastqs"
    shell:
        "zcat {input} | sed 's/ .*//' > {output}"

rule detect_revcomp:
    """
    Detect whether to reverse complement barcodes
    """
    input:
        lambda w: S3.remote(sample_data[w.sample]["fastq"]["BC"], keep_local=config["keep_inputs"]) 
    output:
        out = temp("temp/{sample}/fastqs/bc_revcomp.txt"),
        qc = "results/{sample}/fastqs/barcode_revcomp.txt"
    conda:
        "envs/fastqs.yaml"
    group: 
        "fastqs"
    script:
        "scripts/barcode_revcomp_detection.py"

rule match_barcodes: 
    """
    Barcode correction and filtering
    """
    input: 
        fq_R1 = "temp/{sample}/fastqs/stripped_R1.fastq",
        fq_R2 = "temp/{sample}/fastqs/stripped_R2.fastq",
        fq_BC = "temp/{sample}/fastqs/stripped_BC.fastq",
        whitelist = lambda w: config[w.modality]["bc_whitelist"],
        revcomp = "temp/{sample}/fastqs/bc_revcomp.txt"
    output: 
        fastq1_bc = pipe("temp/{sample}/fastqs/R1_bc_full.fastq"),
        fastq2_bc = pipe("temp/{sample}/fastqs/R2_bc_full.fastq"),
        qc_matching = "results/{sample}/fastqs/barcode_matching.tsv"
    params:
        barcode_dist = lambda w: config["max_barcode_dist"],
        modality = lambda w: sample_data[w.sample]["modality"]
    threads:
        max_threads
    conda:
        "envs/fastqs.yaml"
    group: 
        "fastqs"
    script:
        "scripts/match_barcodes.py"

rule fetch_ren:
    """
    Entry point for Bing Ren data
    """
    input:
        lambda w: S3.remote(sample_data[w.sample]["fastq"][w.read], keep_local=config["keep_inputs"]) 
    output:
        fastq = pipe("temp/{sample}/fastqs/{read}_bc_ren.fastq"),
        revcomp = touch("results/{sample}/fastqs/barcode_revcomp.txt"),
        qc_matching = touch("results/{sample}/fastqs/barcode_matching.tsv")
    conda:
        "envs/fastqs.yaml"
    group: 
        "fastqs"
    shell:
        "zcat {input} | sed 's/ .*//' > {output.fastq}"

rule trim_adapter:
    """
    Read adapter trimming
    """
    input:
        fastq1_bc = lambda w: f"temp/{w.sample}/fastqs/R1_bc_{'ren' if sample_data[w.sample]['technology'] == 'ren' else 'full'}.fastq",
        fastq2_bc = lambda w: f"temp/{w.sample}/fastqs/R2_bc_{'ren' if sample_data[w.sample]['technology'] == 'ren' else 'full'}.fastq"
    output:
        fastq1_trim = "results/{sample}/fastqs/R1_trim.fastq.gz",
        fastq2_trim = "results/{sample}/fastqs/R2_trim.fastq.gz",
        stats = "results/{sample}/fastqs/trim_adapters.txt"
    log:
        html = "logs/{sample}/fastqs/fastp.html",
        json = "logs/{sample}/fastqs/fastp.json"
    threads:
        max_threads
    conda:
        "envs/fastqs.yaml"
    group: 
        "fastqs"
    shell:
        "fastp -i {input.fastq1_bc} -I {input.fastq2_bc} -o {output.fastq1_trim} -O {output.fastq2_trim}"
        " -h {log.html} -j {log.json} -G -Q -L -w {threads} 2> {output.stats}"

rule metadata_fastq:
    """
    Write FASTQ metadata
    """
    input: 
        r1 = "results/{sample}/fastqs/R1_trim.fastq.gz",
        r2 = "results/{sample}/fastqs/R2_trim.fastq.gz"
    output: 
        r1 = "results/{sample}/fastqs/R1_trim_metadata.json",
        r2 = "results/{sample}/fastqs/R2_trim_metadata.json"
    params:
        output_group = "fastqs",
        sample_data = lambda w: sample_data[w.sample]
    conda:
        "envs/fastqs.yaml"
    group: 
        "fastqs"
    script: 
        "scripts/write_file_metadata.py"

rule fastqs_done:
    """
    Touch flag file upon group completion
    """
    input: 
        "results/{sample}/fastqs/R1_trim.fastq.gz",
        "results/{sample}/fastqs/R2_trim.fastq.gz", 
        "results/{sample}/fastqs/barcode_matching.tsv", 
        "results/{sample}/fastqs/trim_adapters.txt",
        "results/{sample}/fastqs/R1_trim_metadata.json", 
        "results/{sample}/fastqs/R2_trim_metadata.json"
    output:
        touch("results/{sample}/fastqs/fastqs.done")
    group: 
        "fastqs"