"""
FASTQ processing
"""

def get_fastq_inputs(w):
    return HTTP.remote(sample_data[w.sample]["fastq"][w.read][0], keep_local=config["keep_inputs"], username=os.environ["DCC_API_KEY"], password=os.environ["DCC_SECRET_KEY"]) 

rule strip_fastq:
    """
    Strip FASTQ read descriptions
    """
    input:
        get_fastq_inputs 
    output:
        pipe("temp/{sample}/fastqs/stripped_{read}.fastq")
    conda:
        "../envs/fastqs.yaml"
    group: 
        "fastqs"
    shell:
        "zcat {input} | sed 's/ .*//' > {output}"

rule detect_revcomp:
    """
    Detect whether to reverse complement barcodes
    """
    input:
        lambda w: HTTP.remote(sample_data[w.sample]["fastq"]["BC"][0], keep_local=config["keep_inputs"], username=os.environ["DCC_API_KEY"], password=os.environ["DCC_SECRET_KEY"]) 
    output:
        out = temp("temp/{sample}/fastqs/revcomp_indicator.txt"),
        qc = temp("temp/{sample}/fastqs/barcode_revcomp_full.txt")
    conda:
        "../envs/fastqs.yaml"
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
        whitelist = lambda w: config["bc_whitelist"][sample_config[w.sample]['modality']],
        revcomp = "temp/{sample}/fastqs/revcomp_indicator.txt"
    output: 
        fastq1_bc = pipe("temp/{sample}/fastqs/R1_bc_full.fastq"),
        fastq2_bc = pipe("temp/{sample}/fastqs/R2_bc_full.fastq"),
        qc_matching = temp("temp/{sample}/fastqs/barcode_matching_full.tsv")
    params:
        barcode_dist = lambda w: config["max_barcode_dist"],
        modality = lambda w: sample_data[w.sample]["modality"]
    threads:
        max_threads
    conda:
        "../envs/fastqs.yaml"
    group: 
        "fastqs"
    script:
        "scripts/match_barcodes.py"

rule fetch_ren:
    """
    Entry point for Bing Ren data
    """
    input:
        "temp/{sample}/fastqs/stripped_{read}.fastq"
    output:
        pipe("temp/{sample}/fastqs/{read}_bc_ren.fastq"),
    conda:
        "../envs/fastqs.yaml"
    group: 
        "fastqs"
    shell:
        "awk -v FS=':' "
        "'{if (NR%4==1) {s=\"@\"$2; for (i=3 ; i<=NF ; i++) {s = s \":\" $i } ; s = s \"\\tCB:Z:\" substr($1,2) ; print s} else {print $0}}' "
        "{input} > {output}"

rule dummy_qc_ren:
    """
    Create empty dummy QC files for Bing Ren data
    """
    output:
        touch("temp/{sample}/fastqs/barcode_revcomp_ren.txt"),
        touch("temp/{sample}/fastqs/barcode_matching_ren.tsv")
    group: 
        "fastqs"

rule move_fastq_qc:
    """
    Move QC files to final location
    """
    input:
        revcomp = lambda w: f"temp/{w.sample}/fastqs/barcode_revcomp_{'ren' if sample_config[w.sample]['modality'] == 'ren' else 'full'}.txt",
        qc_matching = lambda w: f"temp/{w.sample}/fastqs/barcode_matching_{'ren' if sample_config[w.sample]['modality'] == 'ren' else 'full'}.tsv"
    output:
        revcomp = "results/{sample}/fastqs/barcode_revcomp.txt",
        qc_matching = "results/{sample}/fastqs/barcode_matching.tsv"
    conda:
        "../envs/fastqs.yaml"
    group: 
        "fastqs"
    shell:
        "cp {input.revcomp} {output.revcomp}; "
        "cp {input.qc_matching} {output.qc_matching}"

rule trim_adapter:
    """
    Read adapter trimming
    """
    input:
        fastq1_bc = lambda w: f"temp/{w.sample}/fastqs/R1_bc_{'ren' if sample_config[w.sample]['modality'] == 'ren' else 'full'}.fastq",
        fastq2_bc = lambda w: f"temp/{w.sample}/fastqs/R2_bc_{'ren' if sample_config[w.sample]['modality'] == 'ren' else 'full'}.fastq"
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
        "../envs/fastqs.yaml"
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
        "../envs/fastqs.yaml"
    group: 
        "fastqs"
    script: 
        "../scripts/write_file_metadata.py"

rule metadata_qc_reads:
    """
    Write reads QC metadata
    """
    input: 
        data_file = "results/{sample}/fastqs/R1_trim.fastq.gz",
        barcode_matching = "results/{sample}/fastqs/barcode_matching.tsv",
        adapter_trimming = "results/{sample}/fastqs/trim_adapters.txt",
        barcode_revcomp = "results/{sample}/fastqs/barcode_revcomp.txt"
    output: 
        read_stats = "results/{sample}/fastqs/reads_qc_metadata.json",
    params:
        output_group = "fastqs",
        sample_data = lambda w: sample_data[w.sample]
    conda:
        "../envs/fastqs.yaml"
    group: 
        "fastqs"
    script: 
        "../scripts/write_qc_metadata.py"

rule fastqs_done:
    """
    Touch flag file upon group completion
    """
    input: 
        "results/{sample}/fastqs/R1_trim.fastq.gz",
        "results/{sample}/fastqs/R2_trim.fastq.gz", 
        "results/{sample}/fastqs/barcode_revcomp.txt",
        "results/{sample}/fastqs/barcode_matching.tsv", 
        "results/{sample}/fastqs/trim_adapters.txt",
        "results/{sample}/fastqs/R1_trim_metadata.json", 
        "results/{sample}/fastqs/R2_trim_metadata.json"
    output:
        touch("results/{sample}/fastqs/fastqs.done")
    group: 
        "fastqs"