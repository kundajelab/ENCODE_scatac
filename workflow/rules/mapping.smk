"""
Read mapping
"""

rule fetch_index:
    """
    Fetch Bowtie2 hg38 index
    """
    input:
        HTTP.remote(os.path.join(config["bwt2_idx"]), keep_local=True)
    output:
        idx = directory("bwt2_idx"),
        flag = touch(os.path.join("bwt2_idx", config["bwt2_idx_name"], config["bwt2_idx_prefix"]))
    conda:
        "../envs/mapping.yaml"
    group: 
        "mapping"
    shell:
        "unzip {input} -d {output.idx}"

rule bowtie2:
    """
    Read mapping (Bowtie2 aligner)
    """
    input:
        fastq1 = "results/{sample}/fastqs/R1_trim.fastq.gz",
        fastq2 = "results/{sample}/fastqs/R2_trim.fastq.gz",
        idx = os.path.join("bwt2_idx", config["bwt2_idx_name"], config["bwt2_idx_prefix"])
    output:
        bam_raw = pipe("temp/{sample}/mapping/raw.bam"),
    params:
        k = 1 + config["multimapping"],
    log:
        "logs/{sample}/mapping/bwt2.log"
    threads:
        max_threads
    conda:
        "../envs/mapping.yaml"
    group: 
        "mapping"
    shell:
        "bowtie2 -X 2000 --threads {threads} -x {input.idx} "
        "-1 {input.fastq1} -2 {input.fastq2} --sam-append-comment -k {params.k} 2> {log} | "
        "samtools view -u -S -o {output.bam_raw} -"

rule filter_multimappers:
    """
    Remove multimapping reads above threshold
    """
    input:
        "temp/{sample}/mapping/raw.bam"
    output:
        pipe("temp/{sample}/mapping/de-multimap.bam")
    params:
        multimapping = config["multimapping"],
        mmp_path = script_path("scripts/assign_multimappers.py")
    conda:
        "../envs/mapping.yaml"
    group: 
        "mapping"
    shell:
        "samtools view -h -f 2 {input} | "
        "python {params.mmp_path} --paired-end -k {params.multimapping} | "
        "samtools view -u -o {output} - "

rule sort_alignments:
    """
    Sort and index alignments
    """
    input: 
        "temp/{sample}/mapping/de-multimap.bam"
    output: 
        "results/{sample}/mapping/raw.bam"
    log:
        "logs/{sample}/mapping/sort.log"
    threads:
        max_threads
    conda:
        "../envs/mapping.yaml"
    shadow: 
        "minimal"
    group: 
        "mapping"
    shell: 
        "samtools sort -T . -@ {threads} -o {output} {input} 2> {log}; "
        "samtools index {output.bam};"

rule index_bam_raw:
    """
    Index raw BAM
    """
    input: 
        "results/{sample}/mapping/raw.bam"
    output: 
        "results/{sample}/mapping/raw.bam.bai"
    conda:
        "../envs/mapping.yaml"
    group: 
        "mapping"
    shell: 
        "samtools index {input}"

rule samstats_raw:
    """
    Run SAMStats on raw alignments
    """
    input:
        "results/{sample}/mapping/raw.bam"
    output:
        "results/{sample}/mapping/samstats_raw.txt"
    log:
        "logs/{sample}/mapping/samstats_raw.log"
    threads:
        max_threads
    conda:
        "../envs/filtering.yaml"
    group: 
        "mapping"
    shadow: 
        "minimal"
    shell:
        "samtools sort -T . -n -@ {threads} -O SAM {input} | " 
        "SAMstats --sorted_sam_file -  --outf {output} > {log}"

rule metadata_bam_raw:
    """
    Write raw BAM metadata
    """
    input: 
        bam = "results/{sample}/mapping/raw.bam",
        fq_R1 = "results/{sample}/fastqs/R1_trim.fastq.gz",
        fq_R2 = "results/{sample}/fastqs/R2_trim.fastq.gz"
    output: 
        "results/{sample}/mapping/raw_bam_metadata.json"
    params:
        output_group = "mapping",
        sample_data = lambda w: samples[w.sample]
    conda:
        "../envs/mapping.yaml"
    group: 
        "mapping"
    script: 
        "../scripts/write_file_metadata.py"

rule mapping_done:
    """
    Touch flag file upon group completion
    """
    input: 
        "results/{sample}/mapping/raw.bam",
        "results/{sample}/mapping/raw.bam.bai", 
        "results/{sample}/mapping/samstats_raw.txt", 
        "results/{sample}/mapping/raw_bam_metadata.json"
    output:
        touch("results/{sample}/mapping/mapping.done")
    group: 
        "mapping"