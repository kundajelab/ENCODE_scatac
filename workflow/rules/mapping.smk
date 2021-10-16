"""
Read mapping
"""

def get_idx_prefix(w):
    genome = sample_config[w.sample]["genome"]
    entry = config["genome"][genome]
    return os.path.join("bwt2_idx", "unpacked", genome, entry["bwt2_idx_prefix"])

def get_idx_files(w):
    genome = sample_config[w.sample]["genome"]
    entry = config["genome"][genome]
    return [
        os.path.join("bwt2_idx", "unpacked", genome, f"{entry['bwt2_idx_prefix']}.{s}.bt2") 
        for s in ["1", "2", "3", "4", "rev.1", "rev.2"]
    ]

rule bowtie2:
    """
    Read mapping (Bowtie2 aligner)
    """
    input:
        fastq1 = "results/{sample}/fastqs/R1_trim.fastq.gz",
        fastq2 = "results/{sample}/fastqs/R2_trim.fastq.gz",
        prefix = get_idx_prefix,
        files = get_idx_files
    output:
        bam_raw = "results/{sample}/mapping/raw_unsorted.bam",
        qc = "results/{sample}/mapping/bwt2_stats.txt"
    params:
        k = 1 + config["multimapping"]
    threads:
        max_threads
    conda:
        "../envs/mapping.yaml"
    group: 
        "mapping"
    shell:
        "bowtie2 -X 2000 --threads {threads} -x {input.prefix} "
        "-1 {input.fastq1} -2 {input.fastq2} --sam-append-comment -k {params.k} 2> {output.qc} | "
        "samtools view -b -S -o {output.bam_raw} -"

rule sort_alignments:
    """
    Sort alignments
    """
    input: 
        "results/{sample}/mapping/raw_unsorted.bam"
    output: 
        "results/{sample}/mapping/raw.bam"
    log:
        "logs/{sample}/mapping/sort.log"
    threads:
        max_threads
    conda:
        "../envs/mapping.yaml"
    group: 
        "mapping"
    shell: 
        "samtools sort -T . -@ {threads} -o {output} {input} 2> {log} "

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
        bam = "results/{sample}/mapping/raw_unsorted.bam",
        dep = "results/{sample}/mapping/raw.bam" # Force dependency
    output:
        "results/{sample}/mapping/samstats_raw.txt"
    log:
        "logs/{sample}/mapping/samstats_raw.log"
    threads:
        max_threads
    conda:
        "../envs/mapping.yaml"
    group: 
        "mapping"
    shell:
        "samtools sort -T . -n -@ {threads} -O SAM {input.bam} | " 
        "SAMstats --sorted_sam_file -  --outf {output} > {log}"

rule mapping_done:
    """
    Touch flag file upon group completion
    """
    input: 
        "results/{sample}/mapping/raw.bam",
        "results/{sample}/mapping/raw.bam.bai", 
        "results/{sample}/mapping/samstats_raw.txt"
    output:
        touch("results/{sample}/mapping/mapping.done")
    group: 
        "mapping"