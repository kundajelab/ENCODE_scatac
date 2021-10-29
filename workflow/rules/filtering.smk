
"""
Alignment filtering
"""

rule filter_mito:
    """
    Filter and count mitochondrial reads
    """
    input: 
        "results/{sample}/mapping/raw_collated.bam",
    output: 
        bam = temp("temp/{sample}/filtering/no_mito.bam"),
        qc = "results/{sample}/filtering/frac_mito.tsv"
    conda:
        "../envs/filtering.yaml"
    group: 
        "filtering"
    script:
        "../scripts/filter_mito.py"

rule filter_multimappers:
    """
    Remove or assign multimapping reads
    """
    input:
        "temp/{sample}/filtering/no_mito.bam"
    output:
        temp("temp/{sample}/filtering/primary_align.bam")
    params:
        multimapping = config["multimapping"],
        mmp_path = script_path("scripts/assign_multimappers.py")
    conda:
        "../envs/filtering.yaml"
    group: 
        "filtering"
    shell:
        "samtools view -F 524 -f 2 -h {input} | "
        "python {params.mmp_path} --paired-end -k {params.multimapping} | "
        "samtools view -u - | "
        "samtools fixmate -r - {output}"

rule remove_duplicates:
    """
    Mark and remove PCR duplicates
    """
    input:
        "temp/{sample}/filtering/primary_align.bam"
    output:
        bam_markdup = temp("temp/{sample}/filtering/markdup.bam"),
        bam_nodup = temp("temp/{sample}/filtering/dedup.bam"),
        markdup_stats = "results/{sample}/filtering/markdup.txt"
    log:
        "logs/{sample}/filtering/picard.log"
    conda:
        "../envs/filtering.yaml"
    group: 
        "filtering"
    shell:
        "picard MarkDuplicates --INPUT {input} --OUTPUT /dev/stdout --METRICS_FILE {output.markdup_stats} "
        "--VALIDATION_STRINGENCY LENIENT --ASSUME_SORT_ORDER queryname --REMOVE_DUPLICATES false --BARCODE_TAG CB 2> {log} | "
        "tee {output.bam_markdup} | "
        "samtools view -F 1804 -f 2 -b -o {output.bam_nodup} -"

rule library_complexity:
    """
    Calculate PBC library complexity stats
    """
    input:
        "temp/{sample}/filtering/markdup.bam"
    output:
        "results/{sample}/filtering/pbc_stats.tsv",
    params:
        pbc_script = srcdir("../scripts/pbc_stats.py")
    conda:
        "../envs/filtering.yaml"
    group: 
        "filtering"
    shell:
        "samtools view {input} | python {params.pbc_script} {output}"

rule sort_filtered_alignments:
    """
    Sort filtered alignments
    """
    input: 
        "temp/{sample}/filtering/dedup.bam"
    output: 
        "results/{sample}/filtering/filtered.bam"
    log:
        "logs/{sample}/filtering/sort.log"
    threads:
        max_threads
    conda:
        "../envs/filtering.yaml"
    group: 
        "filtering"
    shell: 
        "samtools sort -T . -@ {threads} -o {output} {input} 2> {log} "

rule index_bam_filtered:
    """
    Index filtered BAM
    """
    input: 
        "results/{sample}/filtering/filtered.bam"
    output: 
        "results/{sample}/filtering/filtered.bam.bai"
    conda:
        "../envs/filtering.yaml"
    group: 
        "filtering"
    shell: 
        "samtools index {input}"

rule samstats_filtered:
    """
    SAMstats on filtered alignments
    """
    input:
        bam = "temp/{sample}/filtering/dedup.bam",
        dep = "results/{sample}/filtering/filtered.bam"
    output:
        "results/{sample}/filtering/samstats_filtered.txt"
    log:
        "logs/{sample}/filtering/samstats_filtered.log"
    conda:
        "../envs/filtering.yaml"
    group: 
        "filtering"
    shell:
        "samtools view -o - {input.bam} | " 
        "SAMstats --sorted_sam_file -  --outf {output} > {log}"

rule filtering_done:
    """
    Touch flag file upon group completion
    """
    input: 
        "results/{sample}/filtering/filtered.bam",
        "results/{sample}/filtering/filtered.bam.bai", 
        "results/{sample}/filtering/frac_mito.tsv", 
        "results/{sample}/filtering/markdup.txt",
        "results/{sample}/filtering/pbc_stats.tsv",
        "results/{sample}/filtering/samstats_filtered.txt",
    output:
        touch("results/{sample}/filtering/filtering.done")
    # group: 
    #     "filtering"