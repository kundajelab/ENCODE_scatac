
"""
Alignment filtering
"""

rule filter_mito:
    """
    Filter and count mitochondrial reads
    """
    input: 
        "results/{sample}/mapping/raw_unsorted.bam",
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
        "samtools view -F 524 -f 2 -u {input} |"
        "python {params.mmp_path} --paired-end -k {params.multimapping} | "
        "samtools fixmate -r - {output}"

rule sort_filtered_alignments:
    """
    Sort filtered alignments
    """
    input: 
        "temp/{sample}/filtering/primary_align.bam"
    output: 
        "temp/{sample}/filtering/primary_align_sorted.bam"
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

rule remove_duplicates:
    """
    Mark and remove PCR duplicates
    """
    input:
        "temp/{sample}/filtering/primary_align_sorted.bam"
    output:
        bam_markdup = temp("temp/{sample}/filtering/markdup.bam"),
        bam_nodup = "results/{sample}/filtering/filtered.bam",
        markdup_stats = "results/{sample}/filtering/markdup.txt"
    log:
        "logs/{sample}/filtering/picard.log"
    conda:
        "../envs/filtering.yaml"
    group: 
        "filtering"
    shell:
        "picard MarkDuplicates --INPUT {input} --OUTPUT /dev/stdout --METRICS_FILE {output.markdup_stats} "
        "--VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --REMOVE_DUPLICATES false --BARCODE_TAG CB 2> {log} | "
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
        "results/{sample}/filtering/filtered.bam"
    output:
        "results/{sample}/filtering/samstats_filtered.txt"
    log:
        "logs/{sample}/filtering/samstats_filtered.log"
    threads:
        max_threads
    conda:
        "../envs/filtering.yaml"
    group: 
        "filtering"
    # shadow: 
    #     "minimal"
    shell:
        "samtools sort -T . -n -@ {threads} -O SAM {input} | " 
        "SAMstats --sorted_sam_file -  --outf {output} > {log}"

rule metadata_bam_filtered:
    """
    Write filtered BAM metadata
    """
    input: 
        bam = "results/{sample}/filtering/filtered.bam",
        fq_R1 = "results/{sample}/fastqs/R1_trim.fastq.gz", # attach to fastqs
        fq_R2 = "results/{sample}/fastqs/R2_trim.fastq.gz",
        input_data = "results/{sample}/input_data.json"
    output: 
        "results/{sample}/filtering/filtered_bam_metadata.json",
    params:
        output_group = "filtering",
        sample_data = lambda w: sample_data[w.sample]
    conda:
        "../envs/filtering.yaml"
    group: 
        "filtering"
    script: 
        "../scripts/write_file_metadata.py"

rule metadata_qc_alignments_filtered:
    """
    Write filtered alignments qc metadata
    """
    input: 
        data_file = "results/{sample}/filtering/filtered.bam",
        samstats_filtered = "results/{sample}/filtering/samstats_filtered.txt",
        picard_markdup = "results/{sample}/filtering/markdup.txt",
        pbc_stats = "results/{sample}/filtering/pbc_stats.tsv",
        frac_mito = "results/{sample}/filtering/frac_mito.tsv",
        input_data = "results/{sample}/input_data.json"
    output: 
        alignment_stats = "results/{sample}/filtering/alignments_filtered_qc_metadata.json",
        lib_comp_stats = "results/{sample}/filtering/alignments_lib_comp_qc_metadata.json"
    params:
        output_group = "filtering",
        sample_data = lambda w: sample_data[w.sample]
    conda:
        "../envs/filtering.yaml"
    group: 
        "filtering"
    script: 
        "../scripts/write_qc_metadata.py"

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
        "results/{sample}/filtering/filtered_bam_metadata.json",
        "results/{sample}/filtering/alignments_filtered_qc_metadata.json",
        "results/{sample}/filtering/alignments_lib_comp_qc_metadata.json"
    output:
        touch("results/{sample}/filtering/filtering.done")
    group: 
        "filtering"