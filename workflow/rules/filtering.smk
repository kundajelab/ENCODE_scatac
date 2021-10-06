
"""
Alignment filtering
"""

rule count_mito:
    """
    Count mitochondrial reads
    """
    input: 
        bam = "results/{sample}/mapping/raw.bam",
        bai = "results/{sample}/mapping/raw.bam.bai"
    output: 
        temp("temp/{sample}/filtering/count_mito.txt"),
    conda:
        "../envs/filtering.yaml"
    group: 
        "filtering"
    shell: 
        "samtools view -F 264 -c -u -o {output} {input.bam} chrM "

rule remove_mito:
    """
    Keep and count non-mitochondrial reads
    """
    input: 
        bam = "results/{sample}/mapping/raw.bam",
        bai = "results/{sample}/mapping/raw.bam.bai"
    output: 
        bam = temp("temp/{sample}/filtering/no_mito.bam"),
        count_no_mito = temp("temp/{sample}/filtering/count_no_mito.txt")
    conda:
        "../envs/filtering.yaml"
    group: 
        "filtering"
    shell: 
        "samtools idxstats {input.bam} | cut -f 1 | grep -v chrM | "
        "xargs samtools view -u {input.bam} | "
        "tee {output.bam} | samtools view -F 264 -c -o {output.count_no_mito} -"

rule frac_mito:
    """
    Calculate fraction of mitochondrial reads
    """
    input: 
        count_mito = "temp/{sample}/filtering/count_mito.txt",
        count_no_mito = "temp/{sample}/filtering/count_no_mito.txt"
    output: 
        "results/{sample}/filtering/frac_mito.tsv"
    conda:
        "../envs/filtering.yaml"
    group: 
        "filtering"
    shell: 
        "rm=$(<{input.count_mito}); "
        "rn=$(<{input.count_no_mito}); "
        "printf \"Non-Mitochondrial\\tMitochondrial\\n%d\\t%d\" \"$rn\" \"$rm\" > {output}"

rule assign_primary:
    """
    Assign multimapping reads to primary alignment
    """
    input:
        "temp/{sample}/filtering/no_mito.bam"
    output:
        temp("temp/{sample}/filtering/primary_align.bam")
    conda:
        "../envs/filtering.yaml"
    group: 
        "filtering"
    shell:
        "samtools view -F 1804 -b -o {output} {input} "

rule remove_duplicates:
    """
    Mark and remove PCR duplicates
    """
    input:
        "temp/{sample}/filtering/primary_align.bam"
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
        "picard MarkDuplicates --INPUT {input} --OUTPUT /dev/stdout --METRICS_FILE {output.markdup_stats} --COMPRESSION_LEVEL 0 "
        "--VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --REMOVE_DUPLICATES false --BARCODE_TAG CB 2> {log} | "
        "tee {output.bam_markdup} | "
        "samtools view -F 1024 -f 2 -b -o {output.bam_nodup} -"

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
    shadow: 
        "minimal"
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