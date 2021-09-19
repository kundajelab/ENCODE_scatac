"""
Fragment file generation
"""

rule bam_to_fragments: 
    """
    Convert BAM to fragment file
    """
    input:
        bam = "results/{sample}/filtering/filtered.bam",
        bai = "results/{sample}/filtering/filtered.bam.bai"
    output:
        pipe("temp/{sample}/fragments/fragments_unsorted.tsv")
    log:
        "logs/{sample}/fragments/sinto.log"
    threads:
        max_threads
    conda:
        "../envs/fragments.yaml"
    group: 
        "fragments"
    shell:
        "sinto fragments -b {input.bam} -f {output} " 
        "--shift_plus 4 --shift_minus -4 --min_mapq 0 "
        "--max_distance 2000 --min_distance 10 --barcodetag CB --nproc {threads} > {log}"

rule sort_fragments:
    """
    Sort and compress fragments
    """
    input: 
        "temp/{sample}/fragments/fragments_unsorted.tsv"
    output: 
        "results/{sample}/fragments/fragments.tsv.gz",
    threads: 
        max_threads
    conda:
        "../envs/fragments.yaml"
    group: 
        "fragments"
    shell: 
        "LC_ALL=C sort -k1,1 -k2,2n -k3,3n -k4,4 -t$'\\t' --parallel={threads} {input} | " # Sort the file by chr, start, end, then barcode_id
        "bgzip > {output}"

rule index_fragments:
    """
    Index fragments file
    """
    input: 
        "results/{sample}/fragments/fragments.tsv.gz"
    output: 
        "results/{sample}/fragments/fragments.tsv.gz.tbi"
    conda:
        "../envs/fragments.yaml"
    group: 
        "fragments"
    shell: 
        "tabix --zero-based --preset bed {input}"

rule tarball_fragments: 
    """
    Create fragments + index tarball
    """
    input:
        frag = "results/{sample}/fragments/fragments.tsv.gz",
        frag_ind = "results/{sample}/fragments/fragments.tsv.gz.tbi"
    output:
        "results/{sample}/fragments/fragments.tar.gz"
    conda:
        "../envs/fragments.yaml"
    group: 
        "fragments"
    shell:
        "tar -zcf {output} {input.frag} {input.frag_ind}"

rule metadata_fragments:
    """
    Write fragment file metadata
    """
    input: 
        fragments = "results/{sample}/fragments/fragments.tar.gz",
        bam = "results/{sample}/filtering/filtered.bam"
    output: 
        "results/{sample}/fragments/fragments_metadata.json"
    params:
        output_group = "fragments",
        sample_data = lambda w: samples[w.sample]
    conda:
        "../envs/fragments.yaml"
    group: 
        "fragments"
    script: 
        "../scripts/write_file_metadata.py"

rule fragments_done:
    """
    Touch flag file upon group completion
    """
    input: 
        "results/{sample}/fragments/fragments.tsv.gz",
        "results/{sample}/fragments/fragments.tsv.gz.tbi",
        "results/{sample}/fragments/fragments.tar.gz",
        "results/{sample}/fragments/fragments_metadata.json",
    output:
        touch("results/{sample}/fragments/fragments.done")
    group: 
        "fragments"