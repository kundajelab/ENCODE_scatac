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
        pipe("temp/{sample}/fragments/fragments_raw.tsv")
    params:
        shift_plus = config["tn5_shift_plus"],
        shift_minus = config["tn5_shift_minus"]
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
        "--shift_plus {params.shift_plus} --shift_minus {params.shift_minus} --min_mapq 0 "
        "--max_distance 2000 --min_distance 10 --barcodetag CB --nproc {threads} > {log}"

rule sort_fragments:
    """
    Sort and compress fragments
    """
    input: b
        "temp/{sample}/fragments/fragments_raw.tsv"
    output: 
        "results/{sample}/fragments/fragments.tsv.gz"
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

rule detect_multiplets:
    """
    Detect gel bead and barcode multiplets
    """
    input: 
        frag = "results/{sample}/fragments/fragments.tsv.gz",
        frag_ind = "results/{sample}/fragments/fragments.tsv.gz.tbi"
    output: 
        barcodes_strict = "results/{sample}/fragments/multiplet_barcodes_strict.tsv",
        barcodes_expanded = "results/{sample}/fragments/multiplet_barcodes_expanded.tsv",
        qc = "results/{sample}/fragments/multiplet_stats.txt"
    conda:
        "../envs/fragments.yaml"
    group: 
        "fragments"
    script: 
        "../scripts/detect_multiplets.py"

# rule sort_fragments_filtered:
#     """
#     Compress filtered fragments
#     """
#     input: 
#         "temp/{sample}/fragments/fragments_unsorted.tsv"
#     output: 
#         temp("temp/{sample}/fragments/fragments_filtered.tsv.gz")
#     threads: 
#         max_threads
#     conda:
#         "../envs/fragments.yaml"
#     group: 
#         "fragments"
#     shell: 
#         "bgzip {input} > {output}"

# rule index_fragments:
#     """
#     Index filtered fragments file
#     """
#     input: 
#         "temp/{sample}/fragments/fragments_filtered.tsv.gz"
#     output: 
#         temp("temp/{sample}/fragments/fragments_filtered.tsv.gz.tbi")
#     conda:
#         "../envs/fragments.yaml"
#     group: 
#         "fragments"
#     shell: 
#         "tabix --zero-based --preset bed {input}"

# rule output_fragments:
#     """
#     Move fragments file to final location
#     """
#     input: 
#         frag = lambda w: f"temp/{w.sample}/fragments/fragments_{'filtered' if sample_config[w.sample]['modality'] == '10x' else 'unfiltered'}.tsv.gz",
#         ind = lambda w: f"temp/{w.sample}/fragments/fragments_{'filtered' if sample_config[w.sample]['modality'] == '10x' else 'unfiltered'}.tsv.gz.tbi"
#     output: 
#         frag = "results/{sample}/fragments/fragments.tsv.gz",
#         ind = "results/{sample}/fragments/fragments.tsv.gz.tbi"
#     conda:
#         "../envs/fragments.yaml"
#     group: 
#         "fragments"
#     shell: 
#         "cp {input.frag} {output.frag}; "
#         "cp {input.ind} {output.ind}"

# rule placeholder_fragments_qc:
#     """
#     Create placeholder QC files if no filtering is done
#     """
#     output: 
#         barcodes = temp(touch("temp/{sample}/fragments/excluded_barcodes_unfiltered.tsv")),
#         qc = temp(touch("temp/{sample}/fragments/multiplet_stats_unfiltered.txt"))

# rule output_fragments_qc:
#     """
#     Move fragments qc files to final location
#     """
#     input: 
#         bc = lambda w: f"temp/{w.sample}/fragments/excluded_barcodes_{'filtered' if sample_config[w.sample]['modality'] == '10x' else 'unfiltered'}.tsv",
#         qc = lambda w: f"temp/{w.sample}/fragments/multiplet_stats_{'filtered' if sample_config[w.sample]['modality'] == '10x' else 'unfiltered'}.txt",
#     output: 
#         bc = "results/{sample}/fragments/excluded_barcodes.tsv",
#         qc = "results/{sample}/fragments/multiplet_stats.txt"
#     conda:
#         "../envs/fragments.yaml"
#     group: 
#         "fragments"
#     shell: 
#         "cp {input.bc} {output.bc}; "
#         "cp {input.qc} {output.qc}"

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
        bam = "results/{sample}/filtering/filtered.bam",
        input_data = "results/{sample}/input_data.json"
    output: 
        "results/{sample}/fragments/fragments_metadata.json"
    params:
        output_group = "fragments",
        sample_data = lambda w: sample_data[w.sample]
    conda:
        "../envs/fragments.yaml"
    group: 
        "fragments"
    script: 
        "../scripts/write_file_metadata.py"

rule metadata_qc_fragments:
    """
    Write fragments QC metadata
    """
    input: 
        data_file = "results/{sample}/fragments/fragments.tar.gz",
        multiplets_strict = "results/{sample}/fragments/multiplet_barcodes_strict.tsv",
        multiplets_expanded = "results/{sample}/fragments/multiplet_barcodes_expanded.tsv",
        multiplet_stats = "results/{sample}/fragments/multiplet_stats.txt"
    output: 
        fragments_stats = "results/{sample}/fragments/fragments_qc_metadata.json"
    params:
        output_group = "fragments",
        sample_data = lambda w: sample_data[w.sample]
    conda:
        "../envs/fragments.yaml"
    group: 
        "fragments"
    script: 
        "../scripts/write_qc_metadata.py"

rule fragments_done:
    """
    Touch flag file upon group completion
    """
    input: 
        "results/{sample}/fragments/fragments.tsv.gz",
        "results/{sample}/fragments/fragments.tsv.gz.tbi",
        "results/{sample}/fragments/fragments.tar.gz",
        "results/{sample}/fragments/multiplet_barcodes_strict.tsv",
        "results/{sample}/fragments/multiplet_barcodes_expanded.tsv",
        "results/{sample}/fragments/multiplet_stats.txt",
        "results/{sample}/fragments/fragments_metadata.json",
        "results/{sample}/fragments/fragments_qc_metadata.json"
    output:
        touch("results/{sample}/fragments/fragments.done")
    group: 
        "fragments"