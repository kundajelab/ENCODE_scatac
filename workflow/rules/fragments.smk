"""
Fragment file generation
"""

rule bam_to_fragments: 
    """
    Convert BAM to fragment file
    """
    input:
        "results/{sample}/filtering/filtered.bam",
    output:
        temp("temp/{sample}/fragments/fragments_raw.tsv")
    params:
        shift_plus = config["tn5_shift_plus"],
        shift_minus = config["tn5_shift_minus"]
    conda:
        "../envs/fragments.yaml"
    group: 
        "fragments"
    script:
        "../scripts/bam_to_fragments.py"

rule compress_fragments:
    """
    Compress fragment file
    """
    input:
        "temp/{sample}/fragments/fragments_raw.tsv"
    output: 
        "results/{sample}/fragments/fragments.tsv.gz"
    conda:
        "../envs/fragments.yaml"
    group: 
        "fragments"
    shell: 
        "bgzip -c {input} > {output}"

rule index_fragments:
    """
    Index fragment file
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
        barcode_pairs_strict = "results/{sample}/fragments/barcode_pairs_multiplets.tsv",
        barcode_pairs_expanded = "results/{sample}/fragments/barcode_pairs_expanded.tsv.gz",
        barcodes_status = "results/{sample}/fragments/multiplet_barcodes_status.tsv",
        qc = "results/{sample}/fragments/multiplet_stats.txt",
        multiplets_thresh = "results/{sample}/fragments/multiplets_threshold_plot.png"
    conda:
        "../envs/fragments.yaml"
    group: 
        "fragments"
    script: 
        "../scripts/detect_multiplets.py"

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

# rule fragments_done:
#     """
#     Touch flag file upon group completion
#     """
#     input: 
#         "results/{sample}/fragments/fragments.tsv.gz",
#         "results/{sample}/fragments/fragments.tsv.gz.tbi",
#         "results/{sample}/fragments/fragments.tar.gz",
#         "results/{sample}/fragments/barcode_pairs_multiplets.tsv",
#         "results/{sample}/fragments/barcode_pairs_expanded.tsv.gz",
#         "results/{sample}/fragments/multiplet_barcodes_status.tsv",
#         "results/{sample}/fragments/multiplets_threshold_plot.png",
#         "results/{sample}/fragments/multiplet_stats.txt"
#     output:
#         touch("results/{sample}/fragments/fragments.done")
#     # group: 
#     #     "fragments"