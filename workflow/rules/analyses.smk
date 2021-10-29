"""
First-pass analyses
"""

rule archr_build:
    """
    ArchR analyses
    """
    input:
        frag = "results/{sample}/fragments/fragments.tsv.gz",
        frag_ind = "results/{sample}/fragments/fragments.tsv.gz.tbi"
    output:
        project_tar = "results/{sample}/analyses/archr_project.tar.gz",
        qc_dir = temp(directory("temp/{sample}/analyses/qc")),
        project_dir = temp(directory("temp/{sample}/analyses/archr_project")),
        qc_ds_pdf = temp("temp/{sample}/analyses/qc/{sample}/{sample}-Doublet-Summary.pdf"),
        qc_ds_rds = temp("temp/{sample}/analyses/qc/{sample}/{sample}-Doublet-Summary.rds"),
        qc_frag = temp("temp/{sample}/analyses/qc/{sample}/{sample}-Fragment_Size_Distribution.pdf"),
        qc_meta = temp("temp/{sample}/analyses/qc/{sample}/{sample}-Pre-Filter-Metadata.rds"),
        qc_tss = temp("temp/{sample}/analyses/qc/{sample}/{sample}-TSS_by_Unique_Frags.pdf")
    params:
        sample_name = lambda w: w.sample,
        genome = lambda w: config["genome"][sample_config[w.sample]["genome"]]["archr_genome"],
        seed = config["archr_seed"],
    log:
        console = "logs/{sample}/analyses/archr/console.log",
        arrow_create = "logs/{sample}/analyses/archr/arrow_create.log",
        doublets = "logs/{sample}/analyses/archr/doublets.log",
        lsi = "logs/{sample}/analyses/archr/lsi.log",
        cluster = "logs/{sample}/analyses/archr/cluster.log",
        marker_genes = "logs/{sample}/analyses/archr/marker_genes.log",
        pseudobulk_rep = "logs/{sample}/analyses/archr/pseudobulk_rep.log",
        peak_call = "logs/{sample}/analyses/archr/peak_call.log",
        peak_matrix = "logs/{sample}/analyses/archr/peak_matrix.log",
        marker_peaks = "logs/{sample}/analyses/archr/marker_peaks.log",
        fetch_motif = "logs/{sample}/analyses/archr/fetch_motif.log",
        enr_motif = "logs/{sample}/analyses/archr/enr_motif.log",
        fetch_tf = "logs/{sample}/analyses/archr/fetch_tf.log",
        enr_tf = "logs/{sample}/analyses/archr/enr_tf.log",
        save = "logs/{sample}/analyses/archr/save.log"
    threads:
        max_threads
    conda:
        "../envs/analyses.yaml"
    group:
        "analyses"
    script:
        "../scripts/build_archr_project.R"

rule write_archr_qc_pdf:
    """
    Move ArchR QC PDFs
    """
    input:
        qc_ds_pdf = "temp/{sample}/analyses/qc/{sample}/{sample}-Doublet-Summary.pdf",
        qc_frag = "temp/{sample}/analyses/qc/{sample}/{sample}-Fragment_Size_Distribution.pdf",
        qc_tss = "temp/{sample}/analyses/qc/{sample}/{sample}-TSS_by_Unique_Frags.pdf"
    output:
        qc_ds_pdf = "results/{sample}/analyses/archr_doublet_summary.pdf",
        qc_frag = "results/{sample}/analyses/archr_fragment_size_distribution.pdf",
        qc_tss = "results/{sample}/analyses/archr_tss_by_unique_frags.pdf"
    conda:
        "../envs/analyses.yaml"
    group:
        "analyses"
    shell:
        "cp {input.qc_ds_pdf} {output.qc_ds_pdf}; "
        "cp {input.qc_frag} {output.qc_frag}; "
        "cp {input.qc_tss} {output.qc_tss}; "

rule parse_archr_qc:
    """
    Parse ArchR QC RDSs
    """
    input:
        qc_ds_data = "temp/{sample}/analyses/qc/{sample}/{sample}-Doublet-Summary.rds",
        qc_meta = "temp/{sample}/analyses/qc/{sample}/{sample}-Pre-Filter-Metadata.rds",
    output:
        qc_ds_data = "results/{sample}/analyses/archr_doublet_summary.tsv",
        qc_meta = "results/{sample}/analyses/archr_pre_filter_metadata.tsv",
    log:
        "logs/{sample}/analyses/archr_qc_parse.log"
    conda:
        "../envs/analyses.yaml"
    group:
        "analyses"
    script:
        "../scripts/parse_archr_qc.R"


rule analyses_done:
    """
    Touch flag file upon group completion
    """
    input: 
        "results/{sample}/analyses/archr_project.tar.gz",
        "results/{sample}/analyses/archr_doublet_summary.pdf",
        "results/{sample}/analyses/archr_doublet_summary.tsv", 
        "results/{sample}/analyses/archr_fragment_size_distribution.pdf",
        "results/{sample}/analyses/archr_pre_filter_metadata.tsv",
        "results/{sample}/analyses/archr_tss_by_unique_frags.pdf",
    output:
        touch("results/{sample}/analyses/analyses.done")
    # group: 
    #     "analyses"