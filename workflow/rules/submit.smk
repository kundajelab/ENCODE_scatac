
"""
ENCODE portal submission
"""
# Note: these rules are not run by default

rule submit_fastq_1:
    """
    Submit FASTQ pair 1
    """
    input: 
        json = "metadata/{sample}/R1_trim_metadata.json",
        r1 = "results/{sample}/fastqs/R1_trim.fastq.gz"
    output: 
        touch("submit/{sample}/R1_trim_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/R1_trim_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_fastq_2:
    """
    Submit FASTQ pair 2
    """
    input: 
        json = "metadata/{sample}/R2_trim_metadata.json",
        prev = "submit/{sample}/R1_trim_submit.done",
        r2 = "results/{sample}/fastqs/R2_trim.fastq.gz"
    output: 
        touch("submit/{sample}/R2_trim_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/R2_trim_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_reads_qc: 
    """
    Submit reads QC
    """
    input: 
        json = "metadata/{sample}/reads_qc_metadata.json",
        prev = "submit/{sample}/R1_trim_submit.done",
        data_file = "results/{sample}/fastqs/R1_trim.fastq.gz",
        barcode_matching = "results/{sample}/fastqs/barcode_matching.tsv",
        adapter_trimming = "results/{sample}/fastqs/trim_adapters.txt",
        barcode_revcomp = "results/{sample}/fastqs/barcode_revcomp.txt"
    output: 
        touch("submit/{sample}/reads_qc_metadata_submit.done")
    params:
        schema = "sc_atac_read_quality_metric",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/reads_qc_metadata_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_bam_raw:
    """
    Submit raw BAM
    """
    input: 
        json = "metadata/{sample}/raw_bam_metadata.json",
        prev1 = "submit/{sample}/R1_trim_submit.done",
        prev2 = "submit/{sample}/R2_trim_submit.done",
        bam = "results/{sample}/mapping/raw.bam",
        fq_R1 = "results/{sample}/fastqs/R1_trim.fastq.gz",
        fq_R2 = "results/{sample}/fastqs/R2_trim.fastq.gz"
    output: 
        touch("submit/{sample}/raw_bam_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/raw_bam_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_alignments_raw_qc: 
    """
    Submit raw alignments QC
    """
    input: 
        json = "metadata/{sample}/alignments_raw_qc_metadata.json",
        prev = "submit/{sample}/raw_bam_submit.done",
        data_file = "results/{sample}/mapping/raw.bam",
        samstats_raw = "results/{sample}/mapping/samstats_raw.txt"
    output: 
        touch("submit/{sample}/alignments_raw_qc_metadata_submit.done")
    params:
        schema = "sc_atac_alignment_quality_metric",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/alignments_raw_qc_metadata_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_bam_filtered:
    """
    Submit filtered BAM
    """
    input: 
        json = "metadata/{sample}/filtered_bam_metadata.json",
        prev1 = "submit/{sample}/R1_trim_submit.done",
        prev2 = "submit/{sample}/R2_trim_submit.done",
        bam = "results/{sample}/filtering/filtered.bam",
        fq_R1 = "results/{sample}/fastqs/R1_trim.fastq.gz", # attach to fastqs
        fq_R2 = "results/{sample}/fastqs/R2_trim.fastq.gz",
    output: 
        touch("submit/{sample}/filtered_bam_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/filtered_bam_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_alignments_filtered_qc: 
    """
    Submit filtered alignments QC
    """
    input: 
        json = "metadata/{sample}/alignments_filtered_qc_metadata.json",
        prev = "submit/{sample}/filtered_bam_submit.done",
        data_file = "results/{sample}/filtering/filtered.bam",
        samstats_filtered = "results/{sample}/filtering/samstats_filtered.txt",
        picard_markdup = "results/{sample}/filtering/markdup.txt",
        pbc_stats = "results/{sample}/filtering/pbc_stats.tsv",
        frac_mito = "results/{sample}/filtering/frac_mito.tsv"
    output: 
        touch("submit/{sample}/alignments_filtered_qc_metadata_submit.done")
    params:
        schema = "sc_atac_alignment_quality_metric",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/alignments_filtered_qc_metadata_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_lib_comp_qc: 
    """
    Submit library complexity QC
    """
    input: 
        json = "metadata/{sample}/alignments_lib_comp_qc_metadata.json",
        prev = "submit/{sample}/filtered_bam_submit.done",
        data_file = "results/{sample}/filtering/filtered.bam",
        samstats_filtered = "results/{sample}/filtering/samstats_filtered.txt",
        picard_markdup = "results/{sample}/filtering/markdup.txt",
        pbc_stats = "results/{sample}/filtering/pbc_stats.tsv",
        frac_mito = "results/{sample}/filtering/frac_mito.tsv",
    output: 
        touch("submit/{sample}/alignments_lib_comp_qc_metadata_submit.done")
    params:
        schema = "sc_atac_library_complexity_quality_metric", 
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/alignments_lib_comp_qc_metadata_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_fragments:
    """
    Submit fragment file
    """
    input: 
        json = "metadata/{sample}/fragments_metadata.json",
        prev = "submit/{sample}/filtered_bam_submit.done",
        fragments = "results/{sample}/fragments/fragments.tar.gz",
        bam = "results/{sample}/filtering/filtered.bam",
    output: 
        touch("submit/{sample}/fragments_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/filtered_bam_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_fragments_qc: 
    """
    Submit fragment file QC
    """
    input: 
        json = "metadata/{sample}/fragments_qc_metadata.json",
        prev = "submit/{sample}/fragments_submit.done",
        data_file = "results/{sample}/fragments/fragments.tar.gz",
        barcode_pairs_strict = "results/{sample}/fragments/barcode_pairs_multiplets.tsv",
        barcode_pairs_expanded = "results/{sample}/fragments/barcode_pairs_expanded.tsv.gz",
        barcodes_status = "results/{sample}/fragments/multiplet_barcodes_status.tsv",
        multiplet_stats = "results/{sample}/fragments/multiplet_stats.txt",
        multiplets_thresh = "results/{sample}/fragments/multiplets_threshold_plot.png"
    output: 
        touch("submit/{sample}/fragments_qc_metadata_submit.done")
    params:
        schema = "sc_atac_multiplet_quality_metric", 
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/fragments_qc_metadata_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_analyses:
    """
    Submit analyses
    """
    input: 
        json = "metadata/{sample}/analyses_metadata.json",
        prev = "submit/{sample}/fragments_submit.done",
        archr = "results/{sample}/analyses/archr_project.tar.gz",
        fragments = "results/{sample}/fragments/fragments.tar.gz",
    output: 
        touch("submit/{sample}/analyses_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/analyses_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_analyses_qc: 
    """
    Submit analyses QC
    """
    input: 
        json = "metadata/{sample}/analyses_qc_metadata.json",
        prev = "submit/{sample}/analyses_submit.done",
        data_file = "results/{sample}/fragments/fragments.tar.gz", # attach to fragments
        archr_doublet_summary_text = "results/{sample}/analyses/archr_doublet_summary.tsv",
        archr_doublet_summary_figure =  "results/{sample}/analyses/archr_doublet_summary.pdf",
        archr_fragment_size_distribution = "results/{sample}/analyses/archr_fragment_size_distribution.pdf",
        archr_pre_filter_metadata = "results/{sample}/analyses/archr_pre_filter_metadata.tsv",
        archr_tss_by_unique_frags = "results/{sample}/analyses/archr_tss_by_unique_frags.pdf",
    output: 
        touch("submit/{sample}/analyses_qc_metadata_submit.done")
    params:
        schema = "sc_atac_analysis_quality_metric", #TODO
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/analyses_qc_metadata_submit")
    conda:
        "../envs/portal.yaml"
    group: 
        "submit"
    script: 
        "../scripts/encode_submit.py"

rule submit_done:
    """
    Touch flag file upon group completion
    """
    input: 
        "submit/{sample}/R1_trim_submit.done",
        "submit/{sample}/R2_trim_submit.done",
        "submit/{sample}/reads_qc_metadata_submit.done",
        "submit/{sample}/raw_bam_submit.done",
        "submit/{sample}/filtered_bam_submit.done",
        "submit/{sample}/fragments_submit.done",
        "submit/{sample}/fragments_qc_metadata_submit.done",
        "submit/{sample}/analyses_submit.done",
        "submit/{sample}/analyses_qc_metadata_submit.done",
    output:
        touch("submit/{sample}/submit.done")
    group: 
        "submit"
