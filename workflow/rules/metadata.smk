rule metadata_fastq:
    """
    Write FASTQ metadata
    """
    input: 
        r1 = "results/{sample}/fastqs/R1_trim.fastq.gz",
        r2 = "results/{sample}/fastqs/R2_trim.fastq.gz",
        input_data = "results/{sample}/input_data.json"
    output: 
        r1 = "metadata/{sample}/R1_trim_metadata.json",
        r2 = "metadata/{sample}/R2_trim_metadata.json"
    params:
        output_group = "fastqs",
        sample_data = lambda w: sample_data[w.sample]
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata"
    script: 
        "../scripts/write_file_metadata.py"

rule collate_R1_trim:
    """
    Collate metadata files across samples into tables
    """
    input: 
        expand("metadata/{sample}/R1_trim_metadata.json", sample=samples)
    output:
        "metadata/R1_trim_metadata_all.tsv"
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata_collate"
    script: 
        "../scripts/collate_metadata.py"

rule collate_R2_trim:
    """
    Collate metadata files across samples into tables
    """
    input: 
        expand("metadata/{sample}/R2_trim_metadata.json", sample=samples)
    output:
        "metadata/R2_trim_metadata_all.tsv"
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata_collate"
    script: 
        "../scripts/collate_metadata.py"

rule metadata_qc_reads:
    """
    Write reads QC metadata
    """
    input: 
        data_file = "results/{sample}/fastqs/R1_trim.fastq.gz",
        barcode_matching = "results/{sample}/fastqs/barcode_matching.tsv",
        adapter_trimming = "results/{sample}/fastqs/trim_adapters.txt",
        barcode_revcomp = "results/{sample}/fastqs/barcode_revcomp.txt",
        input_data = "results/{sample}/input_data.json"
    output: 
        read_stats = "metadata/{sample}/reads_qc_metadata.json",
    params:
        output_group = "fastqs",
        sample_data = lambda w: sample_data[w.sample]
     conda:
        "../envs/portal.yaml"
    group: 
        "metadata"
    script: 
        "../scripts/write_qc_metadata.py"

rule collate_qc_reads:
    """
    Collate metadata files across samples into tables
    """
    input: 
        expand("metadata/{sample}/reads_qc_metadata.json", sample=samples)
    output:
        "metadata/reads_qc_metadata_all.tsv"
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata_collate"
    script: 
        "../scripts/collate_metadata.py"

rule metadata_bam_raw:
    """
    Write raw BAM metadata
    """
    input: 
        bam = "results/{sample}/mapping/raw.bam",
        fq_R1 = "results/{sample}/fastqs/R1_trim.fastq.gz",
        fq_R2 = "results/{sample}/fastqs/R2_trim.fastq.gz",
        input_data = "results/{sample}/input_data.json"
    output: 
        "metadata/{sample}/raw_bam_metadata.json"
    params:
        output_group = "mapping",
        sample_data = lambda w: sample_data[w.sample]
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata"
    script: 
        "../scripts/write_file_metadata.py"

rule collate_bam_raw:
    """
    Collate metadata files across samples into tables
    """
    input: 
        expand("metadata/{sample}/raw_bam_metadata.json", sample=samples)
    output:
        "metadata/raw_bam_metadata_all.tsv"
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata_collate"
    script: 
        "../scripts/collate_metadata.py"

rule metadata_qc_alignments_raw:
    """
    Write raw alignments qc metadata
    """
    input: 
        data_file = "results/{sample}/mapping/raw.bam",
        samstats_raw = "results/{sample}/mapping/samstats_raw.txt",
        input_data = "results/{sample}/input_data.json"
    output: 
        alignment_stats = "metadata/{sample}/alignments_raw_qc_metadata.json"
    params:
        output_group = "mapping",
        sample_data = lambda w: sample_data[w.sample]
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata"
    script: 
        "../scripts/write_qc_metadata.py"

rule collate_qc_alignments_raw:
    """
    Collate metadata files across samples into tables
    """
    input: 
        expand("metadata/{sample}/alignments_raw_qc_metadata.json", sample=samples)
    output:
        "metadata/alignments_raw_qc_metadata_all.tsv"
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata_collate"
    script: 
        "../scripts/collate_metadata.py"

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
        "metadata/{sample}/filtered_bam_metadata.json",
    params:
        output_group = "filtering",
        sample_data = lambda w: sample_data[w.sample]
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata"
    script: 
        "../scripts/write_file_metadata.py"

rule collate_bam_filtered:
    """
    Collate metadata files across samples into tables
    """
    input: 
        expand("metadata/{sample}/filtered_bam_metadata.json", sample=samples)
    output:
        "metadata/filtered_bam_metadata_all.tsv"
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata_collate"
    script: 
        "../scripts/collate_metadata.py"

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
        alignment_stats = "metadata/{sample}/alignments_filtered_qc_metadata.json",
        lib_comp_stats = "metadata/{sample}/alignments_lib_comp_qc_metadata.json"
    params:
        output_group = "filtering",
        sample_data = lambda w: sample_data[w.sample]
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata"
    script: 
        "../scripts/write_qc_metadata.py"

rule collate_qc_alignments_filtered:
    """
    Collate metadata files across samples into tables
    """
    input: 
        expand("metadata/{sample}/alignments_filtered_qc_metadata.json", sample=samples)
    output:
        "metadata/alignments_filtered_qc_metadata_all.tsv"
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata_collate"
    script: 
        "../scripts/collate_metadata.py"

rule collate_qc_lib_comp:
    """
    Collate metadata files across samples into tables
    """
    input: 
        expand("metadata/{sample}/alignments_lib_comp_qc_metadata.json", sample=samples)
    output:
        "metadata/alignments_lib_comp_qc_metadata_all.tsv"
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata_collate"
    script: 
        "../scripts/collate_metadata.py"

rule metadata_fragments:
    """
    Write fragment file metadata
    """
    input: 
        fragments = "results/{sample}/fragments/fragments.tar.gz",
        bam = "results/{sample}/filtering/filtered.bam",
        input_data = "results/{sample}/input_data.json"
    output: 
        "metadata/{sample}/fragments_metadata.json"
    params:
        output_group = "fragments",
        sample_data = lambda w: sample_data[w.sample]
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata"
    script: 
        "../scripts/write_file_metadata.py"

rule collate_fragments:
    """
    Collate metadata files across samples into tables
    """
    input: 
        expand("metadata/{sample}/fragments_metadata.json", sample=samples)
    output:
        "metadata/fragments_metadata_all.tsv"
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata_collate"
    script: 
        "../scripts/collate_metadata.py"

rule metadata_qc_fragments:
    """
    Write fragments QC metadata
    """
    input: 
        data_file = "results/{sample}/fragments/fragments.tar.gz",
        barcode_pairs_strict = "results/{sample}/fragments/barcode_pairs_multiplets.tsv",
        barcode_pairs_expanded = "results/{sample}/fragments/barcode_pairs_expanded.tsv.gz",
        barcodes_status = "results/{sample}/fragments/multiplet_barcodes_status.tsv",
        multiplet_stats = "results/{sample}/fragments/multiplet_stats.txt",
        multiplets_thresh = "results/{sample}/fragments/multiplets_threshold_plot.png"
    output: 
        fragments_stats = "metadata/{sample}/fragments_qc_metadata.json"
    params:
        output_group = "fragments",
        sample_data = lambda w: sample_data[w.sample]
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata"
    script: 
        "../scripts/write_qc_metadata.py"

rule collate_qc_fragments:
    """
    Collate metadata files across samples into tables
    """
    input: 
        expand("metadata/{sample}/fragments_qc_metadata.json", sample=samples)
    output:
        "metadata/fragments_qc_metadata_all.tsv"
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata_collate"
    script: 
        "../scripts/collate_metadata.py"

rule metadata_analyses:
    """
    Write analyses metadata
    """
    input: 
        archr = "results/{sample}/analyses/archr_project.tar.gz",
        fragments = "results/{sample}/fragments/fragments.tar.gz",
        input_data = "results/{sample}/input_data.json"
    output: 
        "metadata/{sample}/analyses_metadata.json",
    params:
        output_group = "analyses",
        sample_data = lambda w: sample_data[w.sample]
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata"
    script: 
        "../scripts/write_file_metadata.py"

rule collate_analyses:
    """
    Collate metadata files across samples into tables
    """
    input: 
        expand("metadata/{sample}/analyses_metadata.json", sample=samples)
    output:
        "metadata/analyses_metadata_all.tsv"
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata_collate"
    script: 
        "../scripts/collate_metadata.py"

rule metadata_qc_analyses:
    """
    Write filtered alignments qc metadata
    """
    input: 
        data_file = "results/{sample}/fragments/fragments.tar.gz", # attach to fragments
        archr_doublet_summary_text = "results/{sample}/analyses/archr_doublet_summary.tsv",
        archr_doublet_summary_figure =  "results/{sample}/analyses/archr_doublet_summary.pdf",
        archr_fragment_size_distribution = "results/{sample}/analyses/archr_fragment_size_distribution.pdf",
        archr_pre_filter_metadata = "results/{sample}/analyses/archr_pre_filter_metadata.tsv",
        archr_tss_by_unique_frags = "results/{sample}/analyses/archr_tss_by_unique_frags.pdf",
        input_data = "results/{sample}/input_data.json"
    output: 
        analyses_stats = "metadata/{sample}/analyses_qc_metadata.json",
    params:
        output_group = "analyses",
        sample_data = lambda w: sample_data[w.sample]
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata"
    script: 
        "../scripts/write_qc_metadata.py"

rule collate_qc_analyses:
    """
    Collate metadata files across samples into tables
    """
    input: 
        expand("metadata/{sample}/analyses_qc_metadata.json", sample=samples)
    output:
        "metadata/analyses_qc_metadata_all.tsv"
    conda:
        "../envs/portal.yaml"
    group: 
        "metadata_collate"
    script: 
        "../scripts/collate_metadata.py"

