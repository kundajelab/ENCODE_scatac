rule all:
    """
    Generate all outputs (default)
    """
    input: 
        expand("results/{sample}/fastqs/R1_trim.fastq.gz", sample=samples.keys()),
        expand("results/{sample}/fastqs/R2_trim.fastq.gz", sample=samples.keys()),
        expand("results/{sample}/fastqs/barcode_matching.tsv", sample=samples.keys()),
        expand("results/{sample}/fastqs/trim_adapters.txt", sample=samples.keys()),
        expand("results/{sample}/fastqs/R1_trim_metadata.json", sample=samples.keys()),
        expand("results/{sample}/fastqs/R2_trim_metadata.json", sample=samples.keys()),
        expand("results/{sample}/mapping/raw.bam", sample=samples.keys()),
        expand("results/{sample}/mapping/raw.bam.bai", sample=samples.keys()),
        expand("results/{sample}/mapping/samstats_raw.txt", sample=samples.keys()),
        expand("results/{sample}/mapping/raw_bam_metadata.json", sample=samples.keys()),
        expand("results/{sample}/filtering/filtered.bam", sample=samples.keys()),
        expand("results/{sample}/filtering/filtered.bam.bai", sample=samples.keys()),
        expand("results/{sample}/filtering/frac_mito.txt", sample=samples.keys()),
        expand("results/{sample}/filtering/markdup.txt", sample=samples.keys()),
        expand("results/{sample}/filtering/pbc_stats.tsv", sample=samples.keys()),
        expand("results/{sample}/filtering/samstats_filtered.txt", sample=samples.keys()),
        expand("results/{sample}/filtering/filtered_bam_metadata.json", sample=samples.keys()),
        expand("results/{sample}/fragments/fragments.tsv.gz", sample=samples.keys()),
        expand("results/{sample}/fragments/fragments.tsv.gz.tbi", sample=samples.keys()),
        expand("results/{sample}/fragments/fragments.tar.gz", sample=samples.keys()),
        expand("results/{sample}/fragments/fragments_metadata.json", sample=samples.keys()),
        expand("results/{sample}/analyses/archr_project.tar.gz", sample=samples.keys()),
        expand("results/{sample}/analyses/archr_doublet_summary.pdf", sample=samples.keys()),
        expand("results/{sample}/analyses/archr_doublet_summary.tsv", sample=samples.keys()),
        expand("results/{sample}/analyses/archr_fragment_size_distribution.pdf", sample=samples.keys()),
        expand("results/{sample}/analyses/archr_pre_filter_metadata.tsv", sample=samples.keys()),
        expand("results/{sample}/analyses/archr_tss_by_unique_frags.pdf", sample=samples.keys())

rule submit:
    """
    Submit outputs to ENCODE portal
    """
    input: 
        expand("results/{sample}/submit/R1_trim_submit.done", sample=samples.keys()),
        expand("results/{sample}/submit/R2_trim_submit.done", sample=samples.keys()),
        expand("results/{sample}/submit/raw_bam_submit.done", sample=samples.keys()),
        expand("results/{sample}/submit/filtered_bam_submit.done", sample=samples.keys()),
        expand("results/{sample}/submit/fragments_submit.done", sample=samples.keys())


"""
# ######################
# Download ENCODE input data
# ######################
# """

checkpoint query_portal:
    """
    Query ENCODE portal by sample accession
    """
    output:
        "results/{sample}/input_data.json"
    params:
        sample = lambda w: w.sample,
        modality = lambda w: modalities[w.sample],
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    conda:
        "envs/query.yaml"
    group: 
        "query"
    script:
        "scripts/encode_query.py"

def get_input_data(sample):
    """
    Function for reading queried input data
    """
    with checkpoints.query_portal.get(sample=sample).output[0].open() as f:
        data = json.read(f)
    return data

sample_data = defaultdict(get_input_data) # Structure for caching input data


"""
######################
FASTQ processing
######################
"""

rule strip_fastq:
    """
    Strip FASTQ read descriptions
    """
    input:
        get_input_fastq 
    output:
        pipe("temp/{sample}/fastqs/stripped_{read}.fastq")
    conda:
        "envs/fastqs.yaml"
    group: 
        "fastqs"
    shell:
        "zcat {input} | sed 's/ .*//' > {output}"

rule match_barcodes:  
    """
    Barcode correction and filtering
    """
    # base rule across technologies
    threads:
        max_threads
    conda:
        "envs/fastqs.yaml"
    group: 
        "fastqs"
    script:
        "scripts/match_barcodes.py"

use rule match_barcodes as match_barcodes_10x with:
    # config["technology"] == "10x"
    input: 
        fq_R1 = "temp/{sample}/fastqs/stripped_R1.fastq",
        fq_R2 = "temp/{sample}/fastqs/stripped_R2.fastq",
        fq_R3 = "temp/{sample}/fastqs/stripped_R3.fastq",
        wl_R2 = lambda w: samples[w.sample]["bc_whitelist"]["R2"]
    output: 
        fastq1_bc = temp("temp/{sample}/fastqs/R1_bc_10x.fastq"),
        fastq2_bc = temp("temp/{sample}/fastqs/R2_bc_10x.fastq"),
        qc_matching = temp("temp/{sample}/fastqs/barcode_matching_10x.tsv")
    params:
        rc_R2 = lambda w: samples[w.sample]["bc_revcomp"]["R2"],
        barcode_dist = lambda w: samples[w.sample]["max_barcode_dist"],
        technology = "10x"

use rule match_barcodes as match_barcodes_multiome with:
    # config["technology"] == "multiome"
    input: 
        fq_R1 = "temp/{sample}/fastqs/stripped_R1.fastq",
        fq_R2 = "temp/{sample}/fastqs/stripped_R2.fastq",
        fq_R3 = "temp/{sample}/fastqs/stripped_R3.fastq",
        wl_R2 = lambda w: samples[w.sample]["bc_whitelist"]["R2"]
    output: 
        fastq1_bc = temp("temp/{sample}/fastqs/R1_bc_multiome.fastq"),
        fastq2_bc = temp("temp/{sample}/fastqs/R2_bc_multiome.fastq"),
        qc_matching = temp("temp/{sample}/fastqs/barcode_matching_multiome.tsv")
    params:
        rc_R2 = lambda w: samples[w.sample]["bc_revcomp"]["R2"],
        barcode_dist = lambda w: samples[w.sample]["max_barcode_dist"],
        technology = "multiome"

use rule match_barcodes as match_barcodes_ren with:
    # config["technology"] == "ren"
    input:
        fq_R1 = "temp/{sample}/fastqs/stripped_R1.fastq",
        fq_R2 = "temp/{sample}/fastqs/stripped_R2.fastq",
        fq_I1 = "temp/{sample}/fastqs/stripped_I1.fastq",
        fq_I2 = "temp/{sample}/fastqs/stripped_I2.fastq",
        wl_I1 = lambda w: samples[w.sample]["bc_whitelist"]["I1"],
        wl_I2 = lambda w: samples[w.sample]["bc_whitelist"]["I2"]
    output: 
        fastq1_bc = temp("temp/{sample}/fastqs/R1_bc_ren.fastq"),
        fastq2_bc = temp("temp/{sample}/fastqs/R2_bc_ren.fastq"),
        qc_matching = temp("temp/{sample}/fastqs/barcode_matching_ren.tsv")
    params:
        rc_I1 = lambda w: samples[w.sample]["bc_revcomp"]["I1"],
        rc_I2 = lambda w: samples[w.sample]["bc_revcomp"]["I2"],
        barcode_dist = lambda w: samples[w.sample]["max_barcode_dist"],
        technology = "ren"

rule output_matching_stats:
    """
    Write barcode matching stats
    """
    input:
        lambda w: f"temp/{w.sample}/fastqs/barcode_matching_{samples[w.sample]['technology']}.tsv",
    output:
        "results/{sample}/fastqs/barcode_matching.tsv"
    conda:
        "envs/fastqs.yaml"
    group: 
        "fastqs"
    shell:
        "cp {input} {output}"

rule trim_adapter:
    """
    Read adapter trimming
    """
    input:
        fastq1_bc = lambda w: f"temp/{w.sample}/fastqs/R1_bc_{samples[w.sample]['technology']}.fastq",
        fastq2_bc = lambda w: f"temp/{w.sample}/fastqs/R2_bc_{samples[w.sample]['technology']}.fastq"
    output:
        fastq1_trim = "results/{sample}/fastqs/R1_trim.fastq.gz",
        fastq2_trim = "results/{sample}/fastqs/R2_trim.fastq.gz",
        stats = "results/{sample}/fastqs/trim_adapters.txt"
    log:
        html = "logs/{sample}/fastqs/fastp.html",
        json = "logs/{sample}/fastqs/fastp.json"
    threads:
        max_threads
    conda:
        "envs/fastqs.yaml"
    group: 
        "fastqs"
    shell:
        "fastp -i {input.fastq1_bc} -I {input.fastq2_bc} -o {output.fastq1_trim} -O {output.fastq2_trim}"
        " -h {log.html} -j {log.json} -G -Q -L -w {threads} 2> {output.stats}"

rule metadata_fastq:
    """
    Write FASTQ metadata
    """
    input: 
        r1 = "results/{sample}/fastqs/R1_trim.fastq.gz",
        r2 = "results/{sample}/fastqs/R2_trim.fastq.gz"
    output: 
        r1 = "results/{sample}/fastqs/R1_trim_metadata.json",
        r2 = "results/{sample}/fastqs/R2_trim_metadata.json"
    params:
        output_group = "fastqs",
        sample_data = lambda w: samples[w.sample]
    conda:
        "envs/fastqs.yaml"
    group: 
        "fastqs"
    script: 
        "scripts/write_file_metadata.py"


"""
######################
Read mapping
######################
"""

rule bowtie2:
    """
    Read mapping (Bowtie2 aligner)
    """
    input:
        fastq1 = "results/{sample}/fastqs/R1_trim.fastq.gz",
        fastq2 = "results/{sample}/fastqs/R2_trim.fastq.gz"
    output:
        bam_raw = pipe("temp/{sample}/mapping/raw.bam"),
    params:
        k = 1 + config["multimapping"],
        bwt2_idx = config["mapping_index"]
    log:
        "logs/{sample}/mapping/bwt2.log"
    threads:
        max_threads
    conda:
        "envs/mapping.yaml"
    group: 
        "mapping"
    shell:
        "bowtie2 -X 2000 --threads {threads} -x {params.bwt2_idx} "
        "-1 {input.fastq1} -2 {input.fastq2} --sam-append-comment -k {params.k} 2> {log} | "
        "samtools view -u -S -o {output.bam_raw} -"

rule filter_multimappers:
    """
    Remove multimapping reads above threshold
    """
    input:
        "temp/{sample}/mapping/raw.bam"
    output:
        pipe("temp/{sample}/mapping/de-multimap.bam")
    params:
        multimapping = config["multimapping"],
        mmp_path = script_path("scripts/assign_multimappers.py")
    conda:
        "envs/mapping.yaml"
    group: 
        "mapping"
    shell:
        "samtools view -h -f 2 {input} | "
        "python {params.mmp_path} --paired-end -k {params.multimapping} | "
        "samtools view -u -o {output} - "

rule sort_alignments:
    """
    Sort and index alignments
    """
    input: 
        "temp/{sample}/mapping/de-multimap.bam"
    output: 
        "results/{sample}/mapping/raw.bam"
    log:
        "logs/{sample}/mapping/sort.log"
    threads:
        max_threads
    conda:
        "envs/mapping.yaml"
    shadow: 
        "minimal"
    group: 
        "mapping"
    shell: 
        "samtools sort -T . -@ {threads} -o {output} {input} 2> {log}; "
        "samtools index {output.bam};"

rule index_bam_raw:
    """
    Index raw BAM
    """
    input: 
        "results/{sample}/mapping/raw.bam"
    output: 
        "results/{sample}/mapping/raw.bam.bai"
    conda:
        "envs/mapping.yaml"
    group: 
        "mapping"
    shell: 
        "samtools index {input}"

rule samstats_raw:
    """
    Run SAMStats on raw alignments
    """
    input:
        "results/{sample}/mapping/raw.bam"
    output:
        "results/{sample}/mapping/samstats_raw.txt"
    log:
        "logs/{sample}/mapping/samstats_raw.log"
    threads:
        max_threads
    conda:
        "envs/filtering.yaml"
    group: 
        "mapping"
    shadow: 
        "minimal"
    shell:
        "samtools sort -T . -n -@ {threads} -O SAM {input} | " 
        "SAMstats --sorted_sam_file -  --outf {output} > {log}"

rule metadata_bam_raw:
    """
    Write raw BAM metadata
    """
    input: 
        bam = "results/{sample}/mapping/raw.bam",
        fq_R1 = "results/{sample}/fastqs/R1_trim.fastq.gz",
        fq_R2 = "results/{sample}/fastqs/R2_trim.fastq.gz"
    output: 
        "results/{sample}/mapping/raw_bam_metadata.json"
    params:
        output_group = "mapping",
        sample_data = lambda w: samples[w.sample]
    conda:
        "envs/mapping.yaml"
    group: 
        "mapping"
    script: 
        "scripts/write_file_metadata.py"

"""
######################
Alignment filtering
######################
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
        "envs/filtering.yaml"
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
        bam = pipe("temp/{sample}/filtering/no_mito.bam"),
        count_no_mito = temp("temp/{sample}/filtering/count_no_mito.txt")
    conda:
        "envs/filtering.yaml"
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
        "results/{sample}/filtering/frac_mito.txt"
    conda:
        "envs/filtering.yaml"
    group: 
        "filtering"
    shell: 
        "rm=$(<{input.count_mito}); "
        "rn=$(<{input.count_no_mito}); "
        "frac=$(($rm / ($rm + $rn))); "
        "printf \"%d\\t%d\\t%f\" \"$rn\" \"$rm\" \"$frac\" > {output}"

rule assign_primary:
    """
    Assign multimapping reads to primary alignment
    """
    input:
        "temp/{sample}/filtering/no_mito.bam"
    output:
        temp("temp/{sample}/filtering/primary_align.bam")
    conda:
        "envs/filtering.yaml"
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
        bam_markdup = pipe("temp/{sample}/filtering/markdup.bam"),
        bam_nodup = "results/{sample}/filtering/filtered.bam",
        markdup_stats = "results/{sample}/filtering/markdup.txt"
    log:
        "logs/{sample}/filtering/picard.log"
    conda:
        "envs/filtering.yaml"
    group: 
        "filtering"
    shell:
        "picard MarkDuplicates --INPUT {input} --OUTPUT /dev/stdout --METRICS_FILE {output.markdup_stats} --COMPRESSION_LEVEL 0 "
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
        pbc_script = srcdir("scripts/pbc_stats.py")
    conda:
        "envs/filtering.yaml"
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
        "envs/filtering.yaml"
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
        "envs/filtering.yaml"
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
        fq_R1 = "results/{sample}/fastqs/R1_trim.fastq.gz",
        fq_R2 = "results/{sample}/fastqs/R2_trim.fastq.gz"
    output: 
        "results/{sample}/filtering/filtered_bam_metadata.json",
    params:
        output_group = "filtering",
        sample_data = lambda w: samples[w.sample]
    conda:
        "envs/filtering.yaml"
    group: 
        "filtering"
    script: 
        "scripts/write_file_metadata.py"


"""
######################
Fragment file generation
######################
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
        "envs/fragments.yaml"
    group: 
        "fragments"
    shell:
        "sinto fragments -b {input.bam} -f {output} " 
        "--min_mapq 0 --max_distance 2000 --min_distance 10 --barcodetag CB --nproc {threads} > {log}"

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
        "envs/fragments.yaml"
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
        "envs/fragments.yaml"
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
        "envs/fragments.yaml"
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
        "envs/fragments.yaml"
    group: 
        "fragments"
    script: 
        "scripts/write_file_metadata.py"


"""
######################
First-pass analyses
######################
"""

rule archr_build:
    """
    ArchR analyses
    """
    input:
        frag = "results/{sample}/fragments/fragments.tsv.gz",
        frag_ind = "results/{sample}/fragments/fragments.tsv.gz.tbi"
    output:
        qc_dir =temp(directory("temp/{sample}/analyses/qc")),
        project_dir = temp(directory("temp/{sample}/analyses/archr_project")),
        flag = touch("temp/{sample}/analyses/archr_flag.txt"),
        qc_ds_pdf = temp("temp/{sample}/analyses/qc/{sample}/{sample}-Doublet-Summary.pdf"),
        qc_ds_rds = temp("temp/{sample}/analyses/qc/{sample}/{sample}-Doublet-Summary.rds"),
        qc_frag = temp("temp/{sample}/analyses/qc/{sample}/{sample}-Fragment_Size_Distribution.pdf"),
        qc_meta = temp("temp/{sample}/analyses/qc/{sample}/{sample}-Pre-Filter-Metadata.rds"),
        qc_tss = temp("temp/{sample}/analyses/qc/{sample}/{sample}-TSS_by_Unique_Frags.pdf")
    params:
        sample_name = lambda w: w.sample,
        seed = config["seed"],
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
        "envs/analyses.yaml"
    group:
        "analyses"
    shadow: 
        "shallow"
    script:
        "scripts/build_archr_project.R"

rule tar_archr_results:
    """
    Create ArchR results archive
    """
    input:
        "temp/{sample}/analyses/archr_flag.txt"
    params:
        project_dir = lambda w: f"temp/{w.sample}/analyses/archr_project",
    output:
        "results/{sample}/analyses/archr_project.tar.gz"
    conda:
        "envs/analyses.yaml"
    group:
        "analyses"
    shell:
        "tar -zcf {output} {params.project_dir}"

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
        "envs/analyses.yaml"
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
        "envs/analyses.yaml"
    group:
        "analyses"
    script:
        "scripts/parse_archr_qc.R"


"""
######################
ENCODE portal submission
######################
Note: this group is not run by default
"""

rule submit_fastq_1:
    """
    Submit FASTQ pair 1
    """
    input: 
        json = "results/{sample}/fastqs/R1_trim_metadata.json"
    output: 
        touch("results/{sample}/submit/R1_trim_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/R1_trim_submit")
    conda:
        "envs/submit.yaml"
    group: 
        "submit"
    script: 
        "scripts/encode_submit.py"

rule submit_fastq_2:
    """
    Submit FASTQ pair 2
    """
    input: 
        json = "results/{sample}/fastqs/R2_trim_metadata.json",
        fastq1 = "results/{sample}/submit/R1_trim_submit.done"
    output: 
        touch("results/{sample}/submit/R2_trim_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/R2_trim_submit")
    conda:
        "envs/submit.yaml"
    group: 
        "submit"
    script: 
        "scripts/encode_submit.py"

rule submit_bam_raw:
    """
    Submit raw BAM
    """
    input: 
        json = "results/{sample}/mapping/raw_bam_metadata.json",
        fastq1 = "results/{sample}/submit/R1_trim_submit.done",
        fastq2 = "results/{sample}/submit/R2_trim_submit.done"
    output: 
        touch("results/{sample}/submit/raw_bam_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/raw_bam_submit")
    conda:
        "envs/submit.yaml"
    group: 
        "submit"
    script: 
        "scripts/encode_submit.py"

rule submit_bam_filtered:
    """
    Submit filtered BAM
    """
    input: 
        json = "results/{sample}/filtering/filtered_bam_metadata.json",
        fastq1 = "results/{sample}/submit/R1_trim_submit.done",
        fastq2 = "results/{sample}/submit/R2_trim_submit.done"
    output: 
        touch("results/{sample}/submit/filtered_bam_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/filtered_bam_submit")
    conda:
        "envs/submit.yaml"
    group: 
        "submit"
    script: 
        "scripts/encode_submit.py"

rule submit_fragments:
    """
    Submit fragment file
    """
    input: 
        json = "results/{sample}/fragments/fragments_metadata.json",
        bam = "results/{sample}/submit/filtered_bam_submit.done"
    output: 
        touch("results/{sample}/submit/fragments_submit.done")
    params:
        schema = "file",
        dcc_api_key = os.environ["DCC_API_KEY"], 
        dcc_secret_key = os.environ["DCC_SECRET_KEY"]
    log:
        directory("logs/{sample}/submit/filtered_bam_submit")
    conda:
        "envs/submit.yaml"
    group: 
        "submit"
    script: 
        "scripts/encode_submit.py"