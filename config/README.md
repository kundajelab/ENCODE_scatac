# Pipeline Configuration

🚧 🚧 🚧 Under construction 🚧 🚧 🚧

## Per-Sample Configuration

Per-sample configuration is stored in `samples.tsv`, where each row corresponds to one run of the pipeline. 

**Experiment** is the ENCODE accession ID of the sample's dataset.

**Replicate** is the biological replicate number of the sample on the ENCODE portal.

**Modality** specifies the experimental protocol used to generate the sample's raw data. Currently supported modalities are:
- 10x single-cell ATAC (**10x**)
- 10x Single Cell Multiome ATAC + Gene Expression (**multiome**)
- Bing Ren lab split-pool scATAC (**ren**)

**Genome** specifies the target genome for alignment and analysis. The pipeline currently supports **GRCh38** (human) and **mm10** (mouse).