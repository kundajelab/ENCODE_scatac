# Pipeline Configuration

## Per-Sample Configuration

Per-sample configuration is stored in `samples.tsv`, where each row is a single sample. Here, a sample is defined as a single biological replicate in a scATAC datset. The pipeline will pull all necessary data files from the ENCODE portal.

**Experiment** is the ENCODE accession ID of the sample's dataset.

**Replicate** is the biological replicate number of the sample on the ENCODE portal.

**Modality** specifies the experimental protocol used to generate the sample's raw data. Currently supported modalities are:
- 10x single-cell ATAC (**10x**)
- 10x Single Cell Multiome ATAC + Gene Expression (**multiome**)
- Bing Ren lab split-pool scATAC (**ren**)

**Genome** specifies the target genome for alignment and analysis. The pipeline currently supports **GRCh38** (human) and **mm10** (mouse).

Note: Lines starting with `#` will be ignored.

## Global Configuration

Global configuration parameters are stored in `config.yaml`. (Visit [here](https://quickref.me/yaml) for a YAML cheat sheet.)

🚧 Under construction 🚧

### Compute configuration

### ENCODE DCC information

### Barcode correction parameters and whitelists

### Genome information

### Fragment file parameters

### ArchR parameters
