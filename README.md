# ENCODE scATAC Pipeline

**Note: This pipeline is currently a work in progress.**

A pipeline for processing raw single-cell ATAC datasets on the ENCODE portal.

Information on the specific analysis steps can be found in the [pipeline specification document](https://docs.google.com/document/u/2/d/e/2PACX-1vTlgtT4WeXbvRicybUHXnhZs8RKyB4EkTbcWooQ6qBxxQ_zIHpFEVHy38D5lC_s8_YDGfUTsyomJcs3/pub).

## Requirements

- A Linux-based OS
- A conda-based Python 3 installation
- [Snakemake v6.6.1+](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (full installation)
- An ENCODE DCC account with access to the necessary datasets

Additional requirements for cloud execution:
- [Kubectl](https://kubernetes.io/docs/tasks/tools/install-kubectl-linux/)
- A cloud provider CLI for Kubernetes cluster creation
- A cloud provider CLI for remote storage (if different from above)

All other dependencies are handled by the pipeline itself

## Local Execution

1. Install any necessary requirements above
2. Download the pipeline
    ```
    git clone https://github.com/kundajelab/ENCODE_scatac
    ```
3. Activate the `snakemake` conda environment:
    ```
    conda activate snakemake
    ```
4. Configure the pipeline in the `/config` directory. Detailed information can be found [here](config/README.md).
5. Run the pipeline:
    ```
    snakemake --use-conda --cores $NCORES 
    ```
    Here, `$NCORES` is the number of cores to utilize

Note: When run for the first time, the pipeline will take some time to install conda packages. 

## Cloud Execution with Kubernetes

1. Install and configure the pipeline as specified above
2. Create a cloud cluster. Note that setup specifics may differ depending on the cloud provider. Example setup instructions for [GCP] (https://snakemake.readthedocs.io/en/stable/executing/cloud.html#setup-kubernetes-on-google-cloud-engine) and for [Azure] (https://snakemake.readthedocs.io/en/stable/executor_tutorial/azure_aks.html#create-an-auto-scaling-kubernetes-cluster).
3. Configure remote storage. Instructions for each provider can be found [here] (https://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html). For our purpose, only the environment variables and command line configuration are needed.
4. Run the pipeline:
    ```
    snakemake --kubernetes --use-conda --default-remote-provider $REMOTE --default-remote-prefix $PREFIX --jobs $NJOBS --envvars $VARS
    ```
    Here:
        * `$REMOTE` is the cloud storage provider, and should be one of `{S3,GS,FTP,SFTP,S3Mocked,gfal,gridftp,iRODS,AzBlob,XRootD}`
        * `$PREFIX` is the target bucket name or subfolder in storage
        * `$NJOBS` is the maximum number of jobs to be run in parallel
        * `$VARS` is a list of environment variables for accessing remote storage. The `--envvars` flag can be omitted if no variables are required.

## Additional Execution Modes

This pipeline has been tested locally and on the cloud via Kubernetes. However, Snakemake offers a number of additional execution modes.

[Documentation on cluster execution](https://snakemake.readthedocs.io/en/stable/executing/cluster.html)

[Documentation on cloud execution](https://snakemake.readthedocs.io/en/stable/executing/cloud.html)