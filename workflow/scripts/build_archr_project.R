Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

# library(devtools)

# devtools::install_github("GreenleafLab/ArchR@2be1294eb1fbff364fb538d71f4f545ee4384d09", ref="master", repos = BiocManager::repositories(), dependencies = FALSE, INSTALL_opts = '--no-lock')
# devtools::install_github("GreenleafLab/chromVARmotifs@38bed559c1f4770b6c91c80bf3f8ea965da26076", ref="master", repos = BiocManager::repositories(), dependencies = FALSE, INSTALL_opts = '--no-lock')

library(ArchR)
library(chromVARmotifs)

addArchRVerbose(verbose = FALSE)

# Disable HDF5 file locking
# Workaround for HDF5 I/O issues on NFS
# https://github.com/GreenleafLab/ArchR/issues/248#issuecomment-789453997
Sys.setenv("HDF5_USE_FILE_LOCKING" = "FALSE")
Sys.setenv("RHDF5_USE_FILE_LOCKING" = "FALSE")

build_archr_project <- function(arrow_sample_name, input_path, output_paths, threads, log_paths, genome, seed) {
    set.seed(seed)

    addArchRThreads(threads = threads)
    addArchRGenome(genome)

    # input_paths = unlist(input_paths)]
    # print(output_paths) ####
    # print(output_paths[["arrows_temp_dir"]]) ####
    arrow_output_dir = "arrowtmp/"
    arrow_output_name = paste(arrow_output_dir, "/", arrow_sample_name,  sep = "")
    dir.create(arrow_output_dir)
    # print(input_paths) ####
    # print(arrow_output_names) ####
    arrows <- createArrowFiles(
        inputFiles = c(input_path),
        sampleNames = c(arrow_sample_name),
        outputNames = c(arrow_output_name),
        # minTSS = 1, ####
        # minFrags = 1, ####
        offsetPlus = 0,
        offsetMinus = 0,
        addTileMat = TRUE,
        addGeneScoreMat = TRUE,
        force = TRUE,
        subThreading = FALSE, # required for no file locking
        logFile = log_paths[["arrow_create"]],
        QCDir = output_paths[["qc_dir"]]
    )
    # print(arrows) ####

    # Calculate doublet scores
    doub_scores <- addDoubletScores(
        input = arrows,
        k = 10, 
        knnMethod = "UMAP", 
        LSIMethod = 1,
        outDir = output_paths[["qc_dir"]],
        logFile = log_paths[["doublets"]],
    )

    # Create project
    dir.create(output_paths[["project_dir"]])
    # setwd(output_paths[["project_dir"]])
    proj <- ArchRProject(
        ArrowFiles = arrows, 
        outputDirectory = output_paths[["project_dir"]],
        copyArrows = FALSE 
    )
    markers_dir = file.path(output_paths[["project_dir"]], "Markers")
    dir.create(markers_dir)
    
    # Filter doublets
    proj <- filterDoublets(proj)

    # Conduct LSI dimensionality reduction
    proj <- addIterativeLSI(
        ArchRProj = proj,
        useMatrix = "TileMatrix", 
        name = "IterativeLSI", 
        logFile = log_paths[["lsi"]]
    )

    # Cluster cells by LSI values
    proj <- addClusters(
        input = proj,
        reducedDims = "IterativeLSI",
        method = "Seurat",
        name = "Clusters",
        resolution = 0.8,
        logFile = log_paths[["cluster"]]
    )

    # Calculate UMAP coordinates from LSI values
    proj <- addUMAP(
        ArchRProj = proj, 
        reducedDims = "IterativeLSI", 
        name = "UMAP", 
        nNeighbors = 30, 
        minDist = 0.5, 
        metric = "cosine",
    )

    # Get marker genes
    marker_genes <- getMarkerFeatures(
        ArchRProj = proj, 
        useMatrix = "GeneScoreMatrix", 
        groupBy = "Clusters",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon",
        logFile = log_paths[["marker_genes"]]
    )
    marker_genes_list <- getMarkers(marker_genes, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
    marker_genes_path <- file.path(markers_dir, "marker_genes.rds")
    saveRDS(marker_genes_list, marker_genes_path)

    # Generate pseudo-bulk replicates
    proj <- addGroupCoverages(
        ArchRProj = proj, 
        groupBy = "Clusters",
        logFile = log_paths[["pseudobulk_rep"]]
    )

    # Call peaks
    proj <- addReproduciblePeakSet(
        ArchRProj = proj, 
        groupBy = "Clusters", 
        pathToMacs2 = findMacs2(),
        logFile = log_paths[["peak_call"]]
    )

    # Build cell-peak matrix
    proj <- addPeakMatrix(
        ArchRProj = proj,
        logFile = log_paths[["peak_matrix"]]
    )

    # Get marker peaks
    marker_peaks <- getMarkerFeatures(
        ArchRProj = proj, 
        useMatrix = "PeakMatrix", 
        groupBy = "Clusters",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon",
        logFile = log_paths[["marker_peaks"]]
    )
    marker_peaks_list <- getMarkers(marker_peaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
    marker_peaks_path <- file.path(markers_dir, "marker_peaks.rds")
    saveRDS(marker_peaks_list, marker_peaks_path)

    # Calculate motif enrichment
    proj <- addMotifAnnotations(
        ArchRProj = proj, 
        motifSet = "cisbp", 
        name = "Motif",
        logFile = log_paths[["fetch_motif"]]
    )
    enrichMotifs <- peakAnnoEnrichment(
        seMarker = marker_peaks,
        ArchRProj = proj,
        peakAnnotation = "Motif",
        cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
        logFile = log_paths[["enr_motif"]]
    )

    # Calculate ENCODE TF binding site enrichment
    proj <- addArchRAnnotations(
        ArchRProj = proj,
        db = "ArchR",
        collection = "EncodeTFBS",
        logFile = log_paths[["fetch_tf"]]
    )
    enrichTF <- peakAnnoEnrichment(
        seMarker = marker_peaks,
        ArchRProj = proj,
        peakAnnotation = "EncodeTFBS",
        cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
        logFile = log_paths[["enr_tf"]]
    )

    saveArchRProject(
        ArchRProj = proj,
        overwrite = TRUE,
        # dropCells = TRUE,
        load = FALSE,
        logFile = log_paths[["save"]],
    )

    tar_cmd = paste("tar -zcf", output_paths[["project_tar"]], output_paths[["project_dir"]])
    system(tar_cmd)

}

build_archr_project(snakemake@params[["sample_name"]], snakemake@input[["frag"]], snakemake@output, snakemake@threads, snakemake@log, snakemake@params[["genome"]], snakemake@params[["seed"]])

sink(type = "message")
sink()