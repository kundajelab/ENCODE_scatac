library(ArchR)

# Disable HDF5 file locking
# Workaround for HDF5 I/O issues on shared filesystems
# https://github.com/GreenleafLab/ArchR/issues/248#issuecomment-789453997
Sys.setenv("HDF5_USE_FILE_LOCKING" = "FALSE")
Sys.setenv("RHDF5_USE_FILE_LOCKING" = "FALSE")

build_archr_project <- function(arrow_sample_names, input_paths, output_paths, threads, log_paths, seed, genome) {
    set.seed(seed)

    addArchRThreads(threads = threads)
    addArchRGenome(genome)

    input_paths = unlist(input_paths)
    arrow_output_names = paste0(output_paths[["arrows_temp_dir"]], arrow_sample_names)
    # print(input_paths) ####
    # print(arrow_output_names) ####
    arrows <- createArrowFiles(
        inputFiles = input_paths,
        sampleNames = arrow_sample_names,
        outputNames = arrow_output_names,
        # minTSS = 1, ####
        # minFrags = 1, ####
        addTileMat = TRUE,
        addGeneScoreMat = TRUE,
        force = TRUE,
        subThreading = FALSE, # required or else file locking gets turned back on
        logFile = log_paths[["arrow_create"]],
        QCDir = output_paths[["qc_dir"]]
    )
    # print(arrows) ####

    proj <- ArchRProject(
        ArrowFiles = arrows, 
        outputDirectory = output_paths[["project_dir"]],
        copyArrows = FALSE 
    )
    print(nCells(input = proj)) ####

    # Calculate doublet scores
    doub_scores <- addDoubletScores(
        input = proj,
        k = 10, 
        knnMethod = "UMAP", 
        LSIMethod = 1,
        logFile = log_paths[["doublets"]],
    )
    print(nCells(input = proj)) ####

    
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

    # Calculate motif enrichment
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
        dropCells = TRUE,
        load = FALSE,
        logFile = log_paths[["save"]],
    )

}

build_archr_project(snakemake@params[["sample_names"]], snakemake@input, snakemake@output, snakemake@threads, snakemake@log, snakemake@params[["seed"]], snakemake@params[["genome"]])