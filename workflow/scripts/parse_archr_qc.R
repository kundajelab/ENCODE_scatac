in_path_doublets <- snakemake@input[["qc_ds_data"]]
in_path_meta <- snakemake@input[["qc_meta"]]
out_path_doublets <- snakemake@output[["qc_ds_data"]]
out_path_meta <- snakemake@output[["qc_meta"]]

r <- readRDS(in_path_meta)
write.table(r, out_path_meta, sep = '\t', row.names = FALSE, quote = FALSE)

r <- readRDS(in_path_doublets)
d <- data.frame(
    barcode = as.vector(names(r[['doubletResults']][['doubletEnrichLSI']])), 
    doubletEnrichLSI = r[['doubletResults']][['doubletEnrichLSI']], 
    doubletScoreLSI = r[['doubletResults']][['doubletScoreLSI']]
)
write.table(d, out_path_doublets, sep = '\t', row.names = FALSE, quote = FALSE)