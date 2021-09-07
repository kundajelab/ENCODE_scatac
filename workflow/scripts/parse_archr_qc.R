console_log <- file(snakemake@log[[1]], open = "wt")
# sink(console_log)
# sink(console_log, type = "message")

in_path_doublets <- snakemake@input[["qc_ds_data"]]
in_path_meta <- snakemake@input[["qc_meta"]]
out_path_doublets <- snakemake@output[["qc_ds_data"]]
out_path_meta <- snakemake@output[["qc_meta"]]

r <- readRDS(in_path_meta)
write.table(r, out_path_meta, sep = '\t', row.names = FALSE, quote = FALSE)

r <- readRDS(in_path_doublets)
print(r$simulatedDoubletUMAP) ####
print(r$originalDataUMAP) ####
res <- r[['doubletResults']]
d <- data.frame(
    # barcode = as.vector(names(res[['doubletEnrichLSI']])), 
    doubletEnrichLSI = res[['doubletEnrichLSI']], 
    doubletScoreLSI = res[['doubletScoreLSI']]
)
# write.table(d, out_path_doublets, sep = '\t', row.names = FALSE, quote = FALSE)

# sink(type = "message")
# sink()