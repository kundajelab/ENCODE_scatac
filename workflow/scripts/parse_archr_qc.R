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
# print(r$simulatedDoubletUMAP) ####
# print(length(r$originalDataUMAP$score)) ####
# print(length(r$originalDataUMAP$enrichment)) ####
# print(length(r$doubletResults$doubletScoreUMAP)) ####
# print(length(r$doubletResults$doubletEnrichUMAP)) ####
# print(length(r$doubletResults$doubletEnrichLSI)) ####
# print(length(r$doubletResults$doubletScoreLSI)) ####
print(names(r$originalDataUMAP)) ####
print(head(r$originalDataUMAP)) ####
print(row.names(r$originalDataUMAP)) ####

res <- r[['doubletResults']]
d <- data.frame(
    # barcode = as.vector(names(r$originalDataUMAP)), 
    doubletScoreUMAP = res$doubletScoreUMAP,
    doubletEnrichUMAP = res$doubletEnrichUMAP,
    doubletEnrichLSI = res$doubletEnrichLSI, 
    doubletScoreLSI = res$doubletScoreLSI
)
# print(d) ####
# write.table(d, out_path_doublets, sep = '\t', row.names = FALSE, quote = FALSE)

# sink(type = "message")
# sink()