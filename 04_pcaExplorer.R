library(pcaExplorer)

# used vst cause normalized
dds_pcaE <- readRDS("files_from_pipeline/dds.rds")

pcaExplorer(dds_pcaE)
# NB: very slow package
# NB: annotation is a 2 col csv with "gene ID" and "other possible gene names"
# For model organism possible to extract extra info with PCA2GO