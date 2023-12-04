library(DESeq2)

# extract genes
dds_ff <- readRDS("files_from_pipeline/dds.rds")

results_ff <- results(dds_ff, 
                      contrast = c("description", "feStarved", "control"))
res_stv_ctrl_ff <- subset(results_ff, log2FoldChange > 2)

# modify names for ff format
gene_vector_ff <- rownames(res_stv_ctrl_ff) 

gene_vector_ff <- gsub("Afu", "AFUA_", gene_vector_ff)
gene_vector_ff <- toupper(gene_vector_ff)

# filter for genes present in GO on ff

names_ff <- read.csv2("results/fungi_fun/fungi_fun_files/go_aspergillus_fumigatus_af293_gene_ids.txt",
                      skip = 2,
                      header = FALSE,
                      sep = "\t")
vector_names_ff <- names_ff[,1]

intersect_ff <- intersect(vector_names_ff, gene_vector_ff)

# print txt
write(intersect_ff, file = "results/fungi_fun/fungi_fun_files/genes_ff.txt")
#NB: file has empty line at the end not read by ff, delete it and it works