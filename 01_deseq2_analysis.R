library(DESeq2)

sampledescription <- read.csv2("files_from_pipeline/sampledescription.csv",
                               header = TRUE, sep = ",", row.names = 1)

counts_tab <- read.table("files_from_pipeline/counts_tab.txt",
                         header = TRUE, sep = " ", row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = counts_tab,
                              colData = sampledescription,
                              design = ~ description)

dds <- DESeq(dds)
vst <- vst(dds)

results <- results(dds, contrast = c("description", "feStarved", "control"))
res_stv_ctrl <- subset(results, log2FoldChange > 1)
res_str_ctrl_order <- res_stv_ctrl[order(res_stv_ctrl$log2FoldChange), ]

head(res_str_ctrl_order)
tail(res_str_ctrl_order)
