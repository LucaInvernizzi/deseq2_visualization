library(DESeq2)
library(VennDiagram)

dds_vd <- readRDS("files_from_pipeline/dds.rds")

# contrast for the 3 condition vs control
res_stv_ctrl_vd <- results(dds_vd, 
                           contrast = c("description", 
                                        "feStarved", "control"))
res_h2o2_ctrl_vd <- results(dds_vd, 
                            contrast = c("description", 
                                         "h2o2Treated", "control"))
res_double_ctrl_vd <- results(dds_vd, 
                              contrast = c("description", 
                                           "feStarved_h2o2Treated", "control"))
# up and downregulated genes
stv_ctrl_up_vd <- subset(res_stv_ctrl_vd, log2FoldChange > 2)
stv_ctrl_up_vd <- rownames(stv_ctrl_up_vd) 
stv_ctrl_dw_vd <- subset(res_stv_ctrl_vd, log2FoldChange < -2)
stv_ctrl_dw_vd <- rownames(stv_ctrl_dw_vd)

h2o2_ctrl_up_vd <- subset(res_h2o2_ctrl_vd, log2FoldChange > 2)
h2o2_ctrl_up_vd <- rownames(h2o2_ctrl_up_vd)
h2o2_ctrl_dw_vd <- subset(res_h2o2_ctrl_vd, log2FoldChange < -2)
h2o2_ctrl_dw_vd <- rownames(h2o2_ctrl_dw_vd)

double_ctrl_up_vd <- subset(res_double_ctrl_vd, log2FoldChange > 2)
double_ctrl_up_vd <- rownames(double_ctrl_up_vd)
double_ctrl_dw_vd <- subset(res_double_ctrl_vd, log2FoldChange < -2)
double_ctrl_dw_vd <- rownames(double_ctrl_dw_vd)

# build diagram of only upregulated genes
venn.diagram(
        x = list(stv_ctrl_up_vd, h2o2_ctrl_up_vd, double_ctrl_up_vd),
        category.names = c("stv_ctrl" , "h2o2_ctrl" , "double_ctrl"),
        filename = 'results/venn/Upregulated Genes Venn Diagram.png',
        output=TRUE
)
