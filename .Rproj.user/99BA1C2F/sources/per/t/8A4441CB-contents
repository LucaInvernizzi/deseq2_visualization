# original github: https://github.com/kevinblighe/E-MTAB-6141

library(DESeq2)
library(ComplexHeatmap)
library(cluster)
library(circlize) # for color ramp in heatmap

dds_hm <- readRDS("files_from_pipeline/dds.rds") # to get overexpressed genes
vst_hm <- readRDS("files_from_pipeline/vst.rds") # vst is homoscedastic

# get list of overexpressed genes, used stv vs ctrl for example
results_hm <- results(dds_hm, 
                      contrast = c("description", "feStarved", "control"))
res_stv_ctrl_hm <- subset(results_hm, log2FoldChange > 1)
deg_hm <- row.names(res_stv_ctrl_hm)

# get expression matrix, homoscedasticity is needed for this analysis
mat_hm <- assay(vst_hm)
mat_deg_hm <- mat_hm[deg_hm,]

# z score normalization
heat <- t(scale(t(mat_deg_hm)))

# some graphicaldetails
myCol <- colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
myBreaks <- seq(-3, 3, length.out = 100)

# some clustering to represent some data
pamClusters <- cluster::pam(heat, k = 4) # pre-select k = 4 centers
pamClusters$clustering <- paste0('Cluster ', pamClusters$clustering)
# fix order of the clusters to have 1 to 4, top to bottom
pamClusters$clustering <- factor(pamClusters$clustering,
                                 levels = c('Cluster 1', 'Cluster 2', 
                                            'Cluster 3', 'Cluster 4'))

# Actual heatmap
hmap <- Heatmap(heat,
                
                # split the genes / rows according to the PAM clusters
                split = pamClusters$clustering,
                cluster_row_slices = FALSE,
                
                name = 'Gene\nZ-\nscore',
                
                col = colorRamp2(myBreaks, myCol),
                
                # parameters for colour-bar representing gradient of expression
                heatmap_legend_param = list(
                        color_bar = 'continuous',
                        legend_direction = 'vertical',
                        legend_width = unit(8, 'cm'),
                        legend_height = unit(5.0, 'cm'),
                        title_position = 'topcenter',
                        title_gp=gpar(fontsize = 12, fontface = 'bold'),
                        labels_gp=gpar(fontsize = 12, fontface = 'bold')),
                
                # row (gene) parameters
                cluster_rows = TRUE,
                show_row_dend = TRUE,
                #row_title = 'Statistically significant genes',
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
                row_title_rot = 90,
                show_row_names = FALSE,
                row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                row_names_side = 'left',
                row_dend_width = unit(25,'mm'),
                
                # column (sample) parameters
                cluster_columns = TRUE,
                show_column_dend = TRUE,
                column_title = '',
                column_title_side = 'bottom',
                column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                column_title_rot = 0,
                # show_column_names = FALSE, (removes col names if hm too big)
                column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                column_names_max_height = unit(10, 'cm'),
                column_dend_height = unit(25,'mm'),
                
                # cluster methods for rows and columns
                clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                clustering_method_columns = 'ward.D2',
                clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                clustering_method_rows = 'ward.D2'
                )

# visualize heatmap
draw(hmap,
     heatmap_legend_side = 'left',
     annotation_legend_side = 'right',
     row_sub_title_side = 'left')
