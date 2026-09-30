#######LckCre mTmG LRH-1 analysis
#THYMUS_TCRb


library(nlme)
library(ggprism)
library(ggpubr)
library(dbplyr)
library(BH)
library(RcppArmadillo)
library(unikn)
library(gridExtra)
library(umap)
library(ggplot2)
library(reshape2)
library(flowCore)
library(igraph)
library(FlowSOM)
library(CATALYST)
library(scales)
library(ConsensusClusterPlus)
library(readxl)
library(grid)
library(ggthemes)
library(ComplexHeatmap)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(gplots)
library(heatmaply)
library(ncdfFlow)
library(tidyr)
library(matrixStats)
library(gghighlight)
library(tidyverse)
library(rstatix)
library(openxlsx)

md <- read_excel("metadata.xlsx")
md

#Define condition variables as named in metadata
md$group <- factor(md$group, levels = c("LRH1 KO","Ctr"))
head(data.frame(md))
tail(data.frame(md))

#Define colors for organ and genotype
color_group <- c("#c07631","#121212")
barplot(1:2, col = color_group)
names(color_group) <- levels(md$group)

  
#Load the content of the .fcs files into R
fcs_raw <- read.ncdfFlowSet(md$file_name, transformation = FALSE, truncate_max_range = FALSE)
colnames(fcs_raw)
fcs_raw

#Load the content of the channel.names.xlsx file into R
channel.names <- read_excel("channel.names.xlsx")
channel.names
new.names <- channel.names$fcs_colname
colnames(fcs_raw) <- new.names
panel_filename <- "channel.names.xlsx"
panel <- read_excel(panel_filename)
panel

#Replace problematic characters (to eliminate potential error source)
panel$fcs_colname <- gsub("-", "_", panel$fcs_colname)
panel_fcs <- pData(parameters(fcs_raw[[1]]))
panel_fcs$desc <- gsub("-", "_", panel_fcs$name)

#Define Lineage markers (or other classification)
lineage_markers <- panel$Antigen[panel$`marker_class`=="Lineage"] 
lineage_markers

#Define Funcional markers
functional_markers <- panel$Antigen[panel$`marker_class`=="Function"]
functional_markers

#Define All markers in your panel (for global analysis)
all_markers <- panel$Antigen

#Column subsetting
all_markers <- channel.names$fcs_colname
fcs <- fsApply((fcs_raw[, c(all_markers)]), exprs)
class(fcs)
head(fcs)
dim(fcs)
min(fcs)
max(fcs)

#Transformation
asinh_scale_200 <- c("CD8a")
cofactor_1 <- 200
fcs[,asinh_scale_200] <- asinh(fcs[,asinh_scale_200] / cofactor_1)

asinh_scale_200 <- c("CD25")
cofactor_2 <- 200
fcs[,asinh_scale_200] <- asinh(fcs[,asinh_scale_200] / cofactor_2)

asinh_scale_150 <- c("GFP")
cofactor_3 <- 150
fcs[,asinh_scale_150] <- asinh(fcs[,asinh_scale_150] / cofactor_3)

asinh_scale_4500 <- c("CD44")
cofactor_4 <- 4500
fcs[,asinh_scale_4500] <- asinh(fcs[,asinh_scale_4500] / cofactor_4)

asinh_scale_2000 <- c("dtTomato")
cofactor_5 <- 2000
fcs[,asinh_scale_2000] <- asinh(fcs[,asinh_scale_2000] / cofactor_5)

asinh_scale_2000 <- c("CD4")
cofactor_6 <- 2000
fcs[,asinh_scale_2000] <- asinh(fcs[,asinh_scale_2000] / cofactor_6)

asinh_scale_300 <- c("CD28")
cofactor_7 <- 300
fcs[,asinh_scale_300] <- asinh(fcs[,asinh_scale_300] / cofactor_7)

asinh_scale_300 <- c("CD27")
cofactor_8 <- 300
fcs[,asinh_scale_300] <- asinh(fcs[,asinh_scale_300] / cofactor_8)

asinh_scale_200 <- c("TCRb")
cofactor_9 <- 200
fcs[,asinh_scale_200] <- asinh(fcs[,asinh_scale_200] / cofactor_9)



# Transform the markers to values between 0-1 (Normalization for FlowSOM algorithm)
rng <- colQuantiles(fcs, probs = c(0.01, 0.99))
expr01 <- t((t(fcs) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1

# Generate sample IDs corresponding to each cell in the 'expr' matrix
sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))
sample_ids
ggdf <- data.frame(sample_id = sample_ids, fcs)
ggdf <- melt(ggdf, id.var = "sample_id",
             value.name = "expression", variable.name = "antigen")
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$group <- md$group[mm]
head(ggdf)
tail(ggdf)

theme_bw2 <- function (base_size = 14, base_family = "", face="bold") 
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.background = element_rect(fill = "white", 
                                          colour = NA), panel.border = element_rect(fill = "NA", 
                                                                                    colour = "NA"), panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "NA", 
                                          size = 0.25), strip.background = element_rect(fill = "NA", 
                                                                                        colour = "NA"), legend.key = element_rect(fill = "white", 
                                                                                                                                  colour = NA), complete = FALSE)
}
plot1 <- ggplot(ggdf, aes(x = expression, color = group, group = sample_id)) +
  geom_density(size = 0.9) +
  facet_wrap(~ antigen, ncol = 3, scales = "free") + theme_bw2() +
  theme(axis.text.x = element_text(angle = 0), 
        axis.title.y.left = element_text(vjust = 2.7, size = 18, face ="bold"), 
        axis.title.x = element_text(size = 18),
        strip.text = element_text(size = 18), 
        axis.text = element_text(size = 18, face = "bold"), 
        legend.position = "right",
        legend.text = element_text(size = 16), 
        legend.title = element_blank() ) +
  guides(color = guide_legend(ncol = 1)) +
  scale_color_manual(values = color_group)
plot1

### data normalization for UMAP analysis
# Adjust data to start from (approximately) 0
q.vector <- apply(fcs[,all_markers], 2, function(x) quantile(x, 0.01, names = F))
q.vector
data.shift <- fcs
data.shift[,all_markers] <- sweep(data.shift[,all_markers], 2, q.vector)
head(data.shift)

# Normalize to have everything from 0 to 1
per.vector <- apply(data.shift[,all_markers], 2, function(x) quantile(x, 0.99999, names = F))
per.vector
data.shift[,all_markers] <- t(t(data.shift[,all_markers]) / as.numeric(per.vector))
dim(data.shift)
tail(data.shift)

# Check whether the range is approximately from 0 to 1
apply(data.shift, 2, min)
apply(data.shift, 2, max)

# Find and skip duplicates
dups <- which(!duplicated(data.shift[, lineage_markers]))
dups

#plot gene expression of all markers to check for negative values
normdf <- data.frame(sample_id = sample_ids, data.shift)
tail(normdf)
normdf <- melt(normdf, id.var = "sample_id",
               value.name = "expression", variable.name = "antigen")
mm <- match(normdf$sample_id, md$sample_id)
normdf$group <- md$group[mm]
tail(normdf)

norm_exp <- ggplot(normdf, aes(x = expression, color = group, group = sample_id)) +
  geom_density(size = 0.9) +
  facet_wrap(~ antigen, ncol = 3, scales = "free") + theme_bw2() +
  theme(axis.text.x = element_text(angle = 0), 
        axis.title.y.left = element_text(vjust = 2.7, size = 16, face ="bold"), 
        axis.title.x = element_text(size = 16),
        strip.text = element_text(size = 16), 
        axis.text = element_text(size = 16, face = "bold"), 
        legend.position = "right",
        legend.text = element_text(size = 16), 
        legend.title = element_blank() ) +
  guides(color = guide_legend(ncol = 1)) +
  scale_color_manual(values = color_group)
norm_exp

# Data subsampling: create indices by sample
inds <- split(1:length(sample_ids), sample_ids)
length(sample_ids)
length(data.shift)

# Subsample (cells per sample) for dimenionality reduction
tsne_ncells <- pmin(table(sample_ids), 5000)
tsne_ncells

# Get subsampled indices without dups
set.seed(1234)
tsne_inds <- lapply(names(inds), function(i){
  s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)
  intersect(s, dups)
})
length(dups)
length(tsne_inds)
length(data.shift)
expr <- as.data.frame(data.shift)
expr$cell_id <- 1:nrow(expr)
tsne_inds <- unlist(tsne_inds)

#Add expression data to subsamples acc. to cell indices
tsne_expr_func <- expr[tsne_inds, functional_markers]
tsne_expr <- expr[tsne_inds, lineage_markers]
tsne_expr
    
#Create UMAP plot  
data_umap <- umap(tsne_expr, random_state=1234, verbose = T)
umap <- as.data.frame(data_umap$layout)
dim(umap) 
head(umap)
tail(umap)
colnames(umap) <- c("UMAP1", "UMAP2")

### plot UMAP black
plot_umap <- ggplot(umap, aes(x = UMAP1, y = UMAP2)) +
  geom_point(size = 0.5) +
  coord_fixed(ratio = 1)
plot_umap 

### plot UMAP with expression overlayed
data_plot_u <- cbind(tsne_expr, umap)
head(data_plot_u)
data_melt_u <- melt(data_plot_u, id.vars = c("UMAP1", "UMAP2"), variable.name = "antigen", value.name = "Expression")
tail(data_melt_u)

color_grad_flow2 <- c("#331820", "#4d4f55", "#55626b", "#5a767e", "#628b8a", "#709b91", "#82aa96", "#98b89a", "#b0c6a2", "#c9d3ab", "#e4e0b6", "#feedc3")

plot_umap_expr <- ggplot(data_melt_u, aes(x = UMAP1, y = UMAP2, color = Expression)) +
  geom_point(size = 0.01) +
  coord_fixed(ratio = 1) +
  scale_colour_gradientn(colours = color_grad_flow2, limits = c(0,1)) +
  facet_wrap(~ antigen, nrow = 2) +
  theme_few() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
plot_umap_expr

#Combine with samples
#combination umap+expression with sample IDs
tsne_expr$cell_id <- expr[tsne_inds,c("cell_id")]
head(tsne_expr)
tsne_expr$sample_id <- sample_ids[tsne_inds]
data_plot_u$sample_id <- sample_ids[tsne_inds]

data_plot_u$cell_id <- expr[tsne_inds,c("cell_id")]
head(data_plot_u)
umap_exp <- cbind(tsne_expr[,c("cell_id", "sample_id")], data_plot_u[,c("UMAP1", "UMAP2")])
head(umap_exp)
#save(umap_exp, file="data_umap_coordinates_Th_4500.RDa")

#Generate ids for sample and cell tracing
dim.vec_surf <- as.numeric(dim(data.shift)[1])
class(dim.vec_surf)
gate.source_surf <- as.vector(x = NULL)
for(i in 1:length(dim.vec_surf)) {temp.source <- rep(i, dim.vec_surf[i])
gate.source_surf <- c(gate.source_surf, temp.source)}
cell_id_surf <- 1:dim(data.shift)[1]
head(cell_id_surf)

#add to master data
data.shift <- cbind(data.shift, gate.source_surf)
data.shift <- cbind(data.shift, cell_id_surf)
head(data.shift)
colnames(data.shift)[7:8] <- c("gate_source", "cell_id")
head(data.shift)

#Color panels
white.black100 <- colorRampPalette(c("white", "black"))(n = 100)

#########CLUSTERING
head(fcs)
fsom <- ReadInput(flowFrame(exprs = fcs, desc = list(FIL = 1)), transform = FALSE, scale = FALSE)
set.seed(1234)
som <- BuildSOM(fsom, colsToUse = lineage_markers)
head(som$map$mapping[,1])
labels <- som$map$mapping[,1]

# Metaclustering with ConsensusClusterPlus
codes <- som$map$codes
plot_outdir <- "consensus_plots"

#Metaclustering with defined number of groups to be generated
nmc <- 20
mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 100,
                           pItem = 0.9, pFeature = 1, title = plot_outdir, 
                           plot = "png", clusterAlg = "hc", innerLinkage = "average", 
                           finalLinkage = "average", distance = "euclidean", seed = 1234)

#Get cluster ids for each cell
code_clustering1 <- mc[[nmc]]$consensusClass
cell_clustering1 <- code_clustering1[som$map$mapping[,1]]

color_clusters <- c("#fc0b03","#fc7f03","#fcd303","#bafc03","#28fc03","#03fce3","#03adfc",
                    "#0349fc","#8803fc","#f003fc","#fc038c","#ce8c33","#657ec5",
                    "#d38977","#3d6b30","#33A02C", "#E6AB02", "#7570B3", "#BEAED4", "#666666", 
                    "#999999", "#f5424b", "#e8bec0", "#bee8d2", "#03592b", "#3a0359", "#592e03",
                    "#1bd1d1", "#fca4c0", "#f7e702", "#a9a1f7")
barplot(1:30, col=color_clusters)

draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.9, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)

plot_clustering_heatmap_wrapper <- function(fcs, expr01, 
                                            cell_clustering, color_clusters, 
                                            cluster_merging = NULL) {
#Calculate the median expression
expr_median <- data.frame(fcs, cell_clustering = cell_clustering) %>%
  group_by(cell_clustering) %>%
  summarize_each(funs(median))
expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>%
  group_by(cell_clustering) %>%
  summarize_each(funs(median))
  
#Calculate cluster frequencies
clustering_table <- as.numeric(table(cell_clustering))
  
#This clustering is based on the markers that were used for the main clustering
d <- dist(expr_median[, colnames(fcs)], method = "euclidean")
cluster_rows <- hclust(d, method = "average")
  
expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
rownames(expr_heat) <- expr01_median$cell_clustering
  
labels_row <- paste0(rownames(expr_heat), " (",
                       round(clustering_table / sum(clustering_table) * 100, 2), "%)")
labels_col <- colnames(expr_heat)
  
#Row annotation for the heatmap
annotation_row <- data.frame(cluster = factor(expr01_median$cell_clustering)) 
rownames(annotation_row) <- rownames(expr_heat)
  
color_clusters <- color_clusters[1:nlevels(annotation_row$cluster)]
names(color_clusters) <- levels(annotation_row$cluster)
annotation_colors <- list(cluster = color_clusters)
annotation_legend <- FALSE
  
if(!is.null(cluster_merging)){
   cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
   annotation_row$cluster_merging <- cluster_merging$new_cluster
   color_clusters <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
   names(color_clusters) <- levels(cluster_merging$new_cluster)
   annotation_colors$cluster_merging <- color_clusters
   annotation_legend <- TRUE
   }  

##Colors for the heatmap
#color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
pheatmap(expr_heat, 
          color = white.black100,
          cluster_cols = FALSE, 
          cluster_rows = cluster_rows,
          #breaks = seq(-0.05, 1, length.out = 100),
          legend_breaks = c(0, 1),
          labels_col = labels_col, 
          labels_row = labels_row,
          cellwidth = 50, 
          cellheight = 49,
          display_numbers = FALSE, 
          number_color = "black", 
          border_color = FALSE,
          fontsize = 25, fontsize_number = 0.005, fontsize_col = 20,
          annotation_row = annotation_row, 
          annotation_colors = annotation_colors,
          annotation_legend = annotation_legend)
}
clusterHM <- plot_clustering_heatmap_wrapper(fcs = fcs[, lineage_markers],
                                             expr01 = expr01[, lineage_markers],
                                             cell_clustering = cell_clustering1, 
                                             color_clusters = color_clusters)

#check heatmap and think about which cluster is what kind of population
mm <- match(umap_exp$sample_id, md$sample_id)
umap_exp$group <- md$group[mm]
umap_exp$exp <-md$Experiment[mm]
umap_exp$cell_clustering1 <- factor(cell_clustering1[tsne_inds], levels = 1:nmc)

#Plot k20 clusters on UMAP to see where the clusters are located and which clusters are close to each other
plotk20 <-ggplot(umap_exp, aes(x = UMAP1, y = UMAP2, color = cell_clustering1)) +
  geom_point(size = 0.05) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.title=element_blank(),
        legend.background = element_rect()) + coord_fixed(ratio = 1)+  
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4)))
plotk20

######### Manual metaclustering
#Load the content of the cluster_mergings.xlsx file into R (see Nowicka et al for description of "cluster_mergings")
cluster_merging_filename <- "cluster_mergings_5000.xlsx" 
cluster_merging1 <- read_excel(cluster_merging_filename) 
cluster_merging1
cluster_merging1$new_cluster

mm <- match(cell_clustering1, cluster_merging1$original_cluster)
cell_clustering1m <- cluster_merging1$new_cluster[mm]
length(cell_clustering1m)

# Update UMAP with the new annotated cell populations.
head(umap_exp)
umap_exp$cell_clustering1m <- factor(cell_clustering1m[tsne_inds])
data_plot_u$cell_clustering1m <- factor(cell_clustering1m[tsne_inds])

# Colors for the heatmap

color_merge <- c( "#a6e1f4",
                  "#398d9f",
                  "#4287f5",
                  "#8672fc",
                  "#324376",
                  "#5a3449",
                  "#111",
                  "#80391e",
                  "#8e2043"
                  
                  )


####Creation of  Heatmap of new defined Metaclustering
head(data_plot_u)
umap_sub <-data_plot_u[,c(lineage_markers, "cell_clustering1m")]
head(umap_sub)
tail(umap_sub)
umap_sub <- dplyr::filter(umap_sub, cell_clustering1m != "NA")
median_subset <- umap_sub %>%
  group_by(cell_clustering1m) %>%
  summarize_all(funs(median(., na.rm=TRUE)))

median_subset <- as.data.frame(median_subset)
median_subset_2 <- as.matrix(sapply(median_subset[, -1], as.numeric))
rownames(median_subset_2) <- median_subset$cell_clustering1m
median_subset_2
max(median_subset_2)
median_subset_2
lineage_markers
colnames(median_subset_2) <- lineage_markers 
colnames(median_subset_2)

annotation_row_1 <- data.frame("Cluster" = factor(median_subset$cell_clustering1m))
annotation_row_1$Cluster
class(annotation_row_1$Cluster)

newcol <- as.vector(annotation_row_1$Cluster)
newcol
annotation_row_1$Cluster <- as.factor(newcol)
annotation_row_1$Cluster
annotation_colors_1 <- list(Cluster = color_merge)
names(annotation_colors_1$Cluster) <- newcol
color_clustersmerge <- color_merge[1:nlevels(annotation_row_1$Cluster)]
color_clustersmerge

rownames(annotation_row_1) = rownames(median_subset_2)

  
#color <- colorRampPalette(brewer.pal(n = 10, name = "Blues", "Grey"))(10)
max(median_subset_2)
min(median_subset_2)
breaks <- seq(-0.0, 0.85,length.out =13)

color_grad_flow2 <- c("#331820","#331820","#331820","#4d4f55", "#55626b", "#5a767e", "#628b8a", "#709b91", "#82aa96", "#98b89a", "#b0c6a2", "#c9d3ab", "#e4e0b6", "#feedc3")



clustermergePH_1 <- pheatmap(median_subset_2, color = color_grad_flow2,
                              breaks = breaks,
                              scale = "none",
                              drop_levels = TRUE,
                              fontsize_number = 15, 
                              cluster_cols = TRUE, 
                              cluster_rows = FALSE, 
                              border_color = TRUE, 
                              legend_breaks = c(0, 0.3, 0.9),
                              fontsize = 30,
                              fontsize_row = 28,
                              fontsize_col = 28,
                              cellwidth = 50, 
                              cellheight = 50, 
                              annotation_row = annotation_row_1 ,
                              annotation_colors = annotation_colors_1,
                              annotation_legend = F,
                              legend =TRUE,
                              treeheight_col = 0,
                              angle_col = "45",
                              main = "Cluster identification")
clustermergePH_1


######### Project defined clustering on UMAP
## Project defined clustering on UMAP
data_umap_compl <- merge(umap_exp, data_plot_u[,c(lineage_markers, "cell_id")], by = "cell_id")
data_umap_compl <- dplyr::filter(data_umap_compl, cell_clustering1m != "NA")
head(data_umap_compl)
head(expr)
mm <- match(data_umap_compl$sample_id, md$sample_id)

plot_umap_clust <- ggplot(data_umap_compl, aes(x = UMAP1, y = UMAP2, color = cell_clustering1m)) +
  geom_point(size = 0.01) + coord_fixed(ratio = 1) +
  scale_color_manual(values = color_merge) +
  theme_light() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(colour=NA, size=0.8, fill="white"),
        strip.text = element_text(color="black", size=18, face="plain"),
        legend.title=element_blank(),
        legend.key.size = unit(0.7, "cm"),
        legend.key = element_rect(size = 5, fill="white", colour = NA),
        legend.text = element_text(size = 17, face = "plain")) +
  guides(color = guide_legend(override.aes = list(size = 4)))

plot_umap_clust


##PLOT GROUPS
#create count table
counts_table <- table(cell_clustering1m, sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table, header = TRUE)
props <- as.data.frame.matrix(props_table, header = TRUE)
props

#Calculate cluster frequency
ggdf1 <- melt(data.frame(cluster = rownames(props), props), id.vars = "cluster", value.name = "Frequency", variable.name = "sample_id")
ggdf1

levels(ggdf1$sample_id) <- sub("^X", "", levels(ggdf1$sample_id))

ggdf1$sample_id <- as.factor(ggdf1$sample_id)
ggdf1$sample_id <-as.numeric(paste(ggdf1$sample_id))

ggdf1 <- dplyr::filter(ggdf1, cluster != "NA")
mm <- match(ggdf1$sample_id, md$sample_id)
ggdf1$group <- factor(md$group[mm])
ggdf1

##PLOT genotype FREQUENCY  ## NOTE: need more than one mice per group
plotsub_genotyp <- ggplot(ggdf1) +
  geom_violin(aes(x = factor(group, level = c("Ctr", "LRH1 KO")), y = Frequency, color = group, fill = group), 
              alpha = 0.2, trim = TRUE, scale = "width") +
  geom_point(aes(x = factor(group, level = c("Ctr", "LRH1 KO")), y = Frequency, color = group),
             alpha = 1, size = 3, stroke = 1.0, 
             position = position_jitter(0.05)) + 
  stat_summary(aes(y=Frequency, x=group, color=group), fun.data = "mean_se", geom = "crossbar", width = 0.15,  colour = "black") + #stat_summary(aes(y=Counts, x=condition, color=condition),fun.data="mean_se", color = "black", size = 0.5) +
  ylab("% of cells") +
  ggtitle("Frequency of clusters in the thymus") +
  facet_wrap(~ cluster, scale = "free_y", nrow=2) + #facet_grid(~cluster) + scales ="free_y",
  theme_light() +
  theme(axis.text.x = element_text(size = 20, hjust = 1, vjust = 1, angle = 45, color = "black"), 
        axis.ticks.x = element_line(size = 1, color = "black"),
        axis.ticks.y = element_line(size = 1, color = "black"),
        axis.ticks.length.x = unit(.30, "cm"),
        axis.ticks.length.y = unit(.30, "cm"),
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y.left = element_text(size = 23, face = "plain", vjust = 2),
        strip.background = element_rect(colour=NA, size = 0.8, fill = "white"),
        panel.border = element_rect(colour = "black", size = 0.8),
        #strip.text = element_blank(), 
        strip.text = element_text(color = "black", size = 20, face="bold"),
        axis.line = element_line(size = 0.6, color = "black"),
        #panel.grid.major.x = element_blank(),
        legend.text = element_blank(), legend.title = element_blank(), legend.position = "none",
        plot.title = element_text(color="black", size=24, hjust = 0.5)) + 
  scale_y_continuous() +
  expand_limits(y=0) +
  scale_color_manual(values = color_group) + 
  scale_fill_manual(values = color_group) +
  NULL

plotsub_genotyp

write.xlsx(ggdf1,"Freq_cluster.xlsx")


####Tota cell number
ggdf1$cell_number <- factor(md$cell_number[mm])
ggdf3 <- ggdf1
head(ggdf3)
class(ggdf3$cell_number)
ggdf3$cell_number <- as.numeric(as.character(ggdf3$cell_number))
class(ggdf3$cell_number)

head(ggdf3)
tail(ggdf3)

ggdf3$total_cell_nr <- (ggdf3$cell_number / 100) *ggdf3$Frequency   # Creating new variable
head(ggdf3) 
tail(ggdf3)

ggdf3[,c("Frequency", "cell_number")] <- NULL
head(ggdf3)
class(ggdf3)

ggdf3 <- do.call(data.frame,                      # Replace Inf in data by NA
                     lapply(ggdf3,
                            function(x) replace(x, is.infinite(x), NA)))

ggdf3$total_cell_nr[is.na(ggdf3$total_cell_nr)] <- 0

ggdf3

##PLOT CELL CLUSTERS THAT SHOW dtTomato EXPRESSION
plotsub_cellCount <- ggplot(ggdf3) +
  geom_violin(aes(x = factor(group, level = c("Ctr", "LRH1 KO")), y = total_cell_nr, color = group, fill = group), 
              alpha = 0.2, trim = TRUE, scale = "width") +
  geom_point(aes(x = factor(group, level = c("Ctr", "LRH1 KO")), y = total_cell_nr, color = group),
             alpha = 1, size = 3, stroke = 1.0, 
             position = position_jitter(0.05)) + 
  stat_summary(aes(y=total_cell_nr, x=group, color=group), fun.data = "mean_se", geom = "crossbar", width = 0.15,  colour = "black") + #stat_summary(aes(y=Counts, x=condition, color=condition),fun.data="mean_se", color = "black", size = 0.5) +
  ylab("% of cells") +
  ggtitle("Total cell number of clusters in the thymus") +
  facet_wrap(~ cluster, scale = "free_y", nrow=2) + #facet_grid(~cluster) + scales ="free_y",
  theme_light() +
  theme(axis.text.x = element_text(size = 20, hjust = 1, vjust = 1, angle = 45, color = "black"), 
        axis.ticks.x = element_line(size = 1, color = "black"),
        axis.ticks.y = element_line(size = 1, color = "black"),
        axis.ticks.length.x = unit(.30, "cm"),
        axis.ticks.length.y = unit(.30, "cm"),
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y.left = element_text(size = 23, face = "plain", vjust = 2),
        strip.background = element_rect(colour=NA, size = 0.8, fill = "white"),
        panel.border = element_rect(colour = "black", size = 0.8),
        #strip.text = element_blank(), 
        strip.text = element_text(color = "black", size = 20, face="bold"),
        axis.line = element_line(size = 0.6, color = "black"),
        #panel.grid.major.x = element_blank(),
        legend.text = element_blank(), legend.title = element_blank(), legend.position = "none",
        plot.title = element_text(color="black", size=24, hjust = 0.5)) + 
  scale_y_continuous() +
  expand_limits(y=0) +
  scale_color_manual(values = color_group) + 
  scale_fill_manual(values = color_group) +
  NULL

plotsub_cellCount

write.xlsx(ggdf3,"Cellcount_cluster.xlsx")

### Expression of dtTomato and GFP in clustered UMAP
#Create UMAP plot 
tsne_expr_func
data_plot_func <- cbind(tsne_expr_func, umap)
head(data_plot_func)

##Combine with samples
#combination umap+expression with sample IDs
tsne_expr_func$cell_id <- expr[tsne_inds,c("cell_id")]
head(tsne_expr_func)
tsne_expr_func$sample_id <- sample_ids[tsne_inds]
expr$sample_id <- sample_ids
tail(expr)
data_plot_func$sample_id <- sample_ids[tsne_inds]
data_plot_func$cell_id <- expr[tsne_inds,c("cell_id")]
head(data_plot_func)
umap_exp_func <- cbind(tsne_expr_func[,c("cell_id", "sample_id")], data_plot_func[,c("UMAP1", "UMAP2")])
head(umap_exp_func)

mm_func <- match(umap_exp_func$sample_id, md$sample_id)
umap_exp_func$group <- md$group[mm_func]
umap_exp_func$cell_clustering1 <- factor(cell_clustering1[tsne_inds], levels = 1:nmc)
head(umap_exp_func)

umap_exp_func$cell_clustering1m <- factor(cell_clustering1m[tsne_inds])
data_plot_func$cell_clustering1m <- factor(cell_clustering1m[tsne_inds])

data_plot_func$genotype <- md$group[mm_func]
head(data_plot_func)
tail(data_plot_func)

##make data longer (tidy)
data_plot_func_tidy <- data_plot_func %>% 
  pivot_longer(GFP:dtTomato, names_to = "antigen", values_to = "expression") %>% 
  arrange(antigen)
head(data_plot_func_tidy)
tail(data_plot_func_tidy)

max(data_plot_func_tidy$expression)


    ### plot all marker vs each other in order to evaluate transformation quality ("facs blots")
    # subsample cells for plotting
    n_cells <- 42000
    n <- nrow(data.shift)
    set.seed(1234)
    ix <- sample(1:n, n_cells)
    length(ix)
    sub_surf <- rbind(data.shift[ix,])
    s_surf <- data.frame(sub_surf)
    color.scale <- c("#1B9E77","gold" , "#7570B3", "#E7298A", "#D95F02" , "lightpink1", "#666666",
                     "#57B0FF", "#66CD00", "#1F78B4", "#B2DF8A", "firebrick1", 
                     "seagreen1", "brown" , "#6A3D9A", "mediumorchid", "grey" ,
                     "#BC80BD", "#80B1D3")
    head(data.shift)
    plot_surf_list <- list()
    for(i in 1:length(colnames(data.shift))){
      plot_surf <- ggplot(s_surf, aes_string(x = colnames(s_surf)[i], y = "GFP")) +
        geom_hex(bins = 100) +
        xlim(-0.1, 1.0) +
        ylim(-0.1, 1.0) +
        coord_fixed(ratio = 1) +
        scale_fill_gradientn(colours = color.scale, trans = "sqrt") +
        theme_bw() + theme(legend.position="none", axis.text = element_text(size = 5), axis.title = element_text(size = 8))
      plot_surf_list[[i]] <-  plot_surf
    }
    
    m1 <- arrangeGrob(grobs = plot_surf_list, ncol = 1, top=textGrob("GFP", gp=gpar(fontsize=6)))
    #for example
 
  
    for(i in 1:length(colnames(data.shift))){
      plot_surf <- ggplot(s_surf, aes_string(x = colnames(s_surf)[i], y = "dtTomato")) +
        geom_hex(bins = 100) +
        xlim(-0.1, 1.0) +
        ylim(-0.1, 1.0) +
        coord_fixed(ratio = 1) +
        scale_fill_gradientn(colours = color.scale, trans = "sqrt") +
        theme_bw() + theme(legend.position="none", axis.text = element_text(size = 5), axis.title = element_text(size = 8))
      plot_surf_list[[i]] <-  plot_surf
    }
    
    
    
    m2 <- arrangeGrob(grobs = plot_surf_list, ncol = 1, top=textGrob("dtTomato", gp=gpar(fontsize=6)))
    #for example
 

#########Plot dtTomate and  GFP expression in UMAP filtered 
###Calculate cluster frequency for functional Markers (dtTomato, GFP)
#Back to expression data of cells
expr
class(expr)
class(expr$cell_id)

#Delete all no-cytokine data + cell_ID and save as Cyto_expr
Cyto_expr <- expr
head(Cyto_expr)
Cyto_expr[,c("CD25", "CD8a","CD44", "CD4", "sample_id", "CD27","CD28", "TCRb")] <- NULL
head(Cyto_expr)

 # #Combine cyto_express info of cells with group info (condition, genotype, stimulation)
## Project defined clustering on UMAP
data_umap_compl <- merge(umap_exp, data_plot_u[,c(lineage_markers, "cell_id")], by = "cell_id")
data_umap_compl <- dplyr::filter(data_umap_compl, cell_clustering1m != "NA")
head(data_umap_compl)
head(expr)
mm <- match(data_umap_compl$sample_id, md$sample_id)

head(data_umap_compl)
data_umap_compl_expr <- merge(Cyto_expr, data_umap_compl, by = "cell_id")
head(data_umap_compl_expr)


#########dtTomato
#Define and filter dtTomato cells (check FACS Blots for dtTomato generated earlier and set threshold, e.g. 0.25)
dtTomato <- data_umap_compl_expr %>%
  dplyr::filter(dtTomato >= 0.15, GFP <= 0.25)
head(dtTomato)
mm <- match(dtTomato$sample_id, md$sample_id)
head(dtTomato)
tail(dtTomato)
dtTomato$group <- factor(md$group[mm])

#Create Count of dtTomato+ cells by making a table out of cell clusters and dtTomato+ samples
dtTomato_count <- table(dtTomato$cell_clustering1m, dtTomato$sample_id)
head(dtTomato_count)
countdtTomato <- as.data.frame.matrix(dtTomato_count, header = TRUE)

#total count without filter

mm <- match(data_umap_compl_expr$sample_id, md$sample_id)
head(data_umap_compl_expr)
tail(data_umap_compl_expr)
data_umap_compl_expr$group <- factor(md$group[mm])

#Create Count of dtTomato+ cells by making a table out of cell clusters and dtTomato+ samples
data_umap_compl_expr_count <- table(data_umap_compl_expr$cell_clustering1m, data_umap_compl_expr$sample_id)
head(data_umap_compl_expr)
countdata_umap_compl_expr <- as.data.frame.matrix(data_umap_compl_expr_count, header = TRUE)
head(countdata_umap_compl_expr)

#Calculate cluster frequency

ggdf_total <- melt(data.frame(cluster = rownames(countdata_umap_compl_expr), countdata_umap_compl_expr), id.vars = "cluster", value.name = "Frequency_total", variable.name = "sample_id")
ggdf_total

ggdf_dtTomato <- melt(data.frame(cluster = rownames(countdtTomato), countdtTomato), id.vars = "cluster", value.name = "Frequency_var", variable.name = "sample_id")
ggdf_dtTomato

ggdf2_dtTomato <- merge(ggdf_total, ggdf_dtTomato, by = c("sample_id", "cluster"))
head(ggdf2_dtTomato)

ggdf2_dtTomato$Frequency <- (ggdf2_dtTomato$Frequency_var / ggdf2_dtTomato$Frequency_total)*100    # Creating new variable
ggdf2_dtTomato 

ggdf2_dtTomato[,c("Frequency_total", "Frequency_var")] <- NULL
head(ggdf2_dtTomato)

levels(ggdf2_dtTomato$sample_id) <- sub("^X", "", levels(ggdf2_dtTomato$sample_id))

ggdf2_dtTomato$sample_id

ggdf2_dtTomato$sample_id <- as.factor(ggdf2_dtTomato$sample_id)
ggdf2_dtTomato$sample_id <-as.numeric(paste(ggdf2_dtTomato$sample_id))

mm <- match(ggdf2_dtTomato$sample_id, md$sample_id)
ggdf2_dtTomato$group <- factor(md$group[mm])

##PLOT CELL CLUSTERS THAT SHOW dtTomato EXPRESSION
plotsubsamples_dtTomato <- ggplot(ggdf2_dtTomato) +
  geom_violin(aes(x = factor(group, level = c("Ctr", "LRH1 KO")), y = Frequency, color = group, fill = group), 
              alpha = 0.2, trim = TRUE, scale = "width") +
  geom_point(aes(x = factor(group, level = c("Ctr", "LRH1 KO")), y = Frequency, color = group),
             alpha = 1, size = 3, stroke = 1.0, 
             position = position_jitter(0.05)) + 
  stat_summary(aes(y=Frequency, x=group, color=group), fun.data = "mean_se", geom = "crossbar", width = 0.15,  colour = "black") + #stat_summary(aes(y=Counts, x=condition, color=condition),fun.data="mean_se", color = "black", size = 0.5) +
  ylab("dtTomato+ cells (% of total cells)") +
  ggtitle("dtTomato expression") +
  facet_wrap(~ cluster, scale = "free_y", nrow=2) + #facet_grid(~cluster) + scales ="free_y",
  theme_light() +
  theme(axis.text.x = element_text(size = 20, hjust = 1, vjust = 1, angle = 45, color = "black"), 
        axis.ticks.x = element_line(size = 1, color = "black"),
        axis.ticks.y = element_line(size = 1, color = "black"),
        axis.ticks.length.x = unit(.30, "cm"),
        axis.ticks.length.y = unit(.30, "cm"),
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y.left = element_text(size = 23, face = "plain", vjust = 2),
        strip.background = element_rect(colour=NA, size = 0.8, fill = "white"),
        panel.border = element_rect(colour = "black", size = 0.8),
        #strip.text = element_blank(), 
        strip.text = element_text(color = "black", size = 20, face="bold"),
        axis.line = element_line(size = 0.6, color = "black"),
        #panel.grid.major.x = element_blank(),
        legend.text = element_blank(), legend.title = element_blank(), legend.position = "none",
        plot.title = element_text(color="black", size=24, hjust = 0.5)) + 
  scale_y_continuous() +
  expand_limits(y=0) +
  scale_color_manual(values = color_group) + 
  scale_fill_manual(values = color_group) 
plotsubsamples_dtTomato

write.xlsx(ggdf2_dtTomato,"dtTomato_cluster.xlsx")

####Tota cell number
ggdf3_dtTomato <- merge(ggdf2_dtTomato, ggdf3, by = c("sample_id", "cluster", "group"))

head(ggdf3_dtTomato)
tail(ggdf3_dtTomato)

ggdf3_dtTomato$total_cell_dtTomato <- (ggdf3_dtTomato$total_cell_nr / 100) *ggdf3_dtTomato$Frequency   # Creating new variable
head(ggdf3_dtTomato) 

ggdf3_dtTomato[,c("Frequency", "cell_number")] <- NULL
head(ggdf3_dtTomato)
class(ggdf3_dtTomato)

ggdf3_dtTomato <- do.call(data.frame,                      # Replace Inf in data by NA
                          lapply(ggdf3_dtTomato,
                                 function(x) replace(x, is.infinite(x), NA)))

ggdf3_dtTomato$total_cell_dtTomato[is.na(ggdf3_dtTomato$total_cell_dtTomato)] <- 0

ggdf3_dtTomato

##PLOT CELL CLUSTERS THAT SHOW dtTomato EXPRESSION
plotsubsamples_dtTomato_cellCount <- ggplot(ggdf3_dtTomato) +
  geom_violin(aes(x = factor(group, level = c("Ctr", "LRH1 KO")), y = total_cell_dtTomato, color = group, fill = group), 
              alpha = 0.2, trim = TRUE, scale = "width") +
  geom_point(aes(x = factor(group, level = c("Ctr", "LRH1 KO")), y = total_cell_dtTomato, color = group),
             alpha = 1, size = 3, stroke = 1.0, 
             position = position_jitter(0.05)) + 
  stat_summary(aes(y=total_cell_dtTomato, x=group, color=group), fun.data = "mean_se", geom = "crossbar", width = 0.15,  colour = "black") + #stat_summary(aes(y=Counts, x=condition, color=condition),fun.data="mean_se", color = "black", size = 0.5) +
  ylab("dtTomato+ total cell number") +
  ggtitle("dtTomato expression") +
  facet_wrap(~ cluster, scale = "free_y", nrow=2) + #facet_grid(~cluster) + scales ="free_y",
  theme_light() +
  theme(axis.text.x = element_text(size = 20, hjust = 1, vjust = 1, angle = 45, color = "black"), 
        axis.ticks.x = element_line(size = 1, color = "black"),
        axis.ticks.y = element_line(size = 1, color = "black"),
        axis.ticks.length.x = unit(.30, "cm"),
        axis.ticks.length.y = unit(.30, "cm"),
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y.left = element_text(size = 23, face = "plain", vjust = 2),
        strip.background = element_rect(colour=NA, size = 0.8, fill = "white"),
        panel.border = element_rect(colour = "black", size = 0.8),
        #strip.text = element_blank(), 
        strip.text = element_text(color = "black", size = 20, face="bold"),
        axis.line = element_line(size = 0.6, color = "black"),
        #panel.grid.major.x = element_blank(),
        legend.text = element_blank(), legend.title = element_blank(), legend.position = "none",
        plot.title = element_text(color="black", size=24, hjust = 0.5)) + 
  scale_y_continuous() +
  expand_limits(y=0) +
  scale_color_manual(values = color_group) + 
  scale_fill_manual(values = color_group) 
plotsubsamples_dtTomato_cellCount

write.xlsx(ggdf3_dtTomato,"Cellcount_dtTomato_cluster.xlsx")

#########GFP
#Define and filter GFP cells (check FACS Blots for GFP generated earlier and set threshold, e.g. 0.25)
GFP <- data_umap_compl_expr %>%
  dplyr::filter(GFP >= 0.25)
head(GFP)
#Note: it may happen that some sample_ids have no GFP+ cells (check whether number of GFP+ sample_id = input sample_id number)
GFP$sample_id

#Combine (based on sample_id) the GFP+ cells with metadata (condition, genotype etc)
mm <- match(GFP$sample_id, md$sample_id)
head(GFP)
tail(GFP)
GFP$group <- factor(md$group[mm])
class(GFP)

#Create Count of GFP+ cells by making a table out of cell clusters and GFP+ samples
GFP_count <- table(GFP$cell_clustering1m, GFP$sample_id)
head(GFP_count)
count_GFP <- as.data.frame.matrix(GFP_count, header = TRUE)

#Calculate cluster frequency

ggdf_GFP <- melt(data.frame(cluster = rownames(count_GFP), count_GFP), id.vars = "cluster", value.name = "Frequency_var", variable.name = "sample_id")
ggdf_GFP


ggdf2_GFP <- merge(ggdf_total, ggdf_GFP, by = c("sample_id", "cluster"))
head(ggdf2_GFP)

ggdf2_GFP$Frequency <- (ggdf2_GFP$Frequency_var / ggdf2_GFP$Frequency_total)*100    # Creating new variable
ggdf2_GFP

ggdf2_GFP[,c("Frequency_total", "Frequency_var")] <- NULL
head(ggdf2_GFP)

levels(ggdf2_GFP$sample_id) <- sub("^X", "", levels(ggdf2_GFP$sample_id))

ggdf2_GFP$sample_id

ggdf2_GFP$sample_id <- as.factor(ggdf2_GFP$sample_id)
ggdf2_GFP$sample_id <-as.numeric(paste(ggdf2_GFP$sample_id))

mm <- match(ggdf2_GFP$sample_id, md$sample_id)
ggdf2_GFP$group <- factor(md$group[mm])
ggdf2_GFP

##PLOT CELL CLUSTERS THAT SHOW GFP EXPRESSION
plotsubsamples_GFP <- ggplot(ggdf2_GFP) +
  geom_violin(aes(x = factor(group, level = c("Ctr", "LRH1 KO")), y = Frequency, color = group, fill = group), 
              alpha = 0.2, trim = TRUE, scale = "width") +
  geom_point(aes(x = factor(group, level = c("Ctr", "LRH1 KO")), y = Frequency, color = group),
             alpha = 1, size = 3, stroke = 1.0, 
             position = position_jitter(0.05)) + 
  stat_summary(aes(y=Frequency, x=group, color=group), fun.data = "mean_se", geom = "crossbar", width = 0.15,  colour = "black") + #stat_summary(aes(y=Counts, x=condition, color=condition),fun.data="mean_se", color = "black", size = 0.5) +
  ylab("GFP+ cells (% of total cells)") +
  ggtitle("GFP expression") +
  facet_wrap(~ cluster, scale = "free_y", nrow=2) + #facet_grid(~cluster) + scales ="free_y",
  theme_light() +
  theme(axis.text.x = element_text(size = 20, hjust = 1, vjust = 1, angle = 45, color = "black"), 
        axis.ticks.x = element_line(size = 1, color = "black"),
        axis.ticks.y = element_line(size = 1, color = "black"),
        axis.ticks.length.x = unit(.30, "cm"),
        axis.ticks.length.y = unit(.30, "cm"),
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y.left = element_text(size = 23, face = "plain", vjust = 2),
        strip.background = element_rect(colour=NA, size = 0.8, fill = "white"),
        panel.border = element_rect(colour = "black", size = 0.8),
        #strip.text = element_blank(), 
        strip.text = element_text(color = "black", size = 20, face="bold"),
        axis.line = element_line(size = 0.6, color = "black"),
        #panel.grid.major.x = element_blank(),
        legend.text = element_blank(), legend.title = element_blank(), legend.position = "none",
        plot.title = element_text(color="black", size=24, hjust = 0.5)) + 
  scale_y_continuous() +
  expand_limits(y=0) +
  scale_color_manual(values = color_group) + 
  scale_fill_manual(values = color_group) 
plotsubsamples_GFP

write.xlsx(ggdf2_GFP,"GFP_cluster.xlsx")

####Tota cell number
ggdf3_GFP <- merge(ggdf2_GFP, ggdf3, by = c("sample_id", "cluster", "group"))

head(ggdf3_GFP)
tail(ggdf3_GFP)

ggdf3_GFP$total_cell_dtTomato <- (ggdf3_GFP$total_cell_nr / 100) *ggdf3_GFP$Frequency   # Creating new variable
head(ggdf3_GFP) 

ggdf3_GFP[,c("Frequency", "cell_number")] <- NULL
head(ggdf3_GFP)
class(ggdf3_GFP)

ggdf3_GFP<- do.call(data.frame,                      # Replace Inf in data by NA
                          lapply(ggdf3_GFP,
                                 function(x) replace(x, is.infinite(x), NA)))

ggdf3_GFP$total_cell_dtTomato[is.na(ggdf3_GFP$total_cell_GFP)] <- 0

ggdf3_GFP

##PLOT CELL CLUSTERS THAT SHOW dtTomato EXPRESSION
plotsubsamples_GFP_cellCount <- ggplot(ggdf3_GFP) +
  geom_violin(aes(x = factor(group, level = c("Ctr", "LRH1 KO")), y = total_cell_dtTomato, color = group, fill = group), 
              alpha = 0.2, trim = TRUE, scale = "width") +
  geom_point(aes(x = factor(group, level = c("Ctr", "LRH1 KO")), y = total_cell_dtTomato, color = group),
             alpha = 1, size = 3, stroke = 1.0, 
             position = position_jitter(0.05)) + 
  stat_summary(aes(y=total_cell_dtTomato, x=group, color=group), fun.data = "mean_se", geom = "crossbar", width = 0.15,  colour = "black") + #stat_summary(aes(y=Counts, x=condition, color=condition),fun.data="mean_se", color = "black", size = 0.5) +
  ylab("GFP+ total cell number") +
  ggtitle("GFP expression") +
  facet_wrap(~ cluster, scale = "free_y", nrow=2) + #facet_grid(~cluster) + scales ="free_y",
  theme_light() +
  theme(axis.text.x = element_text(size = 20, hjust = 1, vjust = 1, angle = 45, color = "black"), 
        axis.ticks.x = element_line(size = 1, color = "black"),
        axis.ticks.y = element_line(size = 1, color = "black"),
        axis.ticks.length.x = unit(.30, "cm"),
        axis.ticks.length.y = unit(.30, "cm"),
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y.left = element_text(size = 23, face = "plain", vjust = 2),
        strip.background = element_rect(colour=NA, size = 0.8, fill = "white"),
        panel.border = element_rect(colour = "black", size = 0.8),
        #strip.text = element_blank(), 
        strip.text = element_text(color = "black", size = 20, face="bold"),
        axis.line = element_line(size = 0.6, color = "black"),
        #panel.grid.major.x = element_blank(),
        legend.text = element_blank(), legend.title = element_blank(), legend.position = "none",
        plot.title = element_text(color="black", size=24, hjust = 0.5)) + 
  scale_y_continuous() +
  expand_limits(y=0) +
  scale_color_manual(values = color_group) + 
  scale_fill_manual(values = color_group) 
plotsubsamples_GFP_cellCount

write.xlsx(ggdf3_GFP,"Cellcount_GFP_cluster.xlsx")


### Expression of TCRb

### plot all marker vs each other in order to evaluate transformation quality ("facs blots")
# subsample cells for plotting

plot_surf_list <- list()
for(i in 1:length(colnames(data.shift))){
  plot_surf <- ggplot(s_surf, aes_string(x = colnames(s_surf)[i], y = "TCRb")) +
    geom_hex(bins = 100) +
    xlim(-0.1, 1.0) +
    ylim(-0.1, 1.0) +
    coord_fixed(ratio = 1) +
    scale_fill_gradientn(colours = color.scale, trans = "sqrt") +
    theme_bw() + theme(legend.position="none", axis.text = element_text(size = 5), axis.title = element_text(size = 8))
  plot_surf_list[[i]] <-  plot_surf
}

m1 <- arrangeGrob(grobs = plot_surf_list, ncol = 1, top=textGrob("TCRb", gp=gpar(fontsize=6)))
#for example

#########Plot TCRb in UMAP filtered 
###Calculate cluster frequency for functional Markers (dtTomato, GFP)
#Back to expression data of cells
expr
class(expr)
class(expr$cell_id)

#Delete all no-cytokine data + cell_ID and save as Cyto_expr
Cyto_expr <- expr
head(Cyto_expr)
Cyto_expr[,c("CD25", "CD8a","CD44", "CD4", "sample_id", "CD27","CD28","GFP", "dtTomato")] <- NULL
head(Cyto_expr)

# #Combine cyto_express info of cells with group info (condition, genotype, stimulation)
## Project defined clustering on UMAP
data_umap_compl <- merge(umap_exp, data_plot_u[,c(lineage_markers, "cell_id")], by = "cell_id")
data_umap_compl <- dplyr::filter(data_umap_compl, cell_clustering1m != "NA")
head(data_umap_compl)
head(expr)
mm <- match(data_umap_compl$sample_id, md$sample_id)

head(Cyto_expr)
head(data_umap_compl)
data_umap_compl_expr <- merge(Cyto_expr, data_umap_compl, by = "cell_id")
head(data_umap_compl_expr)


#########TCRb
#Define and filter TCRb cells (check FACS Blots for TCRb generated earlier and set threshold, e.g. 0.25)
TCRb <- data_umap_compl %>%
  dplyr::filter(TCRb >= 0.13)
head(TCRb)
mm <- match(TCRb$sample_id, md$sample_id)
head(TCRb)
tail(TCRb)
TCRb$group <- factor(md$group[mm])

#Create Count of TCRb+ cells by making a table out of cell clusters and TCRb+ samples
TCRb_count <- table(TCRb$cell_clustering1m, TCRb$sample_id)
head(TCRb_count)
countTCRb <- as.data.frame.matrix(TCRb_count, header = TRUE)

#total count without filter

mm <- match(data_umap_compl_expr$sample_id, md$sample_id)
head(data_umap_compl_expr)
tail(data_umap_compl_expr)
data_umap_compl_expr$group <- factor(md$group[mm])

#Create Count of TCRb+ cells by making a table out of cell clusters and TCRb+ samples
data_umap_compl_expr_count <- table(data_umap_compl_expr$cell_clustering1m, data_umap_compl_expr$sample_id)
head(data_umap_compl_expr)
countdata_umap_compl_expr <- as.data.frame.matrix(data_umap_compl_expr_count, header = TRUE)
head(countdata_umap_compl_expr)

#Calculate cluster frequency

ggdf_total <- melt(data.frame(cluster = rownames(countdata_umap_compl_expr), countdata_umap_compl_expr), id.vars = "cluster", value.name = "Frequency_total", variable.name = "sample_id")
ggdf_total

ggdf_TCRb <- melt(data.frame(cluster = rownames(countTCRb), countTCRb), id.vars = "cluster", value.name = "Frequency_var", variable.name = "sample_id")
ggdf_TCRb

ggdf2_TCRb <- merge(ggdf_total, ggdf_TCRb, by = c("sample_id", "cluster"))
head(ggdf2_TCRb)

ggdf2_TCRb$Frequency <- (ggdf2_TCRb$Frequency_var / ggdf2_TCRb$Frequency_total)*100    # Creating new variable
ggdf2_TCRb

ggdf2_TCRb[,c("Frequency_total", "Frequency_var")] <- NULL
head(ggdf2_TCRb)

levels(ggdf2_TCRb$sample_id) <- sub("^X", "", levels(ggdf2_TCRb$sample_id))

ggdf2_TCRb$sample_id

ggdf2_TCRb$sample_id <- as.factor(ggdf2_TCRb$sample_id)
ggdf2_TCRb$sample_id <-as.numeric(paste(ggdf2_TCRb$sample_id))

mm <- match(ggdf2_TCRb$sample_id, md$sample_id)
ggdf2_TCRb$group <- factor(md$group[mm])

##PLOT CELL CLUSTERS THAT SHOW dtTomato EXPRESSION
plotsubsamples_TCRb <- ggplot(ggdf2_TCRb) +
  geom_violin(aes(x = factor(group, level = c("Ctr", "LRH1 KO")), y = Frequency, color = group, fill = group), 
              alpha = 0.2, trim = TRUE, scale = "width") +
  geom_point(aes(x = factor(group, level = c("Ctr", "LRH1 KO")), y = Frequency, color = group),
             alpha = 1, size = 3, stroke = 1.0, 
             position = position_jitter(0.05)) + 
  stat_summary(aes(y=Frequency, x=group, color=group), fun.data = "mean_se", geom = "crossbar", width = 0.15,  colour = "black") + #stat_summary(aes(y=Counts, x=condition, color=condition),fun.data="mean_se", color = "black", size = 0.5) +
  ylab("TCRb+ cells (% of total cells)") +
  ggtitle("TCRb expression") +
  facet_wrap(~ cluster, scale = "free_y", nrow=2) + #facet_grid(~cluster) + scales ="free_y",
  theme_light() +
  theme(axis.text.x = element_text(size = 20, hjust = 1, vjust = 1, angle = 45, color = "black"), 
        axis.ticks.x = element_line(size = 1, color = "black"),
        axis.ticks.y = element_line(size = 1, color = "black"),
        axis.ticks.length.x = unit(.30, "cm"),
        axis.ticks.length.y = unit(.30, "cm"),
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y.left = element_text(size = 23, face = "plain", vjust = 2),
        strip.background = element_rect(colour=NA, size = 0.8, fill = "white"),
        panel.border = element_rect(colour = "black", size = 0.8),
        #strip.text = element_blank(), 
        strip.text = element_text(color = "black", size = 20, face="bold"),
        axis.line = element_line(size = 0.6, color = "black"),
        #panel.grid.major.x = element_blank(),
        legend.text = element_blank(), legend.title = element_blank(), legend.position = "none",
        plot.title = element_text(color="black", size=24, hjust = 0.5)) + 
  scale_y_continuous() +
  expand_limits(y=0) +
  scale_color_manual(values = color_group) + 
  scale_fill_manual(values = color_group) 
plotsubsamples_TCRb

write.xlsx(ggdf2_TCRb,"TCRb_cluster.xlsx")

####Tota cell number
ggdf3_TCRb <- merge(ggdf2_TCRb, ggdf3, by = c("sample_id", "cluster", "group"))

head(ggdf3_TCRb)
tail(ggdf3_TCRb)

ggdf3_TCRb$total_cell_TCRb <- (ggdf3_TCRb$total_cell_nr / 100) *ggdf3_TCRb$Frequency   # Creating new variable
head(ggdf3_TCRb) 

ggdf3_TCRb[,c("Frequency", "cell_number")] <- NULL
head(ggdf3_TCRb)
class(ggdf3_TCRb)

ggdf3_TCRb <- do.call(data.frame,                      # Replace Inf in data by NA
                          lapply(ggdf3_TCRb,
                                 function(x) replace(x, is.infinite(x), NA)))

ggdf3_TCRb$total_cell_TCRb[is.na(ggdf3_TCRb$total_cell_TCRb)] <- 0

ggdf3_TCRb

##PLOT CELL CLUSTERS THAT SHOW dtTomato EXPRESSION
plotsubsamples_TCRb_cellCount <- ggplot(ggdf3_TCRb) +
  geom_violin(aes(x = factor(group, level = c("Ctr", "LRH1 KO")), y = total_cell_TCRb, color = group, fill = group), 
              alpha = 0.2, trim = TRUE, scale = "width") +
  geom_point(aes(x = factor(group, level = c("Ctr", "LRH1 KO")), y = total_cell_TCRb, color = group),
             alpha = 1, size = 3, stroke = 1.0, 
             position = position_jitter(0.05)) + 
  stat_summary(aes(y=total_cell_TCRb, x=group, color=group), fun.data = "mean_se", geom = "crossbar", width = 0.15,  colour = "black") + #stat_summary(aes(y=Counts, x=condition, color=condition),fun.data="mean_se", color = "black", size = 0.5) +
  ylab("TCRb+ total cell number") +
  ggtitle("TCRb expression") +
  facet_wrap(~ cluster, scale = "free_y", nrow=2) + #facet_grid(~cluster) + scales ="free_y",
  theme_light() +
  theme(axis.text.x = element_text(size = 20, hjust = 1, vjust = 1, angle = 45, color = "black"), 
        axis.ticks.x = element_line(size = 1, color = "black"),
        axis.ticks.y = element_line(size = 1, color = "black"),
        axis.ticks.length.x = unit(.30, "cm"),
        axis.ticks.length.y = unit(.30, "cm"),
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y.left = element_text(size = 23, face = "plain", vjust = 2),
        strip.background = element_rect(colour=NA, size = 0.8, fill = "white"),
        panel.border = element_rect(colour = "black", size = 0.8),
        #strip.text = element_blank(), 
        strip.text = element_text(color = "black", size = 20, face="bold"),
        axis.line = element_line(size = 0.6, color = "black"),
        #panel.grid.major.x = element_blank(),
        legend.text = element_blank(), legend.title = element_blank(), legend.position = "none",
        plot.title = element_text(color="black", size=24, hjust = 0.5)) + 
  scale_y_continuous() +
  expand_limits(y=0) +
  scale_color_manual(values = color_group) + 
  scale_fill_manual(values = color_group) 
plotsubsamples_TCRb_cellCount

write.xlsx(ggdf3_TCRb,"Cellcount_TCRb_cluster.xlsx")

