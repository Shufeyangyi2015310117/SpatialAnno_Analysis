############################ theme ######################################
theme_Publication <- function(base_size=14, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.line = element_line(colour="black"),
            axis.ticks.x= element_blank(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.text = element_text(size=28),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.4, "cm"),
            legend.spacing  = unit(0, "cm"),
            legend.title = element_blank(),
            plot.margin=unit(c(10,5,5,0),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


#######################  Figure 5a  ######################################## 
## scSorter initial
## calculate Kappa f1s acc 
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
load("embryo1_allcelltype_unknown_8markers_SCINA_initial_icmem_lfc=0.1894_8marker_presomitic_8m_Gut.RData")
load("embryo1_allcelltype_unknown_8markers_scSorter_initial_icmem_lfc=_8marker_presomitic_8m_Gut_adjusted.RData")
load("embryo1_allcelltype_unknown_8markers_scSorter_8marker_presomitic_8m_Gut.RData")
load("embryo1_allcelltype_unknown_8markers_cellassign2_8marker_presomitic_8m_Gut.RData")
load("embryo1_allcelltype_unknown_8markers_drsc_measure_8marker_presomitic_8m_Gut.RData")
load("embryo1_allcelltype_unknown_8markers_Seurat_BayesSpace2_8marker_presomitic_8m_Gut.RData")

embryo1_cds = readRDS("embryo1_garnett.rds")
layer_Garnett = pData(embryo1_cds)[,5]

Pre_cell_type = matrix(0, length(y), 7)
Pre_cell_type[,1] = colnames(rho)[fit_s2$type]
Pre_cell_type[,2] = results$cell_labels
Pre_cell_type[,3] = rts$Pred_Type
Pre_cell_type[,4] = fit$cell_type
Pre_cell_type[,5] = layer_drsc
Pre_cell_type[,6] = layer_BayesSpace
Pre_cell_type[,7] = layer_Garnett

source("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/function_f1score.R")
library(psych)
value = matrix(0, 7, 5)
for (i in 1:7){
  idx1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
  idx2 = which(Pre_cell_type[,i]!="Unknown" & Pre_cell_type[,i]!="unknown" & Pre_cell_type[,i]!="other")
  idx3 = intersect(idx1,idx2)
  value[i,1] = cohen.kappa(x=cbind(as.character(y)[idx3], Pre_cell_type[idx3,i]))$kappa
  value[i,2] = mean(evaluate(as.character(y)[idx3], Pre_cell_type[idx3,i])$F1)
  value[i,3] = mean(as.character(y)[idx3] == Pre_cell_type[idx3,i])
  value[i,4] = mclust::adjustedRandIndex(as.character(y)[idx3], Pre_cell_type[idx3,i])
  value[i,5] = evaluate(as.character(y)[idx3], Pre_cell_type[idx3,i])$Acc
}
colnames(value) = c("Kappa", "mF1", "Acc", "ARI", "Acc")
rownames(value) = c("SpatialAnno", "SCINA", "scSorter", "CellAssign", "DR-SC", "BayesSpace", "Garnett")

value2 = value[c(1:4,7),1:3]

library(reshape2) 
df = melt(value2)
colnames(df) = c("method", "measurement", "score")
df$method = factor(df$method, levels = c("SpatialAnno", "scSorter", "SCINA", "Garnett", "CellAssign"))

library(ggsci)
library(systemfonts)
library(ggplot2)

# Basic barplot
p<-ggplot(data=df, aes(x= method, y=score, fill = method)) +
  geom_bar(stat="identity")+
  facet_grid(.~measurement) + 
  scale_fill_simpsons() + 
  theme_Publication(base_size=30) + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_all_celltypes_unknown_8marker_measure5_8marker_presomitic_4m_Gut.pdf"), plot = p, width = 25, height = 20, units = "cm")



#######################  Figure 5b  ########################################

index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y_character = as.character(y)
y_character[-index1]="Unknown"
dat = data.frame(pos[ ,1], pos[ ,2], factor(y_character, levels=sort(unique(y_character))))
names(dat)= c("imagerow", "imagecol", "cluster")

library(ggplot2)
p0 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=cluster)) +
  geom_point(size = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5)))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080"))

ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_ground_truth_without_legend.png"), plot = p0, width = 8.5, height = 10, units = "cm")





#######################  Figure 5b  ########################################

index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y_character = as.character(y)
y_character[-index1]="Unknown"
dat = data.frame(pos[ ,1], pos[ ,2], factor(y_character, levels=sort(unique(y_character))))
names(dat)= c("imagerow", "imagecol", "Cluster")

library(ggplot2)
p0 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=Cluster)) +
  geom_point(size = 0.5) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5)))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080"))

ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_ground_truth.png"), plot = p0, width = 20, height = 10, units = "cm")
# Using the cowplot package
legend <- cowplot::get_legend(p0)
pdf(file = "/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_ground_truth_legend.pdf",width = 5, height = 4)
grid.newpage()
grid.draw(legend)
dev.off()



#######################  Figure 5b  ########################################

### spatialAnno with 8markers of Presomitic Gut from embryo1
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
library(Seurat)
### SCINA without legend
load("embryo1_allcelltype_unknown_8markers_SCINA_initial_icmem_lfc=0.1894_8marker_presomitic_8m_Gut.RData")

celltype = results$cell_labels
celltype_reorder_rename = sort(unique(celltype))
celltype_reorder_rename[celltype_reorder_rename=="unknown"] = "Unknown"
dat = data.frame(pos[,1], pos[,2], factor(celltype, levels=sort(unique(celltype)), labels=celltype_reorder_rename))
names(dat)= c("imagerow", "imagecol", "cluster")

index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y2 = y[index1]

idx_col = match(sort(unique(celltype)), sort(unique(y2)))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=cluster)) +
  geom_point(size = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 2))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080")[idx_col])
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_all_celltypes_unknown_8marker_SCINA_without_legend_8marker_presomitic_8m_Gut.png"), plot = p1, width = 8.5, height = 10, units = "cm")



#######################  Figure 5b  ########################################

### spatialAnno with 8markers of Presomitic from embryo1
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
library(Seurat)
### scSorter
load("embryo1_allcelltype_unknown_8markers_scSorter_8marker_presomitic_8m_Gut.RData")

celltype = rts$Pred_Type
celltype_reorder_rename = sort(unique(celltype))
celltype_reorder_rename[celltype_reorder_rename=="unknown"] = "Unknown"
dat = data.frame(pos[,1], pos[,2], factor(celltype, levels=sort(unique(celltype)), labels=celltype_reorder_rename))
names(dat)= c("imagerow", "imagecol", "cluster")

index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y2 = y[index1]

idx_col = match(sort(unique(celltype)), sort(unique(y2)))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=cluster)) +
  geom_point(size = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 2))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080")[idx_col])
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_all_celltypes_unknown_8marker_scSorter_without_legend_8marker_presomitic_8m_Gut.png"), plot = p1, width = 8.5, height = 10, units = "cm")




#######################  Figure 5b  ########################################

### spatialAnno with 8markers of Presomitic from embryo1
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
library(Seurat)
### cellassign
load("embryo1_allcelltype_unknown_8markers_cellassign2_8marker_presomitic_8m_Gut.RData")

celltype = fit$cell_type
celltype_reorder_rename = sort(unique(celltype))
celltype_reorder_rename[celltype_reorder_rename=="unknown"] = "Unknown"
dat = data.frame(pos[,1], pos[,2], factor(celltype, levels=sort(unique(celltype)), labels=celltype_reorder_rename))
names(dat)= c("imagerow", "imagecol", "cluster")

index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y2 = y[index1]

idx_col = match(sort(unique(celltype)), sort(unique(y2)))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=cluster)) +
  geom_point(size = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 2))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080")[idx_col])
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_all_celltypes_unknown_8marker_cellassign_without_legend_8marker_presomitic_8m_Gut.png"), plot = p1, width = 8.5, height = 10, units = "cm")





#######################  Figure 5b  ########################################

### spatialAnno with 8markers of Presomitic from embryo1
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
library(Seurat)
### Garnett
load("embryo1_allcelltype_unknown_8markers_SCINA_initial_icmem_lfc=0.1894_8marker_presomitic_8m_Gut.RData")

library(garnett)
embryo1_cds = readRDS("embryo1_garnett.rds")
layer_Garnett = pData(embryo1_cds)[,5]
celltype = layer_Garnett
celltype_reorder_rename = sort(unique(celltype))
celltype_reorder_rename[celltype_reorder_rename=="unknown"] = "Unknown"
dat = data.frame(pos[,1], pos[,2], factor(celltype, levels=sort(unique(celltype)), labels=celltype_reorder_rename))
names(dat)= c("imagerow", "imagecol", "cluster")

index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y2 = y[index1]

idx_col = match(sort(unique(celltype)), sort(unique(y2)))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=cluster)) +
  geom_point(size = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 2))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080")[idx_col])
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_all_celltypes_unknown_8marker_garnett_without_legend_8marker_presomitic_8m_Gut.png"), plot = p1, width = 8.5, height = 10, units = "cm")



#######################  Figure 5b  ########################################

### spatialAnno with 8markers of Presomitic from embryo1
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
library(Seurat)
### ICMEM
load("embryo1_allcelltype_unknown_8markers_scSorter_initial_icmem_lfc=_8marker_presomitic_8m_Gut_adjusted.RData")

celltype = colnames(rho)[fit_s2$type]
celltype_reorder_rename = sort(unique(celltype))
celltype_reorder_rename[celltype_reorder_rename=="unknown"] = "Unknown"
dat = data.frame(pos[,1], pos[,2], factor(celltype, levels=sort(unique(celltype)), labels=celltype_reorder_rename))
names(dat)= c("imagerow", "imagecol", "cluster")

index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y2 = y[index1]

idx_col = match(sort(unique(celltype)), sort(unique(y2)))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=cluster)) +
  geom_point(size = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 2))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080")[idx_col])
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_all_celltypes_unknown_8marker_scSorter_initial_icmem_without_legend_lfc=_8marker_presomitic_8m_Gut.png"), plot = p1, width = 8.5, height = 10, units = "cm")



#######################  Figure 4d  ########################################

### dpt_pseudotime for embryo1
library(RColorBrewer)
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
library(Seurat)
Emall <- readRDS('/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo/counts.Rds')
meta = readRDS("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo/metadata.Rds")
head(meta)
idx1 = which(meta$embryo == "embryo1")
### Take out sample 1
X <- t(as.matrix(Emall[,idx1])) 
pos <- cbind(meta$x_global[idx1], meta$y_global[idx1])
rownames(pos) <- rownames(meta)[idx1]

dpt_pseudotime = read.csv("embryo1_8markers_dpt_pseudotime_presomitic_8m_Gut.csv")
idx = match(dpt_pseudotime$X, rownames(pos))
pseudotime = dpt_pseudotime$dpt_pseudotime
pseudotime[pseudotime==Inf] = 1
dat = data.frame(pos[idx,1], pos[idx,2], pseudotime)
names(dat)= c("imagerow", "imagecol", "pseudotime")

index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y2 = y[index1]

idx_col = match(sort(unique(celltype)), sort(unique(y2)))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=pseudotime)) +
  geom_point(size = 0.1) +
  theme(legend.position = "right",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  scale_fill_gradientn(colours = brewer.pal(9,"YlOrRd"))
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_all_celltypes_unknown_8marker_dpt_pseudotime_presomitic_8m_Gut.png"), plot = p1, width = 11, height = 10, units = "cm")



  

#######################  Figure 5c  ########################################

rm(list = ls())
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo")
library(Seurat)
Emall <- readRDS('counts.Rds')
meta = readRDS("metadata.Rds")
head(meta)
idx1 = which(meta$embryo == "embryo1")

### Take out sample 3
y <- meta$celltype_mapped_refined[idx1]
X <- t(as.matrix(Emall[,idx1])) 
pos <- cbind(meta$x_global[idx1], meta$y_global[idx1])


library(SingleCellExperiment)
fit_pca = prcomp(X)
fit_y = fit_pca$x[,1:15]
# -------------------------------------------------
# make BayesSpace metadata used in BayesSpace
counts <- t(X[,])

## Make array coordinates - filled rectangle
cdata <- list()
cdata$row <- pos[,1]
cdata$col <- pos[,2]
cdata <- as.data.frame(do.call(cbind, cdata))
## Scale and jitter image coordinates
#scale.factor <- rnorm(1, 8);  n_spots <- n
#cdata$imagerow <- scale.factor * cdata$row + rnorm(n_spots)
#cdata$imagecol <- scale.factor * cdata$col + rnorm(n_spots)
cdata$imagerow <- cdata$row
cdata$imagecol <- cdata$col 
## Make SCE
## note: scater::runPCA throws warning on our small sim data, so use prcomp
sce <- SingleCellExperiment(assays=list(counts=counts), colData=cdata)
reducedDim(sce, "PCA") <- fit_y
sce$spatial.cluster <- floor(runif(ncol(sce), 1, 3))

metadata(sce)$BayesSpace.data <- list()
metadata(sce)$BayesSpace.data$platform <- "Visium"
metadata(sce)$BayesSpace.data$is.enhanced <- FALSE

set.seed(101)
library(scuttle)
sce <- logNormCounts(sce)
#dec <- scran::modelGeneVar(sce)
#top <- scran::getTopHVGs(dec, n = 2000)
#sce <- scater::runPCA(sce, subset_row = top)
X = logcounts(sce)

# $Endothelium
# [1] "Cdh5"   "Plvap"  "Eng"    "Cd34"   "Cldn5"  "Kdr"    "Sox18"  "Pecam1"
# $`Haematoendothelial progenitors`
# [1] "Cldn5" "Sox18" "Cdh5"  "Plvap" "Cd34"  "Eng"   "Kdr"   "Sox7"
# $`Spinal cord`
# [1] "Hoxb9" "Hoxd4" "Sox2"  "Hoxb8" "Foxb1" "Hoxc6" "Hoxc8" "Foxa2"
# $`Presomitic mesoderm`
# [1] "Tbx6"   "Dll1"   "Dll3"   "Pcdh19" "Lef1"   "Hoxb1"  "Hes7"   "Mesp2" 
# $Dermomyotome
# [1] "Meox1"   "Aldh1a2" "Six1"    "Hoxb3"   "Col26a1" "Foxc2"   "Hoxb4"  
# [8] "Fst" 
  
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2")
load("embryo1_allcelltype_unknown_8markers_SCINA_initial_icmem_lfc=0.1894_8marker_presomitic_8m_Gut.RData")
dim(pos)
celltype = sort(unique(colnames(rho)[fit_s2$type]))
genes = c("Meox1", "Cdh5","Cldn5", "Dll1", "Hoxb9")
for (i in 1:5){
  if (i == 1){
    idx = match(genes[i], rownames(X))
    dat = data.frame(pos[,1], pos[,2], X[idx,])
    names(dat)= c("imagerow", "imagecol", "cluster")
    dat$gene = genes[i]
  }else{
    idx = match(genes[i], rownames(X))
    dat2 = data.frame(pos[,1], pos[,2], X[idx,])
    names(dat2)= c("imagerow", "imagecol", "cluster")
    dat2$gene = genes[i]
    dat = rbind(dat, dat2)
  }
}
cols <- c("#0571B0",  "#CA0020")
quant = 0.5
med <- apply(X[genes,], 1, quantile, quant)
dat$gene = factor(dat$gene, levels = genes)
dat$marker = "marker" 
p1 <- ggplot(dat, aes(x=imagerow, y=3000 - imagecol, color=cluster))+
  geom_point(size = 0.1) + 
  facet_grid(marker~gene, scales = "free")+
  scale_colour_gradient2(
    low = cols[1],
    mid = "white",
    high = cols[2], midpoint = 0.5)+
  theme_Publication()+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_text(size = 8))

ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_allcelltype_unknown_8marker_celltype_specific_marker.pdf"), plot = p1, width = 27.6, height = 7, units = "cm")



######### marker gene plot with legend #########################
rm(list = ls())
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo")
library(Seurat)
Emall <- readRDS('counts.Rds')
meta = readRDS("metadata.Rds")
head(meta)
idx1 = which(meta$embryo == "embryo1")

### Take out sample 3
y <- meta$celltype_mapped_refined[idx1]
X <- t(as.matrix(Emall[,idx1])) 
pos <- cbind(meta$x_global[idx1], meta$y_global[idx1])


library(SingleCellExperiment)
fit_pca = prcomp(X)
fit_y = fit_pca$x[,1:15]
# -------------------------------------------------
# make BayesSpace metadata used in BayesSpace
counts <- t(X[,])

## Make array coordinates - filled rectangle
cdata <- list()
cdata$row <- pos[,1]
cdata$col <- pos[,2]
cdata <- as.data.frame(do.call(cbind, cdata))
## Scale and jitter image coordinates
#scale.factor <- rnorm(1, 8);  n_spots <- n
#cdata$imagerow <- scale.factor * cdata$row + rnorm(n_spots)
#cdata$imagecol <- scale.factor * cdata$col + rnorm(n_spots)
cdata$imagerow <- cdata$row
cdata$imagecol <- cdata$col 
## Make SCE
## note: scater::runPCA throws warning on our small sim data, so use prcomp
sce <- SingleCellExperiment(assays=list(counts=counts), colData=cdata)
reducedDim(sce, "PCA") <- fit_y
sce$spatial.cluster <- floor(runif(ncol(sce), 1, 3))

metadata(sce)$BayesSpace.data <- list()
metadata(sce)$BayesSpace.data$platform <- "Visium"
metadata(sce)$BayesSpace.data$is.enhanced <- FALSE

set.seed(101)
library(scuttle)
sce <- logNormCounts(sce)
#dec <- scran::modelGeneVar(sce)
#top <- scran::getTopHVGs(dec, n = 2000)
#sce <- scater::runPCA(sce, subset_row = top)
X = logcounts(sce)

# $Endothelium
# [1] "Cdh5"   "Plvap"  "Eng"    "Cd34"   "Cldn5"  "Kdr"    "Sox18"  "Pecam1"
# $`Haematoendothelial progenitors`
# [1] "Cldn5" "Sox18" "Cdh5"  "Plvap" "Cd34"  "Eng"   "Kdr"   "Sox7"
# $`Spinal cord`
# [1] "Hoxb9" "Hoxd4" "Sox2"  "Hoxb8" "Foxb1" "Hoxc6" "Hoxc8" "Foxa2"
# $`Presomitic mesoderm`
# [1] "Tbx6"   "Dll1"   "Dll3"   "Pcdh19" "Lef1"   "Hoxb1"  "Hes7"   "Mesp2" 
# $Dermomyotome
# [1] "Meox1"   "Aldh1a2" "Six1"    "Hoxb3"   "Col26a1" "Foxc2"   "Hoxb4"  
# [8] "Fst" 

setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2")
load("embryo1_allcelltype_unknown_8markers_SCINA_initial_icmem_lfc=0.1894_8marker_presomitic_8m_Gut.RData")
dim(pos)
celltype = sort(unique(colnames(rho)[fit_s2$type]))
genes = c("Meox1", "Cdh5","Cldn5", "Dll1", "Hoxb9")
for (i in 1:5){
  if (i == 1){
    idx = match(genes[i], rownames(X))
    dat = data.frame(pos[,1], pos[,2], X[idx,])
    names(dat)= c("imagerow", "imagecol", "cluster")
    dat$gene = genes[i]
  }else{
    idx = match(genes[i], rownames(X))
    dat2 = data.frame(pos[,1], pos[,2], X[idx,])
    names(dat2)= c("imagerow", "imagecol", "cluster")
    dat2$gene = genes[i]
    dat = rbind(dat, dat2)
  }
}
cols <- c("#0571B0",  "#CA0020")
quant = 0.5
med <- apply(X[genes,], 1, quantile, quant)
dat$gene = factor(dat$gene, levels = genes)
dat$marker = "marker" 
p1 <- ggplot(dat, aes(x=imagerow, y=3000 - imagecol, color=cluster))+
  geom_point(size = 0.1) + 
  facet_grid(marker~gene, scales = "free")+
  scale_colour_gradient2(
    low = cols[1],
    mid = "white",
    high = cols[2], midpoint = 0.5)+
  theme_Publication()+
  theme(legend.text=element_text(size=5),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_text(size = 8))

ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_allcelltype_unknown_8marker_celltype_specific_marker2.pdf"), plot = p1, width = 27.6, height = 7, units = "cm")



# Using the cowplot package
legend <- cowplot::get_legend(p1)
pdf(file = "/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_allcelltype_unknown_8marker_celltype_specific_marker_legend.pdf",width = 2, height = 1)
grid.newpage()
grid.draw(legend)
dev.off()


#######################  Figure 5c  ########################################

###### 7 celltype #########################################################
##### compare the performance of Ground truth and SpatialAnno cell type by cell type ##############
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2")
load("embryo1_allcelltype_unknown_8markers_scSorter_initial_icmem_lfc=_8marker_presomitic_8m_Gut_adjusted.RData")
method = c("Groundtruth", "SpatialAnno")
dim(pos)
library(scales)
library(RColorBrewer)
colfunc <- colorRampPalette(c("red", "white"))
library(ggplot2)
#celltype = sort(unique(colnames(rho)[fit_s2$type]))
data = list()
for (k in 1:2){
  if (k == 1){
    fit_s2_type = y
    unique_celltype_sort = c("Endothelium",
                             "Splanchnic mesoderm", 
                             "Gut tube",
                             "Presomitic mesoderm",
                             "Neural crest",
                             "Cardiomyocytes",
                             "Erythroid")
  }else{
    fit_s2_type = colnames(rho)[fit_s2$type]
    unique_celltype_sort = c("Endothelium",
                             "Splanchnic mesoderm",
                             "Gut tube",
                             "Presomitic mesoderm",
                             "Neural crest",
                             "Cardiomyocytes",
                             "Erythroid")
  }
  
  plist =list()
  for (i in 1:7){
    if (i == 1){
      fit_s2_type_one = rep("Unknown", length(fit_s2_type))
      fit_s2_type_one[fit_s2_type==unique_celltype_sort[i]] = unique_celltype_sort[i]
      dat = data.frame(pos[,1], pos[,2], factor(fit_s2_type_one))
      names(dat)= c("imagerow", "imagecol", "cluster")
      dat$celltype = unique_celltype_sort[i]
    }else{
      fit_s2_type_one = rep("Unknown", length(fit_s2_type))
      fit_s2_type_one[fit_s2_type==unique_celltype_sort[i]] = unique_celltype_sort[i]
      dat2 = data.frame(pos[,1], pos[,2], factor(fit_s2_type_one))
      names(dat2)= c("imagerow", "imagecol", "cluster")
      dat2$celltype = unique_celltype_sort[i]
      dat = rbind(dat, dat2)
    }
  }
  data[[k]] = dat
}

df1 = data[[1]]
df1$method = method[1]
df2 = data[[2]]
df2$method = method[2]
df = rbind(df1, df2)
df$method = factor(df$method, levels = c("Groundtruth", "SpatialAnno"))
df2 = df[df$cluster!="Unknown",]

p1 <- ggplot(df2, aes(x=imagerow, y=3000-imagecol, color=cluster))+
  geom_point(size = 0.1, alpha = 1) +
  facet_grid(method~celltype)+
  theme_Publication()+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_text(size = 8))
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_allcelltype_unknown_8markers_celltype_separately5_groundtruth_SpatialAnno.pdf"), plot = p1, width = 27.6, height = 10, units = "cm")



#######################  Figure 5c  ########################################
######### marker gene plot with legend #########################
#rm(list = ls())
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo")
library(Seurat)
Emall <- readRDS('counts.Rds')
meta = readRDS("metadata.Rds")
head(meta)
idx1 = which(meta$embryo == "embryo1")

### Take out sample 3
y <- meta$celltype_mapped_refined[idx1]
X <- t(as.matrix(Emall[,idx1])) 
pos <- cbind(meta$x_global[idx1], meta$y_global[idx1])


library(SingleCellExperiment)
fit_pca = prcomp(X)
fit_y = fit_pca$x[,1:15]
# -------------------------------------------------
# make BayesSpace metadata used in BayesSpace
counts <- t(X[,])

## Make array coordinates - filled rectangle
cdata <- list()
cdata$row <- pos[,1]
cdata$col <- pos[,2]
cdata <- as.data.frame(do.call(cbind, cdata))
## Scale and jitter image coordinates
#scale.factor <- rnorm(1, 8);  n_spots <- n
#cdata$imagerow <- scale.factor * cdata$row + rnorm(n_spots)
#cdata$imagecol <- scale.factor * cdata$col + rnorm(n_spots)
cdata$imagerow <- cdata$row
cdata$imagecol <- cdata$col 
## Make SCE
## note: scater::runPCA throws warning on our small sim data, so use prcomp
sce <- SingleCellExperiment(assays=list(counts=counts), colData=cdata)
reducedDim(sce, "PCA") <- fit_y
sce$spatial.cluster <- floor(runif(ncol(sce), 1, 3))

metadata(sce)$BayesSpace.data <- list()
metadata(sce)$BayesSpace.data$platform <- "Visium"
metadata(sce)$BayesSpace.data$is.enhanced <- FALSE

set.seed(101)
library(scuttle)
sce <- logNormCounts(sce)
#dec <- scran::modelGeneVar(sce)
#top <- scran::getTopHVGs(dec, n = 2000)
#sce <- scater::runPCA(sce, subset_row = top)
X = logcounts(sce)

# $Endothelium
# [1] "Cdh5"   "Plvap"  "Eng"    "Cd34"   "Cldn5"  "Kdr"    "Sox18"  "Pecam1"
# $`Splanchnic mesoderm`
# [1] "Foxf1" "Osr1"  "Hoxb1" "Isl1"  "Gata5" "Tbx5"  "Kcng1" "Gata4"
# $`Gut tube`
# [1] "Krt18" "Foxa1" "Cldn4" "Shh"   "Cdh1"  "Myh9"  "Cpm"   "Itga3"
# $`Presomitic mesoderm`
# [1] "Meox1"  "Dll3"   "Foxc2"  "Dll1"   "Lef1"   "Cer1"   "Notch1" "Mesp2" 
# $`Neural crest`
# [1] "Sox10"  "Tfap2b" "Tfap2a" "Nr2f1"  "Prrx1"  "Snai1"  "Msx1"   "Alx1"  
# $Cardiomyocytes
# [1] "Popdc2"  "Atp1b1"  "Tagln"   "Ttn"     "Smarcd3" "Gata5"   "Tbx5"   
# [8] "Hcn4"
# $Erythroid
# [1] "Slc4a1" "Alas2"  "Klf1"   "Gata1"  "Epor"   "Acp5"   "Hemgn"  "Smim1"


setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2")
load("embryo1_allcelltype_unknown_8markers_SCINA_initial_icmem_lfc=0.1894_8marker_presomitic_8m_Gut.RData")
dim(pos)
celltype = sort(unique(colnames(rho)[fit_s2$type]))
genes = c( "Popdc2","Cdh5", "Slc4a1", "Foxa1","Sox10", "Meox1", "Foxf1")
for (i in 1:7){
  if (i == 1){
    idx = match(genes[i], rownames(X))
    dat = data.frame(pos[,1], pos[,2], X[idx,])
    names(dat)= c("imagerow", "imagecol", "cluster")
    dat$gene = genes[i]
  }else{
    idx = match(genes[i], rownames(X))
    dat2 = data.frame(pos[,1], pos[,2], X[idx,])
    names(dat2)= c("imagerow", "imagecol", "cluster")
    dat2$gene = genes[i]
    dat = rbind(dat, dat2)
  }
}
cols <- c("#0571B0",  "#CA0020")
quant = 0.5
med <- apply(X[genes,], 1, quantile, quant)
dat$gene = factor(dat$gene, levels = genes)
dat$marker = "marker" 
p1 <- ggplot(dat, aes(x=imagerow, y=3000 - imagecol, color=cluster))+
  geom_point(size = 0.1) + 
  facet_grid(marker~gene, scales = "free")+
  scale_colour_gradient2(
    low = cols[1],
    mid = "white",
    high = cols[2], midpoint = 0.5)+
  theme_Publication()+
  theme(#legend.position="none",
        legend.text=element_text(size=5),
        legend.title= element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.x = element_text(size = 8,face = "italic"),
        strip.text.y = element_text(size = 8))

ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_allcelltype_unknown_8marker_celltype_specific_marker4.pdf"), plot = p1, width = 27.6, height = 6, units = "cm")


# Using the cowplot package
legend <- cowplot::get_legend(p1)
pdf(file = "/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_allcelltype_unknown_8marker_celltype_specific_marker4_legend.pdf",width = 2, height = 1)
grid.newpage()
grid.draw(legend)
dev.off()



#################### Supplementary Figure 26d  ############################

#############################################################
## boxplot of two genes Otx2 and sfrp1
#############################################################
rm(list = ls())
load("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/emberyo1_brain_slingtime_basedon_drsc_pen.const=1.53_seed100_2_marker_embedding_from_DRSC_scSorter_initial.RData")
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo")
library(Seurat)
gene = c("Otx2", "Sfrp1")
Emall <- readRDS('counts.Rds')
meta = readRDS("metadata.Rds")
head(meta)
idx1 = which(meta$embryo == "embryo1")

### Take out sample 3
y <- meta$celltype_mapped_refined[idx1]
X <- t(as.matrix(Emall[,idx1])) 
pos <- cbind(meta$x_global[idx1], meta$y_global[idx1])


library(SingleCellExperiment)
fit_pca = prcomp(X)
fit_y = fit_pca$x[,1:15]
# -------------------------------------------------
# make BayesSpace metadata used in BayesSpace
counts <- t(X[,])

## Make array coordinates - filled rectangle
cdata <- list()
cdata$row <- pos[,1]
cdata$col <- pos[,2]
cdata <- as.data.frame(do.call(cbind, cdata))
## Scale and jitter image coordinates
#scale.factor <- rnorm(1, 8);  n_spots <- n
#cdata$imagerow <- scale.factor * cdata$row + rnorm(n_spots)
#cdata$imagecol <- scale.factor * cdata$col + rnorm(n_spots)
cdata$imagerow <- cdata$row
cdata$imagecol <- cdata$col 
## Make SCE
## note: scater::runPCA throws warning on our small sim data, so use prcomp
sce <- SingleCellExperiment(assays=list(counts=counts), colData=cdata)
reducedDim(sce, "PCA") <- fit_y
sce$spatial.cluster <- floor(runif(ncol(sce), 1, 3))

metadata(sce)$BayesSpace.data <- list()
metadata(sce)$BayesSpace.data$platform <- "Visium"
metadata(sce)$BayesSpace.data$is.enhanced <- FALSE

set.seed(101)
library(scuttle)
sce <- logNormCounts(sce)
#dec <- scran::modelGeneVar(sce)
#top <- scran::getTopHVGs(dec, n = 2000)
#sce <- scater::runPCA(sce, subset_row = top)
X = t(logcounts(sce))



idx = which(y == "Forebrain/Midbrain/Hindbrain")

df = as.data.frame(pos[idx,])
# Rhombencephalon 后脑
# Tegmentum 中脑背盖部
# Mesencephalon 中脑
# Prosencephalon 前脑
singcluster2 = "Rhombencephalon(Glia2)"
singcluster2[singcluster==1] = "Tegmentum(Glia3)"
singcluster2[singcluster==2] = "Prosencephalon(Neurons)"
singcluster2[singcluster==3] = "Mesencephalon(Glia1)"
singcluster2[singcluster==4] = "Rhombencephalon(Glia2)"


rownames(pos) = rownames(X)
df = as.data.frame(pos[idx,])
df$cluster = singcluster
df$y = 1600 - df$V2
df2 = df[df$y>175 & df$V1 < 1640, ]
df3 = df2[df2$cluster=="3" | df2$cluster=="4" | df2$cluster=="2" ,]
df3[df3$cluster=="3" | df3$cluster=="2", ] = "Mesencephalon+\nProsencephalon"
df3[df3$cluster=="4",] = "Rhombencephalon"

p1 <- ggplot(df2, aes(x=V1, y=y, color=cluster))+
  geom_point(size = 0.1)

idx4 = match(rownames(df3), rownames(X))
data = as.data.frame(X[idx4, gene])
data$cluster = df3$cluster
data2 = melt(data)

colnames(data2) = c("cluster", "gene", "logcount")
p<-ggplot(data2, aes(x=cluster, y=logcount, color=cluster)) +
  geom_boxplot() + 
  facet_grid(.~gene)+
  theme(strip.text.x = element_text(size = 8,face = "italic"))
p
ggsave(p, file = "/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_brain_geneexp_rmoutlier_basedon_drsc_boxplot_pen.const=1.53_seed100_2_scSorter_initial_initial_cluster3.pdf", width = 8, height = 3)



####################  Figure 5d  ############################
#############################################################
#### initial cluster3 branch1 marker embedding from DRSC
load("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/emberyo1_brain_slingtime_basedon_drsc_pen.const=1.53_seed100_2_marker_embedding_from_DRSC_scSorter_initial_initial_cluster3.RData")
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo")
library(Seurat)
Emall <- readRDS('counts.Rds')
meta = readRDS("metadata.Rds")
head(meta)
idx1 = which(meta$embryo == "embryo1")


### Take out sample 3
y <- meta$celltype_mapped_refined[idx1]
X <- t(as.matrix(Emall[,idx1])) 
pos <- cbind(meta$x_global[idx1], meta$y_global[idx1])

idx = which(y == "Forebrain/Midbrain/Hindbrain")


df = as.data.frame(pos[idx,])
df$pseudotime = colData(sce_cancer)[,11]
df$y = 1600 - df$V2
df = df[df$y>175 & df$V1 < 1640, ]
df$branch = "branch1"

df2 = as.data.frame(pos[idx,])
df2$pseudotime = colData(sce_cancer)[,12]
df2$y = 1600 - df2$V2
df2 = df2[df2$y>175 & df2$V1 < 1640, ]
df2$branch = "branch2"

df = rbind(df, df2)

p <- ggplot(df, aes(V1, y, colour=pseudotime)) + 
  geom_point(size = 0.5)+
  facet_grid(.~branch)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
p
ggsave(p, file = paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_brain_pseudotime_rmoutlier_basedon_drsc_pen.const=1.53_seed100_2_marker_embedding_from_DRSC_scSorter_initial_initial_cluster3.pdf"), width = 7, height = 3)




#################### Supplementary Figure 26a  ############################
#############################################################
## Perform slingshot with initial cluster 3 on embeddings obtained from PCA
rm(list = ls())
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
load("emberyo1_brain_PCA_slingtime_startclus_basedon_drsc_pen.const=1.53_seed100_2_initial_cluster3.RData")
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo")
library(Seurat)
Emall <- readRDS('counts.Rds')
meta = readRDS("metadata.Rds")
head(meta)
idx1 = which(meta$embryo == "embryo1")

### Take out sample 3
y <- meta$celltype_mapped_refined[idx1]
X <- t(as.matrix(Emall[,idx1])) 
pos <- cbind(meta$x_global[idx1], meta$y_global[idx1])

idx = which(y == "Forebrain/Midbrain/Hindbrain")

df = as.data.frame(pos[idx,])
df$Pseudotime = colData(sce_cancer)[,12]
df$branch = "branch1"
df$y = 1600 - df$V2

df_2 = as.data.frame(pos[idx,])
df_2$Pseudotime = colData(sce_cancer)[,13]
df_2$branch = "branch2"
df_2$y = 1600 - df_2$V2

df_3 = as.data.frame(pos[idx,])
df_3$Pseudotime = colData(sce_cancer)[,14]
df_3$branch = "branch3"
df_3$y = 1600 - df_3$V2

df = rbind(df, df_2, df_3)

df2 = df[df$y>175 & df$V1 < 1640, ]
p <- ggplot(df2, aes(V1, y, colour=Pseudotime)) + 
  geom_point(size = 0.5)+
  facet_grid(.~branch)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
p
ggsave(p, file = "/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_brain_PCA_pseudotime_rmoutlier_startclus_basedon_drsc_multbranch_pen.const=1.53_seed100_2_initial_cluster3.pdf", width = 10, height = 3)




#################### Supplementary Figure 26b  ############################
#############################################################
## Perform slingshot with initial cluster 3 on embeddings obtained from DRSC
rm(list = ls())
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
load("emberyo1_brain_DR-SC_slingtime_startclus_basedon_drsc_pen.const=1.53_seed100_2_initial_cluster3.RData")
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo")
library(Seurat)
Emall <- readRDS('counts.Rds')
meta = readRDS("metadata.Rds")
head(meta)
idx1 = which(meta$embryo == "embryo1")

### Take out sample 3
y <- meta$celltype_mapped_refined[idx1]
X <- t(as.matrix(Emall[,idx1])) 
pos <- cbind(meta$x_global[idx1], meta$y_global[idx1])

idx = which(y == "Forebrain/Midbrain/Hindbrain")

df = as.data.frame(pos[idx,])
df$Pseudotime = colData(sce_cancer)[,12]
df$branch = "branch1"
df$y = 1600 - df$V2

df_2 = as.data.frame(pos[idx,])
df_2$Pseudotime = colData(sce_cancer)[,13]
df_2$branch = "branch2"
df_2$y = 1600 - df_2$V2

df_3 = as.data.frame(pos[idx,])
df_3$Pseudotime = colData(sce_cancer)[,14]
df_3$branch = "branch3"
df_3$y = 1600 - df_3$V2

df = rbind(df, df_2, df_3)

df2 = df[df$y>175 & df$V1 < 1640, ]
p <- ggplot(df2, aes(V1, y, colour=Pseudotime)) + 
  geom_point(size = 0.5)+
  facet_grid(.~branch)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
p
ggsave(p, file = "/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo1_brain_DR-SC_pseudotime_rmoutlier_startclus_basedon_drsc_multbranch_pen.const=1.53_seed100_2_initial_cluster3.pdf", width = 10, height = 3)



#################### Supplementary Figure 26c  ############################
### DE genes Changes along a trajectory scSorter initial
#############################################################
rm(list = ls())
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
load("emberyo1_brain_slingtime_basedon_drsc_pen.const=1.53_seed100_2_marker_embedding_from_DRSC_scSorter_initial_initial_cluster3.RData")
library(SingleCellExperiment)
sce_cancer
sce_cancer$Clusters = sce_cancer$slingClusters
sce_cancer$Clusters[sce_cancer$slingClusters=="1"] = "Tegmentum"
sce_cancer$Clusters[sce_cancer$slingClusters=="2"] = "Prosencephalon"
sce_cancer$Clusters[sce_cancer$slingClusters=="3"] = "Mesencephalon"
sce_cancer$Clusters[sce_cancer$slingClusters=="4"] = "Rhombencephalon"

library(scran)
library(TSCAN)
pseudotime = apply(colData(sce_cancer)[,c(11,12)], 1, mean, na.rm = TRUE)
pseudo <- testPseudotime(sce_cancer, pseudotime=pseudotime)


pseudo[order(pseudo$p.value),]
sorted <- pseudo[order(pseudo$p.value),]
up.left <- sorted # [sorted$logFC < 0,]
head(up.left, 20)

best <- head(row.names(up.left), 10)
rowData(sce_cancer)$SYMBOL <- row.names(sce_cancer)
sce_cancer$Pseudotime <- pseudotime

features_heat <- rownames(head(up.left, 20))
features_heat <- c("Hes3", "Irx3",  "Fgfr3","Sfrp2","Meis2","Cntfr",
                   "Nr2f1","Lfng","Fst","Irx5", "Irx2",
                   "Irx1", "Fgf15", "Otx2", "Pou3f1","Sfrp1",
                   "Lhx2","Fezf1", "Dusp6","Shh")

# features_heat <- c("Sfrp1", "Sfrp2", "Fgfr3","Meis2", "Irx3","Fst","Foxb1", 
#                    "Shh","Dusp6","Col4a1", "Foxa1", "En1", "Foxa2", 
#                    "Otx2","Irx1","Lhx2",
#                    "Shisa2", "Six3", 
#                    "Dlk1", "Fezf1")

library(scater)
p1 <- plotHeatmap(sce_cancer,  order_columns_by="Pseudotime",
                  color_columns_by="Clusters", features = features_heat,
                  labels_row = as.expression(lapply(features_heat, function(a) bquote(italic(.(a))))),
                  center=TRUE, swap_rownames="SYMBOL", cluster_rows=F, zlim = c(-3,3))
dir.file <- '/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/'
ggsave(paste0(dir.file, 'embryo1_brain_Traject_HeatMap_basedon_drsc_pen.const=1.53_seed100_2_marker_embedding_from_DRSC_scSorter_initial_initial_cluster3.pdf'), plot = p1, width = 12, height = 5, units = "in", dpi = 1000)









## embryo2
#################### Supplementary Figure 25a  ############################
## add garnett
library(garnett)
#rm(list = ls())
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
## calculate Kappa f1s acc 
load("embryo2_allcelltype_unknown_8markers_SCINA_initial_icmem_lfc=_8marker_presomitic_8m_Gut_adjusted.RData")
load("embryo2_allcelltype_unknown_8markers_scSorter_initial_icmem_lfc=_8marker_presomitic_8m_Gut_adjusted.RData")
load("embryo2_allcelltype_unknown_8markers_scSorter_8marker_presomitic_8m_Gut.RData")
load("embryo2_allcelltype_unknown_8markers_cellassign2_8marker_presomitic_8m_Gut.RData")

embryo2_cds = readRDS("embryo2_garnett.rds")
layer_Garnett = pData(embryo2_cds)[,5]

Pre_cell_type = matrix(0, length(y), 7)
Pre_cell_type[,1] = colnames(rho)[fit_s2$type]
Pre_cell_type[,2] = results$cell_labels
Pre_cell_type[,3] = rts$Pred_Type
Pre_cell_type[,4] = fit$cell_type
Pre_cell_type[,7] = layer_Garnett

source("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/function_f1score.R")
library(psych)
value = matrix(0, 7, 5)
for (i in 1:7){
  idx1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
  idx2 = which(Pre_cell_type[,i]!="Unknown" & Pre_cell_type[,i]!="unknown" & Pre_cell_type[,i]!="other")
  idx3 = intersect(idx1,idx2)
  value[i,1] = cohen.kappa(x=cbind(as.character(y)[idx3], Pre_cell_type[idx3,i]))$kappa
  value[i,2] = mean(evaluate(as.character(y)[idx3], Pre_cell_type[idx3,i])$F1)
  value[i,3] = mean(as.character(y)[idx3] == Pre_cell_type[idx3,i])
  value[i,4] = mclust::adjustedRandIndex(as.character(y)[idx3], Pre_cell_type[idx3,i])
  value[i,5] = evaluate(as.character(y)[idx3], Pre_cell_type[idx3,i])$Acc
}
colnames(value) = c("Kappa", "mF1", "Acc", "ARI", "acc")
rownames(value) = c("SpatialAnno", "SCINA", "scSorter", "CellAssign", "imputed2000_pruning", "imputed1000", "Garnett")

value2 = value[c(1:4,7),1:3]

library(reshape2) 
df = melt(value2)
colnames(df) = c("method", "measurement", "score")
df$method = factor(df$method, levels = c("SpatialAnno", "scSorter", "SCINA", "Garnett", "CellAssign"))

library(ggsci)
library(systemfonts)
library(ggplot2)

# Basic barplot
p<-ggplot(data=df, aes(x= method, y=score, fill = method)) +
  geom_bar(stat="identity")+
  facet_grid(.~measurement) + 
  scale_fill_simpsons() + 
  theme_Publication(base_size=30) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo2_all_celltypes_unknown_8marker_measure2_8marker_presomitic_8m_Gut.pdf"), plot = p, width = 25, height = 20, units = "cm")


#################### Supplementary Figure 25b  ############################
library(scales)
index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y_character = as.character(y)
y_character[-index1]="Unknown"
dat = data.frame(pos[ ,1], pos[ ,2], factor(y_character, levels=sort(unique(y_character))))
names(dat)= c("imagerow", "imagecol", "Cluster")

library(ggplot2)
p0 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=Cluster)) +
  geom_point(size = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5)))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080"))

ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo2_ground_truth_without_legend.png"), plot = p0, width = 8.5, height = 10, units = "cm")



#################### Supplementary Figure 25b  ############################
### spatialAnno with 8markers of Presomitic Gut from embryo1
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
library(Seurat)
### SCINA
load("embryo2_allcelltype_unknown_8markers_SCINA_initial_icmem_lfc=_8marker_presomitic_8m_Gut_adjusted.RData")

celltype = results$cell_labels
celltype_reorder_rename = sort(unique(celltype))
celltype_reorder_rename[celltype_reorder_rename=="unknown"] = "Unknown"
dat = data.frame(pos[,1], pos[,2], factor(celltype, levels=sort(unique(celltype)), labels=celltype_reorder_rename))
names(dat)= c("imagerow", "imagecol", "cluster")

index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y2 = y[index1]

idx_col = match(sort(unique(celltype)), sort(unique(y2)))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=cluster)) +
  geom_point(size = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 2))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080")[idx_col])
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo2_all_celltypes_unknown_8marker_SCINA_without_legend_8marker_presomitic_8m_Gut.png"), plot = p1, width = 8.5, height = 10, units = "cm")


#################### Supplementary Figure 25b  ############################
### spatialAnno with 8markers of Presomitic from embryo1
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
library(Seurat)
### ICMEM
load("embryo2_allcelltype_unknown_8markers_scSorter_initial_icmem_lfc=_8marker_presomitic_8m_Gut_adjusted.RData")

celltype = colnames(rho)[fit_s2$type]
celltype_reorder_rename = sort(unique(celltype))
celltype_reorder_rename[celltype_reorder_rename=="unknown"] = "Unknown"
dat = data.frame(pos[,1], pos[,2], factor(celltype, levels=sort(unique(celltype)), labels=celltype_reorder_rename))
names(dat)= c("imagerow", "imagecol", "cluster")

index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y2 = y[index1]

idx_col = match(sort(unique(celltype)), sort(unique(y2)))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=cluster)) +
  geom_point(size = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 2))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080")[idx_col])
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo2_all_celltypes_unknown_8marker_scSorter_initial_icmem_without_legend_lfc=_8marker_presomitic_8m_Gut.png"), plot = p1, width = 8.5, height = 10, units = "cm")



#################### Supplementary Figure 25b  ############################
### spatialAnno with 8markers of Presomitic from embryo1
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
library(Seurat)
### scSorter
load("embryo2_allcelltype_unknown_8markers_scSorter_8marker_presomitic_8m_Gut.RData")

celltype = rts$Pred_Type
celltype_reorder_rename = sort(unique(celltype))
celltype_reorder_rename[celltype_reorder_rename=="unknown"] = "Unknown"
dat = data.frame(pos[,1], pos[,2], factor(celltype, levels=sort(unique(celltype)), labels=celltype_reorder_rename))
names(dat)= c("imagerow", "imagecol", "cluster")

index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y2 = y[index1]

idx_col = match(sort(unique(celltype)), sort(unique(y2)))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=cluster)) +
  geom_point(size = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 2))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080")[idx_col])
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo2_all_celltypes_unknown_8marker_scSorter_without_legend_8marker_presomitic_8m_Gut.png"), plot = p1, width = 8.5, height = 10, units = "cm")


#################### Supplementary Figure 25b  ############################
### spatialAnno with 8markers of Presomitic from embryo1
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
library(Seurat)
### cellassign
load("embryo2_allcelltype_unknown_8markers_cellassign2_8marker_presomitic_8m_Gut.RData")

celltype = fit$cell_type
celltype[celltype=="other"] = "Unknown"
celltype_reorder_rename = sort(unique(celltype))

dat = data.frame(pos[,1], pos[,2], factor(celltype, levels=sort(unique(celltype)), labels=celltype_reorder_rename))
names(dat)= c("imagerow", "imagecol", "cluster")

index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y2 = y[index1]

idx_col = match(sort(unique(celltype)), sort(unique(y2)))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=cluster)) +
  geom_point(size = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 2))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080")[idx_col])
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo2_all_celltypes_unknown_8marker_cellassign_without_legend_8marker_presomitic_8m_Gut.png"), plot = p1, width = 8.5, height = 10, units = "cm")


#################### Supplementary Figure 25b  ############################
### spatialAnno with 8markers of Presomitic from embryo1
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
library(Seurat)
### Garnett
load("embryo2_allcelltype_unknown_8markers_scSorter_initial_icmem_lfc=_8marker_presomitic_8m_Gut_adjusted.RData")

library(garnett)
embryo1_cds = readRDS("embryo2_garnett.rds")
layer_Garnett = pData(embryo1_cds)[,5]
celltype = layer_Garnett
celltype_reorder_rename = sort(unique(celltype))
celltype_reorder_rename[celltype_reorder_rename=="unknown"] = "Unknown"
dat = data.frame(pos[,1], pos[,2], factor(celltype, levels=sort(unique(celltype)), labels=celltype_reorder_rename))
names(dat)= c("imagerow", "imagecol", "cluster")

index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y2 = y[index1]

idx_col = match(sort(unique(celltype)), sort(unique(y2)))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=cluster)) +
  geom_point(size = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 2))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080")[idx_col])
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo2_all_celltypes_unknown_8marker_garnett_without_legend_8marker_presomitic_8m_Gut.png"), plot = p1, width = 8.5, height = 10, units = "cm")





## embryo3
#################### Supplementary Figure 25c  ############################
## add garnett
rm(list = ls())
library(garnett)
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
## calculate Kappa f1s acc 
load("embryo3_allcelltype_unknown_8markers_SCINA_initial_icmem_lfc=_8marker_presomitic_8m_Gut_adjusted.RData")
load("embryo3_allcelltype_unknown_8markers_scSorter_initial_icmem_lfc=_8marker_presomitic_8m_Gut_adjusted.RData")
load("embryo3_allcelltype_unknown_8markers_scSorter_8marker_presomitic_8m_Gut.RData")
load("embryo3_allcelltype_unknown_8markers_cellassign2_8marker_presomitic_8m_Gut.RData")

embryo2_cds = readRDS("embryo3_garnett.rds")
layer_Garnett = pData(embryo2_cds)[,5]

Pre_cell_type = matrix(0, length(y), 7)
Pre_cell_type[,1] = colnames(rho)[fit_s2$type]
Pre_cell_type[,2] = results$cell_labels
Pre_cell_type[,3] = rts$Pred_Type
Pre_cell_type[,4] = fit$cell_type
Pre_cell_type[,7] = layer_Garnett

source("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/function_f1score.R")
library(psych)
value = matrix(0, 7, 5)
for (i in 1:7){
  idx1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
  idx2 = which(Pre_cell_type[,i]!="Unknown" & Pre_cell_type[,i]!="unknown" & Pre_cell_type[,i]!="other")
  idx3 = intersect(idx1,idx2)
  value[i,1] = cohen.kappa(x=cbind(as.character(y)[idx3], Pre_cell_type[idx3,i]))$kappa
  value[i,2] = mean(evaluate(as.character(y)[idx3], Pre_cell_type[idx3,i])$F1)
  value[i,3] = mean(as.character(y)[idx3] == Pre_cell_type[idx3,i])
  value[i,4] = mclust::adjustedRandIndex(as.character(y)[idx3], Pre_cell_type[idx3,i])
  value[i,5] = evaluate(as.character(y)[idx3], Pre_cell_type[idx3,i])$Acc
}
colnames(value) = c("Kappa", "mF1", "Acc", "ARI", "acc")
rownames(value) = c("SpatialAnno", "SCINA", "scSorter", "CellAssign", "imputed2000_pruning", "imputed1000", "Garnett")

value2 = value[c(1:4,7),1:3]

library(reshape2) 
df = melt(value2)
colnames(df) = c("method", "measurement", "score")
df$method = factor(df$method, levels = c("SpatialAnno", "scSorter", "SCINA", "Garnett", "CellAssign"))

library(ggsci)
library(systemfonts)
library(ggplot2)

# Basic barplot
p<-ggplot(data=df, aes(x= method, y=score, fill = method)) +
  geom_bar(stat="identity")+
  facet_grid(.~measurement) + 
  scale_fill_simpsons() + 
  theme_Publication(base_size=30) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo3_all_celltypes_unknown_8marker_measure2_8marker_presomitic_8m_Gut.pdf"), plot = p, width = 25, height = 20, units = "cm")



#################### Supplementary Figure 25d  ############################
library(scales)
index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y_character = as.character(y)
y_character[-index1]="Unknown"
dat = data.frame(pos[ ,1], pos[ ,2], factor(y_character, levels=sort(unique(y_character))))
names(dat)= c("imagerow", "imagecol", "Cluster")

library(ggplot2)
p0 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=Cluster)) +
  geom_point(size = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5)))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080"))

ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo3_ground_truth_without_legend.png"), plot = p0, width = 8.5, height = 10, units = "cm")


#################### Supplementary Figure 25d  ############################
### spatialAnno with 8markers of Presomitic Gut from embryo1
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
library(Seurat)
### SCINA
load("embryo3_allcelltype_unknown_8markers_SCINA_initial_icmem_lfc=_8marker_presomitic_8m_Gut_adjusted.RData")

celltype = results$cell_labels
celltype_reorder_rename = sort(unique(celltype))
celltype_reorder_rename[celltype_reorder_rename=="unknown"] = "Unknown"
dat = data.frame(pos[,1], pos[,2], factor(celltype, levels=sort(unique(celltype)), labels=celltype_reorder_rename))
names(dat)= c("imagerow", "imagecol", "cluster")

index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y2 = y[index1]

idx_col = match(sort(unique(celltype)), sort(unique(y2)))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=cluster)) +
  geom_point(size = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 2))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080")[idx_col])
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo3_all_celltypes_unknown_8marker_SCINA_without_legend_8marker_presomitic_8m_Gut.png"), plot = p1, width = 8.5, height = 10, units = "cm")



#################### Supplementary Figure 25d  ############################
### spatialAnno with 8markers of Presomitic from embryo1
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
library(Seurat)
### ICMEM
load("embryo3_allcelltype_unknown_8markers_scSorter_initial_icmem_lfc=_8marker_presomitic_8m_Gut_adjusted.RData")

celltype = colnames(rho)[fit_s2$type]
celltype_reorder_rename = sort(unique(celltype))
celltype_reorder_rename[celltype_reorder_rename=="unknown"] = "Unknown"
dat = data.frame(pos[,1], pos[,2], factor(celltype, levels=sort(unique(celltype)), labels=celltype_reorder_rename))
names(dat)= c("imagerow", "imagecol", "cluster")

index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y2 = y[index1]

idx_col = match(sort(unique(celltype)), sort(unique(y2)))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=cluster)) +
  geom_point(size = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 2))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080")[idx_col])
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo3_all_celltypes_unknown_8marker_scSorter_initial_icmem_without_legend_lfc=_8marker_presomitic_8m_Gut.png"), plot = p1, width = 8.5, height = 10, units = "cm")




#################### Supplementary Figure 25d  ############################
### spatialAnno with 8markers of Presomitic from embryo1
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
library(Seurat)
### scSorter
load("embryo3_allcelltype_unknown_8markers_scSorter_8marker_presomitic_8m_Gut.RData")

celltype = rts$Pred_Type
celltype_reorder_rename = sort(unique(celltype))
celltype_reorder_rename[celltype_reorder_rename=="unknown"] = "Unknown"
dat = data.frame(pos[,1], pos[,2], factor(celltype, levels=sort(unique(celltype)), labels=celltype_reorder_rename))
names(dat)= c("imagerow", "imagecol", "cluster")

index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y2 = y[index1]

idx_col = match(sort(unique(celltype)), sort(unique(y2)))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=cluster)) +
  geom_point(size = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 2))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080")[idx_col])
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo3_all_celltypes_unknown_8marker_scSorter_without_legend_8marker_presomitic_8m_Gut.png"), plot = p1, width = 8.5, height = 10, units = "cm")


#################### Supplementary Figure 25d  ############################
### spatialAnno with 8markers of Presomitic from embryo1
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
library(Seurat)
### cellassign
load("embryo3_allcelltype_unknown_8markers_cellassign2_8marker_presomitic_8m_Gut.RData")

celltype = fit$cell_type
celltype[celltype=="other"] = "Unknown"
celltype_reorder_rename = sort(unique(celltype))

dat = data.frame(pos[,1], pos[,2], factor(celltype, levels=sort(unique(celltype)), labels=celltype_reorder_rename))
names(dat)= c("imagerow", "imagecol", "cluster")

index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y2 = y[index1]

idx_col = match(sort(unique(celltype)), sort(unique(y2)))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=cluster)) +
  geom_point(size = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 2))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080")[idx_col])
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo3_all_celltypes_unknown_8marker_cellassign_without_legend_8marker_presomitic_8m_Gut.png"), plot = p1, width = 8.5, height = 10, units = "cm")



#################### Supplementary Figure 25d  ############################
### spatialAnno with 8markers of Presomitic from embryo1
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/")
library(Seurat)
### Garnett
load("embryo3_allcelltype_unknown_8markers_scSorter_initial_icmem_lfc=_8marker_presomitic_8m_Gut_adjusted.RData")

library(garnett)
embryo1_cds = readRDS("embryo3_garnett.rds")
layer_Garnett = pData(embryo1_cds)[,5]
celltype = layer_Garnett
celltype_reorder_rename = sort(unique(celltype))
celltype_reorder_rename[celltype_reorder_rename=="unknown"] = "Unknown"
dat = data.frame(pos[,1], pos[,2], factor(celltype, levels=sort(unique(celltype)), labels=celltype_reorder_rename))
names(dat)= c("imagerow", "imagecol", "cluster")

index1 = which(y!="Low quality" &  y!="ExE endoderm" & y!="Blood progenitors")
y2 = y[index1]

idx_col = match(sort(unique(celltype)), sort(unique(y2)))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=3000-imagecol, color=cluster)) +
  geom_point(size = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 2))+ 
  scale_color_manual(values=c(hue_pal()(21),"#808080")[idx_col])
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/embryo2/figure/embryo3_all_celltypes_unknown_8marker_garnett_without_legend_8marker_presomitic_8m_Gut.png"), plot = p1, width = 8.5, height = 10, units = "cm")



