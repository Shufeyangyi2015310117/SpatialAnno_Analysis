#######################  theme  ############################
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





#######################  Figure 3a  ############################ 
## ground truth
## rm(list = ls())
library(Seurat)
library(scales)
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/")
load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc=0.7172.RData")
dim(pos)
fit_s2_type = colnames(rho)[fit_s2$type]
fit_s2_type[fit_s2_type=="GC"] = "GCL"
fit_s2_type[fit_s2_type=="M/TC"] = "MCL"
fit_s2_type[fit_s2_type=="OSNs"] = "ONL"
fit_s2_type[fit_s2_type=="PGC"] = "GL"
fit_s2_type[fit_s2_type=="EPL-IN"] = "GL"
fit_s2_type[fit_s2_type=="Unknown"] = "MCL"

## manual annotation
## mapping the original coordinate to the integer grid
location = round(pos, 0)

for (k in c(22,24,25)){
  fit_s2_type[location[,1] == 10 & location[,2] == k] = "GL"
}
for (k in c(21,27)){
  fit_s2_type[location[,1] == 11 & location[,2] == k] = "GL"
}
for (k in c(17)){
  fit_s2_type[location[,1] == 11 & location[,2] == k] = "Unknown"
}
for (k in c(19,28)){
  fit_s2_type[location[,1] == 12 & location[,2] == k] = "GL"
}
for (k in c(22,26)){
  fit_s2_type[location[,1] == 12 & location[,2] == k] = "MCL"
}
for (k in c(15,17,18)){
  fit_s2_type[location[,1] == 13 & location[,2] == k] = "GL"
}
for (k in c(27)){
  fit_s2_type[location[,1] == 13 & location[,2] == k] = "MCL"
}
for (k in c(20)){
  fit_s2_type[location[,1] == 14 & location[,2] == k] = "MCL"
}
for (k in c(15)){
  fit_s2_type[location[,1] == 14 & location[,2] == k] = "GL"
}
for (k in c(29)){
  fit_s2_type[location[,1] == 15 & location[,2] == k] = "GL"
}
for (k in c(28,16)){
  fit_s2_type[location[,1] == 15 & location[,2] == k] = "MCL"
}
for (k in c(29)){
  fit_s2_type[location[,1] == 16 & location[,2] == k] = "GL"
}
for (k in c(18)){
  fit_s2_type[location[,1] == 17 & location[,2] == k] = "MCL"
}
for (k in c(14)){
  fit_s2_type[location[,1] == 17 & location[,2] == k] = "Unknown"
}
for (k in c(25,18)){
  fit_s2_type[location[,1] == 18 & location[,2] == k] = "GL"
}
for (k in c(20)){
  fit_s2_type[location[,1] == 18 & location[,2] == k] = "Unknown"
}
for (k in c(24,21,20)){
  fit_s2_type[location[,1] == 19 & location[,2] == k] = "ONL"
}
for (k in c(24)){
  fit_s2_type[location[,1] == 20 & location[,2] == k] = "GL"
}
for (k in c(18,20:25)){
  fit_s2_type[location[,1] == 21 & location[,2] == k] = "GL"
}
for (k in c(16)){
  fit_s2_type[location[,1] == 22 & location[,2] == k] = "GL"
}
for (k in c(27)){
  fit_s2_type[location[,1] == 23 & location[,2] == k] = "Unknown"
}
for (k in c(17)){
  fit_s2_type[location[,1] == 26 & location[,2] == k] = "GL"
}
for (k in c(28)){
  fit_s2_type[location[,1] == 27 & location[,2] == k] = "MCL"
}
for (k in c(27,21,20)){
  fit_s2_type[location[,1] == 27 & location[,2] == k] = "GL"
}
for (k in c(25:27,23)){
  fit_s2_type[location[,1] == 28 & location[,2] == k] = "GL"
}
for (k in c(24:26)){
  fit_s2_type[location[,1] == 29 & location[,2] == k] = "GL"
}
for (k in c(19)){
  fit_s2_type[location[,1] == 30 & location[,2] == k] = "ONL"
}
Rep12_MOB_manual_annotation = fit_s2_type
save(Rep12_MOB_manual_annotation, file = "Rep12_MOB_manual_annotation.RData")
dat = data.frame(pos[,1], pos[,2], factor(fit_s2_type, levels = c("GCL", "MCL", "ONL", "GL", "Unknown")))

colnames(dat)= c("imagerow", "imagecol", "Layer")

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=imagecol, color=Layer)) +
  geom_point(size = 3, alpha = 0.5) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) + 
  coord_cartesian(xlim = c(7, 30), ylim = c(13, 30)) + 
  scale_x_continuous(breaks = seq(5, 31, by = 1)) +
  scale_y_continuous(breaks = seq(12, 31, by = 1)) +
  guides(colour = guide_legend(override.aes = list(size = 5)))+ 
  scale_color_manual(values=c(hue_pal()(4),"#808080"))
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_ground_truth.pdf"), plot = p1, width = 15 , height = 10, units = "cm")


library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=imagecol, color=Layer)) +
  geom_point(size = 3, alpha = 0.5) +
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
  scale_color_manual(values=c(c(colfunc(100)[(7)*5], brewer.pal(9, "Greens")[5],hue_pal()(4)[4], brewer.pal(9, "Blues")[7]), "#808080"))
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_ground_truth.pdf"), plot = p1, width = 8, height = 6, units = "cm")








#######################  Figure 3a  ############################ 
## code on mac
## icmem
## rm(list = ls())
library(Seurat)
library(scales)
library(RColorBrewer)
colfunc <- colorRampPalette(c("red", "white"))
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/")
load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc_replace_with_Penk_Th.RData")
dim(pos)
fit_s2_type = colnames(rho)[fit_s2$type]
dat = data.frame(pos[,1], pos[,2], factor(fit_s2_type))
names(dat)= c("imagerow", "imagecol", "Cell type")

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=imagecol, color=`Cell type`)) +
  geom_point(size = 3, alpha = 0.5) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5)))+ 
  scale_color_manual(values=c(c("#FFD700", colfunc(100)[(7)*5], brewer.pal(9, "Greens")[5],hue_pal()(4)[4], brewer.pal(9, "Blues")[7]), "#808080"))
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_4markers_icmem_lfc_replace_with_Th_Penk.pdf"), plot = p1, width = 11, height = 6, units = "cm")
library(mclust)
mclust::adjustedRandIndex(fit_s2_type, Rep12_MOB_manual_annotation)


# Using the cowplot package
legend <- cowplot::get_legend(p1)
pdf(file = "/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_celltype_legend_replace_with_Th_Penk.pdf",width = 1, height = 2)
grid.newpage()
grid.draw(legend)
dev.off()


#######################  Figure 3a  ############################ 
## code on mac
## icmem
## rm(list = ls())
library(Seurat)
library(scales)
library(RColorBrewer)
colfunc <- colorRampPalette(c("red", "white"))
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/")
load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc_replace_with_Penk_Th.RData")
dim(pos)
fit_s2_type = colnames(rho)[fit_s2$type]
dat = data.frame(pos[,1], pos[,2], factor(fit_s2_type))
names(dat)= c("imagerow", "imagecol", "Cell type")

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=imagecol, color=`Cell type`)) +
  geom_point(size = 3, alpha = 0.5) +
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
  scale_color_manual(values=c(c("#FFD700", colfunc(100)[(7)*5], brewer.pal(9, "Greens")[5],hue_pal()(4)[4], brewer.pal(9, "Blues")[7]), "#808080"))
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_4markers_icmem_lfc_replace_with_Th_Penk.pdf"), plot = p1, width = 8, height = 6, units = "cm")
library(mclust)
mclust::adjustedRandIndex(fit_s2_type, Rep12_MOB_manual_annotation)



#######################  Figure 3a  ############################ 
## code on mac
## SCINA
## rm(list = ls())
library(Seurat)
library(scales)
library(RColorBrewer)
colfunc <- colorRampPalette(c("red", "white"))
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/")
load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc_replace_with_Penk_Th.RData")
dim(pos)
fit_s2_type = results$cell_labels
fit_s2_type[fit_s2_type=="unknown"] = "Unknonw"
dat = data.frame(pos[,1], pos[,2], factor(fit_s2_type))
names(dat)= c("imagerow", "imagecol", "cluster")

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=imagecol, color=cluster)) +
  geom_point(size = 3, alpha = 0.5) +
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
  scale_color_manual(values=c(c("#FFD700", colfunc(100)[(7)*5], brewer.pal(9, "Greens")[5],hue_pal()(4)[4], brewer.pal(9, "Blues")[7]), "#808080"))
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_4markers_SCINA_replace_with_Th_Penk.pdf"), plot = p1, width = 8, height = 6, units = "cm")
mclust::adjustedRandIndex(fit_s2_type, Rep12_MOB_manual_annotation)





#######################  Figure 3a  ############################ 
## code on mac
## scSorter
## rm(list = ls())
library(Seurat)
library(scales)
library(RColorBrewer)
colfunc <- colorRampPalette(c("red", "white"))
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/")
load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc_replace_with_Penk_Th.RData")
dim(pos)
fit_s2_type = rts$Pred_Type
dat = data.frame(pos[,1], pos[,2], factor(fit_s2_type))
names(dat)= c("imagerow", "imagecol", "cluster")

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=imagecol, color=cluster)) +
  geom_point(size = 3, alpha = 0.5) +
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
  scale_color_manual(values=c(c("#FFD700", colfunc(100)[(7)*5], brewer.pal(9, "Greens")[5],hue_pal()(4)[4], brewer.pal(9, "Blues")[7]), "#808080"))
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_4markers_scSorter_Th_Penk.pdf"), plot = p1, width = 8, height = 6, units = "cm")




#######################  Figure 3a  ############################ 
## code on mac
## cellassign
## rm(list = ls())
library(Seurat)
library(scales)
library(RColorBrewer)
colfunc <- colorRampPalette(c("red", "white"))
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/")
load("Rep12_MOB_4markers_cellassign2_replace_with_Th_Penk.RData")
dim(pos)
fit_s2_type = fit$cell_type
fit_s2_type[fit_s2_type == "other"] = "Unknown"
dat = data.frame(pos[,1], pos[,2], factor(fit_s2_type))
names(dat)= c("imagerow", "imagecol", "cluster")

idx_col = match(sort(unique(fit_s2_type)), c("EPL-IN", "GC", "M/TC", "OSNs", "PGC", "Unknown"))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=imagecol, color=cluster)) +
  geom_point(size = 3, alpha = 0.5) +
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
  scale_color_manual(values=c(c(brewer.pal(3, "Set3")[2], colfunc(100)[(7)*5], brewer.pal(9, "Greens")[5],hue_pal()(4)[4], brewer.pal(9, "Blues")[7]), "#808080")[idx_col])
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_4markers_cellassign_replace_with_Th_Penk.pdf"), plot = p1, width = 8, height = 6, units = "cm")
mclust::adjustedRandIndex(fit_s2_type, Rep12_MOB_manual_annotation)



#######################  Figure 3a  ############################ 
## code on mac
## Garnett
## rm(list = ls())
library(Seurat)
library(scales)
library(RColorBrewer)
colfunc <- colorRampPalette(c("red", "white"))
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/")
load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc=0.7172.RData")
MOB_cds = readRDS("MOB_garnett_5unknown_Th_Penk.rds")
library(garnett)
layer_Garnett = pData(MOB_cds)[,5]
dim(pos)
fit_s2_type = layer_Garnett
dat = data.frame(pos[,1], pos[,2], factor(fit_s2_type))
names(dat)= c("imagerow", "imagecol", "cluster")

idx_col = match(sort(unique(fit_s2_type)), c("EPL-IN", "GC", "M/TC", "OSNs", "PGC", "Unknown"))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=imagecol, color=cluster)) +
  geom_point(size = 3, alpha = 0.5) +
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
  scale_color_manual(values=c(c("#FFD700", colfunc(100)[(7)*5], brewer.pal(9, "Greens")[5],hue_pal()(4)[4], brewer.pal(9, "Blues")[7]), "#808080")[idx_col])
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_4markers_Garnett_replace_with_Th_Penk.pdf"), plot = p1, width = 8, height = 6, units = "cm")





#######################  Figure 3b  ############################ 
## calculate Kappa f1s acc 
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2")
load("Rep12_MOB_manual_annotation.RData")
load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc_replace_with_Penk_Th.RData")
load("Rep12_MOB_4markers_cellassign2_replace_with_Th_Penk.RData")
load("Rep12_MOB_4markers_adj_DRSC2.RData")
load("Rep12_MOB_4markers_Seurat_BayesSpace2.RData")
MOB_cds = readRDS("MOB_garnett_5unknown_Th_Penk.rds")
library(garnett)
layer_Garnett = pData(MOB_cds)[,5]

Pre_cell_type = matrix(0, length(Rep12_MOB_manual_annotation), 7)
Pre_cell_type[,1] = colnames(rho)[fit_s2$type]
Pre_cell_type[,2] = results$cell_labels
Pre_cell_type[,3] = rts$Pred_Type
Pre_cell_type[,4] = fit$cell_type
Pre_cell_type[,5] = layer_drsc
Pre_cell_type[,6] = layer_BayesSpace
Pre_cell_type[,7] = layer_Garnett
Rep12_MOB_manual_annotation_matched = Rep12_MOB_manual_annotation
Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="GCL"] = "GC"
Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="MCL"] = "M/TC"
Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="ONL"] = "OSNs"
Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="GL"] = "PGC"

source("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/function_f1score.R")
library(psych)
value = matrix(0, 7, 4)
for (i in 1:7){
  idx1 = which(Rep12_MOB_manual_annotation_matched!="Unknown")
  idx2 = which(Pre_cell_type[,i]!="Unknown" & Pre_cell_type[,i]!="unknown" & Pre_cell_type[,i]!="other")
  value[i,1] = cohen.kappa(x=cbind(Rep12_MOB_manual_annotation_matched[intersect(idx1,idx2)], Pre_cell_type[intersect(idx1,idx2),i]))$kappa
  value[i,2] = mean(evaluate(Rep12_MOB_manual_annotation_matched[intersect(idx1,idx2)], Pre_cell_type[intersect(idx1,idx2),i])$F1)
  value[i,3] = mean(Rep12_MOB_manual_annotation_matched[intersect(idx1,idx2)] == Pre_cell_type[intersect(idx1,idx2),i])
  value[i,4] = mclust::adjustedRandIndex(Rep12_MOB_manual_annotation_matched, Pre_cell_type[,i])
}

colnames(value) = c("Kappa", "mF1", "Acc", "ARI")
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
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_4markers_measure3_replace_with_Th_Penk.pdf"), plot = p, width = 25, height = 20, units = "cm")




#################### Supplementary Figure 15c  ############################ 
## calculate Kappa f1s acc 
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2")
load("Rep12_MOB_manual_annotation.RData")
load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc_replace_with_Penk_Th.RData")
load("Rep12_MOB_4markers_cellassign2_replace_with_Th_Penk.RData")
load("Rep12_MOB_4markers_adj_DRSC2.RData")
load("Rep12_MOB_4markers_Seurat_BayesSpace2.RData")
MOB_cds = readRDS("MOB_garnett_5unknown_Th_Penk.rds")
library(garnett)
layer_Garnett = pData(MOB_cds)[,5]

Pre_cell_type = matrix(0, length(Rep12_MOB_manual_annotation), 7)
Pre_cell_type[,1] = colnames(rho)[fit_s2$type]
Pre_cell_type[,2] = results$cell_labels
Pre_cell_type[,3] = rts$Pred_Type
Pre_cell_type[,4] = fit$cell_type
Pre_cell_type[,5] = layer_drsc
Pre_cell_type[,6] = layer_BayesSpace
Pre_cell_type[,7] = layer_Garnett
Pre_cell_type[Pre_cell_type == "EPL-IN"] = "PGC"


Rep12_MOB_manual_annotation_matched = Rep12_MOB_manual_annotation
Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="GCL"] = "GC"
Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="MCL"] = "M/TC"
Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="ONL"] = "OSNs"
Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="GL"] = "PGC"

source("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/function_f1score.R")
library(psych)
value = matrix(0, 7, 4)
for (i in 1:7){
  idx1 = which(Rep12_MOB_manual_annotation_matched!="Unknown")
  idx2 = which(Pre_cell_type[,i]!="Unknown" & Pre_cell_type[,i]!="unknown" & Pre_cell_type[,i]!="other")
  value[i,1] = cohen.kappa(x=cbind(Rep12_MOB_manual_annotation_matched[intersect(idx1,idx2)], Pre_cell_type[intersect(idx1,idx2),i]))$kappa
  value[i,2] = mean(evaluate(Rep12_MOB_manual_annotation_matched[intersect(idx1,idx2)], Pre_cell_type[intersect(idx1,idx2),i])$F1)
  value[i,3] = mean(Rep12_MOB_manual_annotation_matched[intersect(idx1,idx2)] == Pre_cell_type[intersect(idx1,idx2),i])
  value[i,4] = mclust::adjustedRandIndex(Rep12_MOB_manual_annotation_matched, Pre_cell_type[,i])
}

colnames(value) = c("Kappa", "mF1", "acc", "ARI")
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
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_4markers_measure4_replace_with_Th_Penk.pdf"), plot = p, width = 30, height = 20, units = "cm")






#################### Supplementary Figure 15d  ############################ 
## calculate Kappa f1s acc 
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2")
load("Rep12_MOB_manual_annotation.RData")
load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc_replace_with_Penk_Th.RData")
load("Rep12_MOB_4markers_cellassign2_replace_with_Th_Penk.RData")
load("Rep12_MOB_4markers_adj_DRSC2.RData")
load("Rep12_MOB_4markers_Seurat_BayesSpace2.RData")
MOB_cds = readRDS("MOB_garnett_5unknown_Th_Penk.rds")
library(garnett)
layer_Garnett = pData(MOB_cds)[,5]

Pre_cell_type = matrix(0, length(Rep12_MOB_manual_annotation), 7)
Pre_cell_type[,1] = colnames(rho)[fit_s2$type]
Pre_cell_type[,2] = results$cell_labels
Pre_cell_type[,3] = rts$Pred_Type
Pre_cell_type[,4] = fit$cell_type
Pre_cell_type[,5] = layer_drsc
Pre_cell_type[,6] = layer_BayesSpace
Pre_cell_type[,7] = layer_Garnett
Pre_cell_type[Pre_cell_type == "EPL-IN"] = "PGC"
Pre_cell_type[Pre_cell_type == "PGC"] = "M/TC"


Rep12_MOB_manual_annotation_matched = Rep12_MOB_manual_annotation
Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="GCL"] = "GC"
Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="MCL"] = "M/TC"
Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="ONL"] = "OSNs"
Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="GL"] = "PGC"
Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="PGC"] = "M/TC"

source("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/function_f1score.R")
library(psych)
value = matrix(0, 7, 4)
for (i in 1:7){
  idx1 = which(Rep12_MOB_manual_annotation_matched!="Unknown")
  idx2 = which(Pre_cell_type[,i]!="Unknown" & Pre_cell_type[,i]!="unknown" & Pre_cell_type[,i]!="other")
  value[i,1] = cohen.kappa(x=cbind(Rep12_MOB_manual_annotation_matched[intersect(idx1,idx2)], Pre_cell_type[intersect(idx1,idx2),i]))$kappa
  value[i,2] = mean(evaluate(Rep12_MOB_manual_annotation_matched[intersect(idx1,idx2)], Pre_cell_type[intersect(idx1,idx2),i])$F1)
  value[i,3] = mean(Rep12_MOB_manual_annotation_matched[intersect(idx1,idx2)] == Pre_cell_type[intersect(idx1,idx2),i])
  value[i,4] = mclust::adjustedRandIndex(Rep12_MOB_manual_annotation_matched, Pre_cell_type[,i])
}

colnames(value) = c("Kappa", "mF1", "acc", "ARI")
rownames(value) = c("SpatialAnno", "SCINA", "scSorter", "CellAssign", "DR-SC", "BayesSpace", "Garnett")
value2 = value[c(1:4,7),1:3]

library(reshape2) 
df = melt(value2)
colnames(df) = c("method", "measurement", "score")
df$method = factor(df$method, levels = c("SpatialAnno", "scSorter", "SCINA","Garnett", "CellAssign"))

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
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_4markers_measure5_replace_with_Th_Penk.pdf"), plot = p, width = 30, height = 20, units = "cm")




#################### Supplementary Figure 15a  ############################
getdata = function(){
  ################################# KAPPA and M1F #############################3
  ## calculate Kappa f1s acc 
  setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2")
  load("Rep12_MOB_manual_annotation.RData")
  load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc_replace_with_Penk_Th_Mural1_EC1.RData")
  load("Rep12_MOB_4markers_cellassign_replace_with_Th_Penk_Mural1_EC1.RData")
  MOB_cds = readRDS("MOB_garnett_Th_Penk_Mural1_EC1.rds")
  library(garnett)
  layer_Garnett = pData(MOB_cds)[,5]
  
  Pre_cell_type = matrix(0, length(Rep12_MOB_manual_annotation), 5)
  Pre_cell_type[,1] = colnames(rho)[fit_s2$type]
  Pre_cell_type[,2] = results$cell_labels
  Pre_cell_type[,3] = rts$Pred_Type
  Pre_cell_type[,4] = fit$cell_type
  Pre_cell_type[,5] = layer_Garnett
  Rep12_MOB_manual_annotation_matched = Rep12_MOB_manual_annotation
  Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="GCL"] = "GC"
  Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="MCL"] = "M/TC"
  Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="ONL"] = "OSNs"
  Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="GL"] = "PGC"
  
  source("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/function_f1score.R")
  library(psych)
  value = matrix(0, 5, 4)
  for (i in 1:5){
    idx1 = which(Rep12_MOB_manual_annotation_matched!="Unknown")
    idx2 = which(Pre_cell_type[,i]!="Unknown" & Pre_cell_type[,i]!="unknown" & Pre_cell_type[,i]!="other")
    value[i,1] = cohen.kappa(x=cbind(Rep12_MOB_manual_annotation_matched[intersect(idx1,idx2)], Pre_cell_type[intersect(idx1,idx2),i]))$kappa
    value[i,2] = mean(evaluate(Rep12_MOB_manual_annotation_matched[intersect(idx1,idx2)], Pre_cell_type[intersect(idx1,idx2),i])$F1)
    value[i,3] = mean(Rep12_MOB_manual_annotation_matched[intersect(idx1,idx2)] == Pre_cell_type[intersect(idx1,idx2),i])
    value[i,4] = mclust::adjustedRandIndex(Rep12_MOB_manual_annotation_matched, Pre_cell_type[,i])
  }
  
  colnames(value) = c("Kappa", "mF1", "acc", "ARI")
  rownames(value) = c("SpatialAnno", "SCINA", "scSorter", "CellAssign", "Garnett")
  value2 = value[c(1:5),1:3]
  
  library(reshape2) 
  df = melt(value2)
  colnames(df) = c("method", "measurement", "score")
  df$method = factor(df$method, levels = c("SpatialAnno", "scSorter", "SCINA", "Garnett", "CellAssign"))
  return(df)
}

getdata2 = function(){
  ################################# KAPPA and M1F #############################3
  ## calculate Kappa f1s acc 
  setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2")
  load("Rep12_MOB_manual_annotation.RData")
  load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc_replace_with_Penk_Th.RData")
  load("Rep12_MOB_4markers_cellassign2_replace_with_Th_Penk.RData")
  load("Rep12_MOB_4markers_adj_DRSC2.RData")
  load("Rep12_MOB_4markers_Seurat_BayesSpace2.RData")
  MOB_cds = readRDS("MOB_garnett_5unknown_Th_Penk.rds")
  library(garnett)
  layer_Garnett = pData(MOB_cds)[,5]
  
  Pre_cell_type = matrix(0, length(Rep12_MOB_manual_annotation), 7)
  Pre_cell_type[,1] = colnames(rho)[fit_s2$type]
  Pre_cell_type[,2] = results$cell_labels
  Pre_cell_type[,3] = rts$Pred_Type
  Pre_cell_type[,4] = fit$cell_type
  Pre_cell_type[,5] = layer_drsc
  Pre_cell_type[,6] = layer_BayesSpace
  Pre_cell_type[,7] = layer_Garnett
  Rep12_MOB_manual_annotation_matched = Rep12_MOB_manual_annotation
  Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="GCL"] = "GC"
  Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="MCL"] = "M/TC"
  Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="ONL"] = "OSNs"
  Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="GL"] = "PGC"
  
  source("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/function_f1score.R")
  library(psych)
  value = matrix(0, 7, 4)
  for (i in 1:7){
    idx1 = which(Rep12_MOB_manual_annotation_matched!="Unknown")
    idx2 = which(Pre_cell_type[,i]!="Unknown" & Pre_cell_type[,i]!="unknown" & Pre_cell_type[,i]!="other")
    value[i,1] = cohen.kappa(x=cbind(Rep12_MOB_manual_annotation_matched[intersect(idx1,idx2)], Pre_cell_type[intersect(idx1,idx2),i]))$kappa
    value[i,2] = mean(evaluate(Rep12_MOB_manual_annotation_matched[intersect(idx1,idx2)], Pre_cell_type[intersect(idx1,idx2),i])$F1)
    value[i,3] = mean(Rep12_MOB_manual_annotation_matched[intersect(idx1,idx2)] == Pre_cell_type[intersect(idx1,idx2),i])
    value[i,4] = mclust::adjustedRandIndex(Rep12_MOB_manual_annotation_matched, Pre_cell_type[,i])
  }
  
  colnames(value) = c("Kappa", "mF1", "acc", "ARI")
  rownames(value) = c("SpatialAnno", "SCINA", "scSorter", "CellAssign", "DR-SC", "BayesSpace", "Garnett")
  value2 = value[c(1:4,7),1:3]
  
  library(reshape2) 
  df = melt(value2)
  colnames(df) = c("method", "measurement", "score")
  df$method = factor(df$method, levels = c("SpatialAnno", "scSorter", "SCINA", "Garnett", "CellAssign"))
  return(df)
}

df1 = getdata()
df1$num_celltype = "7"
df2 = getdata2()
df2$num_celltype = "5"
df = rbind(df1, df2)
# Basic barplot
library(ggplot2)
library(scales)
library(RColorBrewer)
library(psych)
library(ggsci)
library(systemfonts)
library(ggplot2)
p<-ggplot(data=df, aes(x= method, y=score, fill = method)) +
  geom_bar(stat="identity")+
  facet_grid(num_celltype~measurement) + 
  scale_fill_simpsons() + 
  theme_Publication(base_size=30) + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_4markers_measure12_celltype5-7_replace_with_Th_Penk.pdf"), plot = p, width = 25, height = 20, units = "cm")




#######################  Figure 3c  ############################

setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2")
load("Rep12_MOB_count_matrix-1.RData")

MOB <- CreateSeuratObject(counts = MOB_raw, project = "MOB")
MOB = SCTransform(MOB)
X = MOB@assays$SCT@scale.data
load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc_replace_with_Th_Penk.RData")
dim(pos)
#reducedim(MOB) = pos
celltype = sort(unique(colnames(rho)[fit_s2$type]))
genes = c("Kit", "Penk","Cdhr1", "S100a5", "Th")
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
med <- apply(MOB[["SCT"]]@scale.data[c("Kit", "Penk","Cdhr1", "S100a5", "Th"),], 1, quantile, quant)
dat$gene = factor(dat$gene, levels = genes)
dat$marker = "marker"
p1 <- ggplot(dat, aes(x=imagerow, y=imagecol, color=cluster))+
  geom_point(size = 3) + 
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
        axis.ticks.y=element_blank())

ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_4markers_celltype_specific_marker_with_legend_replace_with_Penk_Th.pdf"), plot = p1, width = 40, height = 7, units = "cm")


# Using the cowplot package
legend <- cowplot::get_legend(p1)
pdf(file = "/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_4markers_celltype_specific_marker_legend_replace_with_Penk_Th.pdf",width = 2, height = 1)
grid.newpage()
grid.draw(legend)
dev.off()



library(Seurat)
library(ggplot2)
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2")
load("Rep12_MOB_count_matrix-1.RData")

MOB <- CreateSeuratObject(counts = MOB_raw, project = "MOB")
MOB = SCTransform(MOB)
X = MOB@assays$SCT@scale.data
load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc_replace_with_Penk_Th.RData")
dim(pos)
#reducedim(MOB) = pos
celltype = sort(unique(colnames(rho)[fit_s2$type]))
genes = c("Kit", "Penk","Cdhr1", "S100a5", "Th")
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
med <- apply(MOB[["SCT"]]@scale.data[c("Kit", "Penk","Cdhr1", "S100a5", "Th"),], 1, quantile, quant)
dat$gene = factor(dat$gene, levels = genes)
dat$marker = "marker"
p1 <- ggplot(dat, aes(x=imagerow, y=imagecol, color=cluster))+
  geom_point(size = 3) + 
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
        strip.text.x = element_text(face = "italic"))

ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_4markers_celltype_specific_marker_replace_with_Penk_Th.pdf"), plot = p1, width = 40, height = 7, units = "cm")



#################### Supplementary Figure 17  ############################

setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2")
load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc_replace_with_Penk_Th.RData")
load("Rep12_MOB_4markers_cellassign2_replace_with_Th_Penk.RData")
MOB_cds = readRDS("MOB_garnett_5unknown_Th_Penk.rds")
method = c("SpatialAnno", "scSorter", "CellAssign", "SCINA", "Garnett")
dim(pos)
library(scales)
library(RColorBrewer)
colfunc <- colorRampPalette(c("red", "white"))
library(ggplot2)
#celltype = sort(unique(colnames(rho)[fit_s2$type]))
data = list()
for (k in 1:5){
  if (k == 1){
    fit_s2_type = colnames(rho)[fit_s2$type]
    unique_celltype_sort = sort(unique(fit_s2_type))
  }
  if (k == 2){
    fit_s2_type = rts$Pred_Type
    unique_celltype_sort = sort(unique(fit_s2_type))
  }
  if (k == 3){
    fit_s2_type = fit$cell_type
    fit_s2_type[fit_s2_type=="other"] = "Unknown"
    unique_celltype_sort = sort(unique(fit_s2_type))
    unique_celltype_sort = c("EPL-IN", "GC", unique_celltype_sort)
  }
  if (k == 4){
    fit_s2_type = results$cell_labels
    unique_celltype_sort = sort(unique(fit_s2_type))
  }
  if (k == 5){
    fit_s2_type = pData(MOB_cds)[,5]
    unique_celltype_sort = sort(unique(fit_s2_type))
  }
  
  plist =list()
  for (i in 1:5){
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
df3 = data[[3]]
df3$method = method[3] 
df4 = data[[4]]
df4$method = method[4] 
df5 = data[[5]]
df5$method = method[5] 

df = rbind(df1, df2, df3, df4, df5)
df$method = factor(df$method, levels = c("SpatialAnno", "scSorter", "CellAssign", "SCINA", "Garnett"))

p1 <- ggplot(df, aes(x=imagerow, y=imagecol, color=cluster))+
  geom_point(size = 3, alpha = 1) +
  facet_grid(method~celltype)+
  theme_Publication()+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+ 
  scale_color_manual(values=c("#FFD700", "#808080", colfunc(100)[(7)*5], brewer.pal(9, "Greens")[5],hue_pal()(4)[4], brewer.pal(9, "Blues")[7]))
#  p1 = cowplot::plot_grid(plotlist = plist, nrow = 1)
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_4markers_celltype_separately2_all_replace_with_Penk_Th.pdf"), plot = p1, width = 40, height = 30, units = "cm")




#######################  Figure 3c  ############################

setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2")
load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc_replace_with_Penk_Th.RData")
MOB_cds = readRDS("MOB_garnett_5unknown_Th_Penk.rds")
method = c("spatialAnno", "scSorter", "Garnett")
dim(pos)
library(scales)
library(RColorBrewer)
colfunc <- colorRampPalette(c("red", "white"))
library(ggplot2)
#celltype = sort(unique(colnames(rho)[fit_s2$type]))
data = list()
for (k in 1:3){
  if (k == 1){
    fit_s2_type = colnames(rho)[fit_s2$type]
    unique_celltype_sort = sort(unique(fit_s2_type))
  }
  if (k == 2){
    fit_s2_type = rts$Pred_Type
    unique_celltype_sort = sort(unique(fit_s2_type))
  }
  if (k == 3){
    fit_s2_type = pData(MOB_cds)[,5]
    unique_celltype_sort = sort(unique(fit_s2_type))
  }
  plist =list()
  for (i in 1:5){
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
df3 = data[[3]]
df3$method = method[3] 
df = df1
df$method = factor(df$method, levels = "spatialAnno")

p1 <- ggplot(df, aes(x=imagerow, y=imagecol, color=cluster))+
  geom_point(size = 3, alpha = 1) +
  facet_grid(method~celltype)+
  theme_Publication()+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+ 
  scale_color_manual(values=c("#FFD700", "#808080", colfunc(100)[(7)*5], brewer.pal(9, "Greens")[5],hue_pal()(4)[4], brewer.pal(9, "Blues")[7]))
#  p1 = cowplot::plot_grid(plotlist = plist, nrow = 1)
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_4markers_celltype_separately_garnett_replace_with_Penk_Th.pdf"), plot = p1, width = 40, height = 7, units = "cm")



#################### Supplementary Figure 15b  ############################
## The effect of the number of non-markers on the performance of SpatialAnno 
rm(list = ls())
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2")
load("Rep12_MOB_manual_annotation.RData")
Pre_cell_type = matrix(0, 282, 9)
load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc_replace_with_Penk_Th.RData")
Pre_cell_type[,1] = colnames(rho)[fit_s2$type]
Pre_cell_type[,2] = rts$Pred_Type
load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc_top300_replace_with_Th_Penk.RData")
Pre_cell_type[,3] = colnames(rho)[fit_s2$type]
load("Rep12_MOB_4markers_scSorter_top300_replace_with_Th_Penk.RData")
Pre_cell_type[,4] = rts$Pred_Type
load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc_top30_replace_with_Th_Penk.RData")
Pre_cell_type[,5] = colnames(rho)[fit_s2$type]
load("Rep12_MOB_4markers_scSorter_top30_replace_with_Th_Penk.RData")
Pre_cell_type[,6] = rts$Pred_Type
MOB_cds_top3000 = readRDS("MOB_garnett_top3000_with_Penk_Th.rds")
Pre_cell_type[,7] = pData(MOB_cds_top3000)[,5]
MOB_cds_top300 = readRDS("MOB_garnett_top300_with_Penk_Th.rds")
Pre_cell_type[,8] = pData(MOB_cds_top300)[,5]
MOB_cds_top30 = readRDS("MOB_garnett_top30_with_Penk_Th.rds")
load("MOB_garnett_top30.RData")
Pre_cell_type[idx,9] = pData(MOB_cds_top30)[,5]
Pre_cell_type[-idx,9] = "Unknown"

Rep12_MOB_manual_annotation_matched = Rep12_MOB_manual_annotation
Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="GCL"] = "GC"
Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="MCL"] = "M/TC"
Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="ONL"] = "OSNs"
Rep12_MOB_manual_annotation_matched[Rep12_MOB_manual_annotation_matched=="GL"] = "PGC"
source("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/function_f1score.R")
library(psych)
value = matrix(0, 9, 4)
for (i in 1:9){
  idx1 = which(Rep12_MOB_manual_annotation_matched!="Unknown")
  idx2 = which(Pre_cell_type[,i]!="Unknown")
  value[i,1] = cohen.kappa(x=cbind(Rep12_MOB_manual_annotation_matched[intersect(idx1,idx2)], Pre_cell_type[intersect(idx1,idx2),i]))$kappa
  value[i,2] = mean(evaluate(Rep12_MOB_manual_annotation_matched[intersect(idx1,idx2)], Pre_cell_type[intersect(idx1,idx2),i])$F1)
  value[i,3] = mean(Rep12_MOB_manual_annotation_matched[intersect(idx1,idx2)] == Pre_cell_type[intersect(idx1,idx2),i])
  value[i,4] = mclust::adjustedRandIndex(Rep12_MOB_manual_annotation_matched, Pre_cell_type[,i])
}

colnames(value) = c("Kappa", "mF1", "acc", "ARI")
rownames(value) = c("SpatialAnno", "scSorter", "SpatialAnno", "scSorter", "SpatialAnno", "scSorter", "Garnett", "Garnett", "Garnett")
value2 = value[,1:3]

library(reshape2) 
df = melt(value2)
colnames(df) = c("method", "measurement", "score")
df$method = factor(df$method, levels = c("SpatialAnno", "scSorter", "Garnett"))
df$nonmarker = rep(c(rep(c("3000", "300", "30"), each = 2),"3000","300","30"), time = 3)
library(ggsci)
library(systemfonts)
library(ggplot2)

# Basic barplot
p<-ggplot(data=df, aes(x= method, y=score, fill = nonmarker)) +
  geom_bar(stat="identity", position=position_dodge())+
  facet_grid(.~measurement, switch = "y")+
  scale_fill_simpsons() + 
  theme_Publication(base_size=30) +
  theme(axis.text.x = element_text(size=16),
        legend.title = element_text(size=28)) + 
  guides(fill=guide_legend(title="# of non-marker"))
ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_4markers_measure8_top30-3000_replace_with_Th_Penk.pdf"), plot = p, width = 30, height = 20, units = "cm")





#######################  Figure 3d  ######################################
## RGB plot
## remove unknowns followed by tSNE
library(scater)
library(DR.SC)
library(iDR.SC)
library(Seurat)
colfunc <- colorRampPalette(c("red", "white"))
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/")
load("Rep12_MOB_manual_annotation.RData")
load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc_replace_with_Penk_Th.RData")
load("Rep12_MOB_4markers_adj_DRSC.RData")
load("Rep12_MOB_4markers_adj_DRSC2.RData")
load("Rep12_MOB_4markers_Seurat_BayesSpace2.RData")
load("Rep12_MOB_4markers_Seurat_PCA.RData")
method = c("ICMEM", "DRSC", "BayesSpace")

ARI = matrix(0,3,1)

for (k in 1){
  if (k == 1){
    celltype = colnames(rho)[fit_s2$type]
    Rep12_MOB_Ez = fit_s2$Ez_u
  }
  if (k == 2){
    celltype = layer_drsc 
    library(DR.SC)
    drsc_out <- DR.SC::selectModel(srsc, criteria = 'MBIC', pen.const=1)
    Rep12_MOB_Ez = drsc_out$hZ
  }
  if (k == 3){
    celltype = layer_BayesSpace
    Rep12_MOB_Ez = MOB_PCA[,1:15]
  }
  
  
  library(mclust)
  
  if (is.null(rownames(Rep12_MOB_Ez)) == TRUE){
    rownames(Rep12_MOB_Ez) = paste0("cell",1:dim(Rep12_MOB_Ez)[1])
  }
  
  if (is.null(colnames(Rep12_MOB_Ez)) == TRUE){
    colnames(Rep12_MOB_Ez) = paste0("PC",1:dim(Rep12_MOB_Ez)[2])
  }
  
  seu = CreateSeuratObject(
    counts = t(Rep12_MOB_Ez),
    project = "CreateSeuratObject",
    assay = "RNA",
    names.field = 1,
    names.delim = "_",
    meta.data = NULL,
  )
  seu[["pca"]] <- CreateDimReducObject(embeddings = Rep12_MOB_Ez, key = "PC_", assay = DefaultAssay(seu))
  
  seu <- FindNeighbors(seu,reduction = "pca", dims = 1:15)
  seu = FindClusters(seu)
  
  idx = which(Rep12_MOB_manual_annotation!="Unknown" & celltype!="Unknown")
  ARI[k] = mclust::adjustedRandIndex(seu$seurat_clusters[idx], Rep12_MOB_manual_annotation[idx])
  
  tsne3dim = scater::calculateTSNE(t(Rep12_MOB_Ez), ncomponents = 3)


  pList <- iDR.SC::plot_RGB(pos, tsne3dim, pointsize = 2) +
    mytheme_graybox()+ xlim(c(7,30)) + ylim(c(12,30)) + geom_point(alpha=0.5)

  ggsave(file=paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_RBG_", method[k], "_replace_with_Th_Penk.pdf"), plot = pList,
         width = 3, height = 2, units = "in", dpi = 500)
}





#################### Supplementary Figure 15e  ############################
## Rep12 MOB
## remove unknowns followed by tSNE
library(scater)
colfunc <- colorRampPalette(c("red", "white"))
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/")
load("Rep12_MOB_4markers_SCINA_initial_icmem_lfc_replace_with_Penk_Th.RData")
load("Rep12_MOB_4markers_adj_DRSC.RData")
load("Rep12_MOB_4markers_adj_DRSC2.RData")
load("Rep12_MOB_4markers_Seurat_BayesSpace2.RData")
load("Rep12_MOB_4markers_Seurat_PCA.RData")
method = c("ICMEM", "DRSC", "BayesSpace")

for (k in 1){
  if (k == 1){
    celltype = colnames(rho)[fit_s2$type]
    Rep12_MOB_Ez = fit_s2$Ez_u
  }
  if (k == 2){
    celltype = layer_drsc 
    library(DR.SC)
    drsc_out <- selectModel(srsc, criteria = 'MBIC', pen.const=1)
    Rep12_MOB_Ez = drsc_out$hZ
  }
  if (k == 3){
    celltype = layer_BayesSpace
    Rep12_MOB_Ez = MOB_PCA[,1:15]
  }
  
  fit = scater::calculateTSNE(t(as.matrix(Rep12_MOB_Ez)))
  
  idx_col = match(sort(unique(celltype)), c("EPL-IN", "GC", "M/TC", "OSNs", "PGC", "Unknown"))
  
  ### tsne plot
  dat = as.data.frame(fit)
  colnames(dat) = c("X", "Y")
  dat$cluster = as.matrix(celltype)
  library(ggplot2)
  p1 <- ggplot(dat, aes(x=X, y=Y, color=cluster)) +
    geom_point(size = 1, alpha = 1) +
    theme_Publication()+
    theme(legend.position="right",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank())+
    guides(colour = guide_legend(override.aes = list(size = 5), ncol=1)) +
    scale_color_manual(values=c(c("#FFD700", colfunc(100)[(7)*5], brewer.pal(9, "Greens")[5],hue_pal()(4)[4], brewer.pal(9, "Blues")[7], "#808080")[idx_col]))
  ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/figure/Rep12_MOB_", method[k], "_followed_by_tsne_replace_with_Th_Penk.pdf"), plot = p1, width = 16, height = 14, units = "cm")
  
}



#################### Supplementary Figure 16  ############################

getdat <- function(filename, method = "SpatialAnno"){
  for (k in 1:12){
    if (k == 1){
      load(paste0("./output/Rep", k, substr(filename, 5, 100)))
      dim(pos)
      if (method == "SpatialAnno"){
        fit_s2_type = colnames(rho)[fit_s2$type]
      }else if (method == "SCINA"){
        fit_s2_type = results$cell_labels
        fit_s2_type[fit_s2_type=="unknown"] = "Unknown"
      }else if(method == "scSorter"){
        fit_s2_type = rts$Pred_Type 
      }else if(method == "cellassign"){
        fit_s2_type = fit$cell_type
        fit_s2_type[fit_s2_type=="other"] = "Unknown"
      }else if(method == "Garnett"){
        fit_s2_type = pData(MOB_cds)[,5]
      }
      dat = data.frame(pos[,1], pos[,2], fit_s2_type)
      names(dat)= c("imagerow", "imagecol", "Cell type")
      dat$Rep = paste0("Rep",k)
    }else{
      load(paste0("./output/Rep", k, substr(filename, 5, 100)))
      dim(pos)
      if (method == "SpatialAnno"){
        fit_s2_type = colnames(rho)[fit_s2$type]
      }else if (method == "SCINA"){
        fit_s2_type = results$cell_labels
        fit_s2_type[fit_s2_type=="unknown"] = "Unknown"
      }else if(method == "scSorter"){
        fit_s2_type = rts$Pred_Type 
      }else if(method == "cellassign"){
        fit_s2_type = fit$cell_type
        fit_s2_type[fit_s2_type=="other"] = "Unknown"
      }else if(method == "Garnett"){
        fit_s2_type = pData(MOB_cds)[,5]
      }
      dat2 = data.frame(pos[,1], pos[,2], fit_s2_type)
      names(dat2)= c("imagerow", "imagecol", "Cell type")
      dat2$Rep = paste0("Rep",k)
      dat = rbind(dat, dat2)
    }
  }
  dat$`Cell type` = factor(dat$`Cell type`)
  dat$Rep = factor(dat$Rep, levels = paste0("Rep",1:12))
  return(dat)
}


################### All Rep MOB ###########
## All Rep MOB
## SpatialAnno scSorter SCINA
###########################################
library(Seurat)
library(scales)
library(RColorBrewer)
colfunc <- colorRampPalette(c("red", "white"))
setwd("/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB3/")
figure_path = "~/Desktop/2020-9-19/2022/Xingjie/MOB3/figure/"

method = c("SpatialAnno", "SCINA", "scSorter", "cellassign", "Garnett")
k=1
filename = "Rep1_MOB_4markers_SCINA_initial_icmem_lfc.RData"
dat = getdat(filename, method = method[k])
dat$Method = "SpatialAnno"
k=2
dat2 = getdat(filename, method = method[k])
dat2$Method = "SCINA"
dat = rbind(dat, dat2)
k=3
dat2 = getdat(filename, method = method[k])
dat2$Method = "scSorter"
dat = rbind(dat, dat2)
k=4
filename = "Rep1_MOB_4markers_cellassign.RData"
dat2 = getdat(filename, method = method[k])
dat2$Method = "CellAssign"
dat = rbind(dat, dat2)
k=5
filename = "Rep1_MOB_4markers_Garnett_lfc.RData"
dat2 = getdat(filename, method = method[k])
dat2$Method = "Garnett"
dat = rbind(dat, dat2)

dat$Method = factor(dat$Method, levels = c("SpatialAnno", "SCINA", "scSorter", "CellAssign", "Garnett"))

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagerow, y=imagecol, color=`Cell type`)) +
  geom_point(size = 2, alpha = 0.5) +
  ggh4x::facet_grid2(Rep~Method , scales = "free", independent = "all", switch = "y")+
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 5)))+ 
  scale_color_manual(values=c(c("#FFD700", colfunc(100)[(7)*5], brewer.pal(9, "Greens")[5],hue_pal()(4)[4], brewer.pal(9, "Blues")[7]), "#808080"))
ggsave(paste0(figure_path, "All_Rep_MOB_4markers.pdf"), plot = p1, width = 25, height = 50, units = "cm")










