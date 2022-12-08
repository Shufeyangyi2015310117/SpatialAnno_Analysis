# before run this code, get input from Brain12summary.R

# 1-------------------------------------------------------------------------
#rm(list=ls())
setwd("./Real_data_results/dataFiles/DLPFC/Brain12cross/")
library(reshape)
library(ggplot2)
library(ggsci)
library(patchwork)
library(systemfonts)    
library(grid)

theme_Publication <- function(base_size=14, legend.position = "bottom", base_family="Helvetica") {
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
               legend.text = element_text(size=18),
               legend.position = legend.position,
               legend.direction = "horizontal",
               legend.key.size= unit(0.4, "cm"),
               legend.spacing  = unit(0, "cm"),
               legend.title = element_blank(),
               plot.margin=unit(c(10,5,5,0),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))
      
}

ss_4 <- readRDS("brain12_semisupervised.rds")
us_2 <- readRDS("brain12_unsupervised.rds")
method <- c("SpatialAnno", "scSorter", "SCINA", "Garnett","CellAssign")#, "BayesSpace", "DR-SC")
name_ID12 <- as.character(c(151507, 151508, 151509, 151510, 151669, 151670,
                              151671, 151672, 151673, 151674, 151675, 151676))
samp <- rep(1:3, each=4)
# ------------------------------------------------------------------------------------------ 
####################################################################################
####################################################################################
##################### Figure 2a

kp <- data.frame(ID=name_ID12, SpatialAnno=ss_4$sp$kp_sp[, 1], scSorter=ss_4$sc$kp_sc[, 1], 
                SCINA = ss_4$scina$kp_scina[, 1], Garnett = ss_4$garnett$kp_garnett[, 1], CellAssign = ss_4$cellassign$kp_cellassign[, 1],
                #BayesSpace = us_2$bs$kp_bs[, 1], 'DR-SC' = us_2$drsc$kp_drsc[, 1], 
                check.names=FALSE)
df_kp <- melt(data = kp, id.vars = c("ID"), variable_name="Method")
#colnames(df_kp)[3] <- "Kappa"
df_kp$measure <- "Kappa"

mF1 <- data.frame(ID=name_ID12, SpatialAnno=ss_4$sp$f1_sp[, 1], scSorter=ss_4$sc$f1_sc[, 1], 
                SCINA = ss_4$scina$f1_scina[, 1], Garnett = ss_4$garnett$f1_garnett[, 1], CellAssign = ss_4$cellassign$f1_cellassign[, 1],
                #BayesSpace = us_2$bs$f1_bs[, 1], 'DR-SC' = us_2$drsc$f1_drsc[, 1], 
                check.names=FALSE)
df_mF1 <- melt(data = mF1, id.vars = c("ID"), variable_name="Method")
#colnames(df_mF1)[3] <- "mF1"
df_mF1$measure <- "mF1"


ACC <- data.frame(ID=name_ID12, SpatialAnno=ss_4$sp$ac_sp[, 1], scSorter=ss_4$sc$ac_sc[, 1], 
                SCINA = ss_4$scina$ac_scina[, 1], Garnett = ss_4$garnett$ac_garnett[, 1], CellAssign = ss_4$cellassign$ac_cellassign[, 1],
                #BayesSpace = us_2$bs$ac_bs[, 1], 'DR-SC' = us_2$drsc$ac_drsc[, 1], 
                check.names=FALSE)
df_ACC <- melt(data = ACC, id.vars = c("ID"), variable_name="Method")
#colnames(df_ACC)[3] <- "ACC"
df_ACC$measure <- "ACC"

# two measures in main text.
df <- rbind(df_kp, df_mF1)
df$measure <- factor(df$measure, level=c("Kappa", "mF1", "ACC")[-3], order=T)

p <- ggplot(data=df, aes(x = Method, y = value, fill = Method)) +   
 geom_boxplot(aes(fill=Method)) +
 coord_cartesian(ylim=c(0, .8)) + 
facet_grid(~measure) +  scale_fill_simpsons() +
theme_Publication(base_size=24) + guides(fill = guide_legend(ncol = 3, byrow = TRUE))
 

pdf("brain_kp_f1.pdf") #保存为pdf
p
dev.off()


## three measures
df <- rbind(df_kp, df_mF1, df_ACC)
df$measure <- factor(df$measure, level=c("Kappa", "mF1", "ACC"), order=T)
library(rstatix) 
library(ggpubr)

stat.test <- df %>% group_by(measure) %>%
                wilcox_test(value ~ Method, ref.group = "SpatialAnno")  
              
stat.test <- stat.test %>% add_xy_position(x="Method")

p_all <- ggplot(data=df, aes(x = Method, y = value)) +   
 geom_boxplot(aes(fill = Method)) +
 #coord_cartesian(ylim=c(0, 1)) + 
facet_grid(~measure) +  scale_fill_simpsons() +
  stat_pvalue_manual(stat.test) +
theme_Publication(base_size=24) + guides(fill = guide_legend(ncol = 3, byrow = TRUE))
 

pdf("brain_kp_f1_acc.pdf") #保存为pdf
p_all 
dev.off()


# ------------------------------------------------------------------------------------------ 
## correct vs over-specify
##################### Supplementary Figure 14b ##################### 
sample2 <- readRDS("dlpfc_sample2_CorrectSpecify.rds")
sample2over <- readRDS("dlpfc_sample2_OverSpecify.rds")
kp_c <- data.frame(layer=5, SpatialAnno=sample2$kpC[, 1], scSorter=sample2$kpC[, 2], 
                SCINA = sample2$kpC[, 3], Garnett = sample2$kpC[, 5], CellAssign = sample2$kpC[, 4],
                check.names=FALSE)
kp_o <- data.frame(layer=7, SpatialAnno=sample2over$kp[, 1], scSorter=sample2over$kp[, 2], 
                SCINA = sample2over$kp[, 3], Garnett = sample2over$kp[, 5], CellAssign = sample2over$kp[, 4],
                check.names=FALSE)
kp <- rbind(kp_c, kp_o)
df_kp <- melt(data = kp, id.vars = c("layer"), variable_name="Method")
df_kp$measure <- "Kappa"

mF1_c <- data.frame(layer=5, SpatialAnno=sample2$f1C[, 1], scSorter=sample2$f1C[, 2], 
                SCINA = sample2$f1C[, 3], Garnett = sample2$f1C[, 5], CellAssign = sample2$f1C[, 4],
                check.names=FALSE)
mF1_o <- data.frame(layer=7, SpatialAnno=sample2over$f1[, 1], scSorter=sample2over$f1[, 2], 
                SCINA = sample2over$f1[, 3], Garnett = sample2over$f1[, 5], CellAssign = sample2over$f1[, 4],
                check.names=FALSE)
mF1 <- rbind(mF1_c, mF1_o)
df_mF1 <- melt(data = mF1, id.vars = c("layer"), variable_name="Method")
df_mF1$measure <- "mF1"

ACC_c <- data.frame(layer=5, SpatialAnno=sample2$acC[, 1], scSorter=sample2$acC[, 2], 
                SCINA = sample2$acC[, 3], Garnett = sample2$acC[, 5], CellAssign = sample2$acC[, 4],
                check.names=FALSE)
ACC_o <- data.frame(layer=7, SpatialAnno=sample2over$ac[, 1], scSorter=sample2over$ac[, 2], 
                SCINA = sample2over$ac[, 3], Garnett = sample2over$ac[, 5], CellAssign = sample2over$ac[, 4],
                check.names=FALSE)
ACC <- rbind(ACC_c, ACC_o)
df_ACC <- melt(data = ACC, id.vars = c("layer"), variable_name="Method")
df_ACC$measure <- "ACC"


df <- rbind(df_kp, df_mF1, df_ACC)
df$measure <- factor(df$measure, level=c("Kappa", "mF1", "ACC"), order=T)
stat.test <- df %>% group_by(measure, layer) %>%
                wilcox_test(value ~ Method, ref.group = "SpatialAnno")  
              
stat.test <- stat.test %>% add_xy_position(x="Method")


p <- ggplot(data=df, aes(x = Method, y = value)) +   
 geom_boxplot(aes(fill=Method)) +
 coord_cartesian(ylim=c(0, 0.8)) + 
facet_grid(measure~layer) +  scale_fill_simpsons() +
 #stat_pvalue_manual(stat.test) +
theme_Publication(base_size=28) + guides(fill = guide_legend(ncol = 3, byrow = TRUE))
 
pdf("brainSample2_kp_f1_acc.pdf") #保存为pdf
p
dev.off()
# ------------------------------------------------------------------------------------------ 
################################################Supplementary Figure 14a####################
# markers numbers: top 5 vs 10 vs 15 vs
# remove (r, c) = (1,1), (2,2), (3,3), (4,4), (9,5), (10,6),(11,7), (12,8)
# (c-1) * 12 + r
# idx_rm <- c(1, 14, 27, 40, 57, 70, 83, 96)
kp_5 <- data.frame(top = 5, SpatialAnno=c(ss_4$sp$kp_sp[, 1:8]), scSorter=c(ss_4$sc$kp_sc[, 1:8]), 
                SCINA = c(ss_4$scina$kp_scina[, 1:8]), Garnett = c(ss_4$garnett$kp_garnett[, 1:8]), CellAssign = c(ss_4$cellassign$kp_cellassign[, 1:8]),
                #BayesSpace = c(us_2$bs$kp_bs[, 1:8]), 'DR-SC' = c(us_2$drsc$kp_drsc[, 1:8]), 
                check.names=FALSE)
kp_10 <- data.frame(top = 10, SpatialAnno=c(ss_4$sp$kp_sp[, 9:16]), scSorter=c(ss_4$sc$kp_sc[, 9:16]), 
                SCINA = c(ss_4$scina$kp_scina[, 9:16]), Garnett = c(ss_4$garnett$kp_garnett[, 9:16]), CellAssign = c(ss_4$cellassign$kp_cellassign[, 9:16]),
                #BayesSpace = c(us_2$bs$kp_bs[, 9:16]), 'DR-SC' = c(us_2$drsc$kp_drsc[, 9:16]), 
                check.names=FALSE)
kp_15 <- data.frame(top = 15, SpatialAnno=c(ss_4$sp$kp_sp[, 17:24]), scSorter=c(ss_4$sc$kp_sc[, 17:24]), 
                SCINA = c(ss_4$scina$kp_scina[, 17:24]), Garnett = c(ss_4$garnett$kp_garnett[, 17:24]), CellAssign = c(ss_4$cellassign$kp_cellassign[, 17:24]),
                #BayesSpace = c(us_2$bs$kp_bs[, 17:24]), 'DR-SC' = c(us_2$drsc$kp_drsc[, 17:24]), 
                check.names=FALSE)
kp <- rbind(kp_5, kp_10, kp_15)
df_kp <- melt(data = kp, id.vars = c("top"), variable_name="Method")
df_kp$measure <- "Kappa"

mF1_5 <- data.frame(top = 5, SpatialAnno=c(ss_4$sp$f1_sp[, 1:8]), scSorter=c(ss_4$sc$f1_sc[, 1:8]), 
                SCINA = c(ss_4$scina$f1_scina[, 1:8]), Garnett = c(ss_4$garnett$f1_garnett[, 1:8]), CellAssign = c(ss_4$cellassign$f1_cellassign[, 1:8]),
                #BayesSpace = c(us_2$bs$f1_bs[, 1:8]), 'DR-SC' = c(us_2$drsc$f1_drsc[, 1:8]), 
                check.names=FALSE)
mF1_10 <- data.frame(top = 10, SpatialAnno=c(ss_4$sp$f1_sp[, 9:16]), scSorter=c(ss_4$sc$f1_sc[, 9:16]), 
                SCINA = c(ss_4$scina$f1_scina[, 9:16]), Garnett = c(ss_4$garnett$f1_garnett[, 9:16]), CellAssign = c(ss_4$cellassign$f1_cellassign[, 9:16]),
                #BayesSpace = c(us_2$bs$f1_bs[, 9:16]), 'DR-SC' = c(us_2$drsc$f1_drsc[, 9:16]), 
                check.names=FALSE)
mF1_15 <- data.frame(top = 15, SpatialAnno=c(ss_4$sp$f1_sp[, 17:24]), scSorter=c(ss_4$sc$f1_sc[, 17:24]), 
                SCINA = c(ss_4$scina$f1_scina[, 17:24]), Garnett = c(ss_4$garnett$f1_garnett[, 17:24]), CellAssign = c(ss_4$cellassign$f1_cellassign[, 17:24]),
                #BayesSpace = c(us_2$bs$f1_bs[, 17:24]), 'DR-SC' = c(us_2$drsc$f1_drsc[, 17:24]), 
                check.names=FALSE)
mF1 <- rbind(mF1_5, mF1_10, mF1_15)
df_mF1 <- melt(data = mF1, id.vars = c("top"), variable_name="Method")
df_mF1$measure <- "mF1"


ACC_5 <- data.frame(top = 5, SpatialAnno=c(ss_4$sp$ac_sp[, 1:8]), scSorter=c(ss_4$sc$ac_sc[, 1:8]), 
                SCINA = c(ss_4$scina$ac_scina[, 1:8]), Garnett = c(ss_4$garnett$ac_garnett[, 1:8]), CellAssign = c(ss_4$cellassign$ac_cellassign[, 1:8]),
                #BayesSpace = c(us_2$bs$ac_bs[, 1:8]), 'DR-SC' = c(us_2$drsc$ac_drsc[, 1:8]), 
                check.names=FALSE)
ACC_10 <- data.frame(top = 10, SpatialAnno=c(ss_4$sp$ac_sp[, 9:16]), scSorter=c(ss_4$sc$ac_sc[, 9:16]), 
                SCINA = c(ss_4$scina$ac_scina[, 9:16]), Garnett = c(ss_4$garnett$ac_garnett[, 9:16]), CellAssign = c(ss_4$cellassign$ac_cellassign[, 9:16]),
                #BayesSpace = c(us_2$bs$ac_bs[, 9:16]), 'DR-SC' = c(us_2$drsc$ac_drsc[, 9:16]), 
                check.names=FALSE)
ACC_15 <- data.frame(top = 15, SpatialAnno=c(ss_4$sp$ac_sp[, 17:24]), scSorter=c(ss_4$sc$ac_sc[, 17:24]), 
                SCINA = c(ss_4$scina$ac_scina[, 17:24]), Garnett = c(ss_4$garnett$ac_garnett[, 17:24]), CellAssign = c(ss_4$cellassign$ac_cellassign[, 17:24]),
                #BayesSpace = c(us_2$bs$ac_bs[, 17:24]), 'DR-SC' = c(us_2$drsc$ac_drsc[, 17:24]), 
                check.names=FALSE)
ACC <- rbind(ACC_5, ACC_10, ACC_15)
df_ACC <- melt(data = ACC, id.vars = c("top"), variable_name="Method")
df_ACC$measure <- "ACC"


df <- rbind(df_kp, df_mF1, df_ACC)
df$measure <- factor(df$measure, level=c("Kappa", "mF1", "ACC"), order=T)
stat.test <- df %>% group_by(measure, top) %>%
                wilcox_test(value ~ Method, ref.group = "SpatialAnno")  
stat.test <- stat.test %>% add_xy_position(x="Method")
stat.test$p.scient <- format(stat.test$p, scientific = TRUE)
    


p <- ggplot(data=df, aes(x = Method, y = value)) +   
 geom_boxplot(aes(fill=Method)) +
 coord_cartesian(ylim=c(0, 1.14)) + 
facet_grid(measure~top) +  scale_fill_simpsons() +  
 stat_pvalue_manual(stat.test, label="p.scient", size = 2.8) +
theme_Publication(base_size=16) #+ guides(fill = guide_legend(ncol = 3, byrow = TRUE))
 
pdf("brainAcross_kp_f1_acc.pdf") #保存为pdf
p
dev.off()
 

########################################################################################################################
# 2. plot  Figure 2b, Supplementary Figure 2b-13b
library(BayesSpace)
library(ggplot2)
library(dplyr)
library(ggsci)
library(patchwork)
K_id <- c(rep(7, 4), rep(5, 4), rep(7, 4))


iter <- 12# 151510 sample
 
  name_ID12 <- as.character(c(151507, 151508, 151509, 151510, 
                            151669, 151670, 151671, 151672, 
                            151673, 151674, 151675, 151676))
  ID <- name_ID12[iter] 
  color <- pal_d3("category10")(8)
  names(color) <- c("Layer1","Layer2","Layer3","Layer4","Layer5","Layer6","WM","Unknown")
  #if(iter %in% 5:8) color <- pal_d3("category10")(8)[3:8]
  dir.exp <- "./Real_data_results/dataFiles/DLPFC/Brain12cross/"
  marker_dir <- "./Real_data_analysis/DLPFC/Datasets/markerfiles/"
  id <- "151507"
  top <- 5
  K <- K_id[iter]

  dlpfcSCINA <- readRDS(paste0(dir.exp, ID, "SCINA_Garnett.rds") )
  
   
  y <- as.character(dlpfcSCINA$layer_guess_reordered)
  y[is.na(y)] <- 'Unknown'
  y_vec <- factor(y, level=c("Layer1","Layer2","Layer3","Layer4","Layer5","Layer6","WM","Unknown"), order=T)

  p0 <- clusterPlot(dlpfcSCINA, label=y_vec, palette=NULL, size=0.05) +
    labs(title=paste0(ID)) + scale_fill_manual(values=color[levels(y_vec)[table(y_vec) != 0]]) + 
    theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), legend.position="none") 

 # plot legend
    loc1=dlpfcSCINA$row; loc2=dlpfcSCINA$col
    datt = data.frame(cluster=y_vec, loc1=dlpfcSCINA$row, loc2=dlpfcSCINA$col)
    p0 = ggplot(datt, aes(x = loc1, y = loc2, color = cluster)) +  geom_point( alpha =1, size=0.45) +
          #scale_color_manual(values = color, drop=FALSE)
           scale_color_manual(values=color[levels(y_vec)[table(y_vec) != 0]]) 
pdf(file = "dlpfc_ground_truth_legend.pdf",width = 5, height = 4)
p <- p0 + theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position="bottom",
        legend.title = element_blank(),
        axis.ticks = element_blank()) + 
          guides(colour = guide_legend(override.aes = list(size = 5))) 
legend <- cowplot::get_legend(p)  
grid.newpage()
grid.draw(legend)
dev.off()



  layer_sp <- colData(dlpfcSCINA)[paste0("icmemlfc", id, "top", top)][[1]]
  sp_vec <- factor(layer_sp, level=c("Layer1","Layer2","Layer3","Layer4","Layer5","Layer6","WM","Unknown"), order=T)
  p1 <- clusterPlot(dlpfcSCINA, label=sp_vec, palette=NULL, size=0.05) +
   labs(title="SpatialAnno") + scale_fill_manual(values=color[levels(sp_vec)[table(sp_vec) != 0]]) + 
    theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), legend.position="none")


  layer_sc <- colData(dlpfcSCINA)[paste0("layer_scSorter", id,"top", top)][[1]]
  sc_vec <- factor(layer_sc, level=c("Layer1","Layer2","Layer3","Layer4","Layer5","Layer6","WM","Unknown"), order=T)
  p2 <- clusterPlot(dlpfcSCINA, label= sc_vec, palette=NULL, size=0.05) +
    labs(title=paste0("scSorter")) + scale_fill_manual(values=color[levels(sc_vec)[table(sc_vec) != 0]]) + 
     theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), legend.position="none")
   


  layer_scina <- colData(dlpfcSCINA)[paste0("layer_scina", id,"top", top)][[1]] 
  layer_scina[layer_scina == "unknown"] <- 'Unknown'
  scina_vec <- factor(layer_scina, level=c("Layer1","Layer2","Layer3","Layer4","Layer5","Layer6","WM","Unknown"), order=T)
  p3 <- clusterPlot(dlpfcSCINA, label= scina_vec, palette=NULL, size=0.05) +
    labs(title= "SCINA") + scale_fill_manual(values=color[levels(scina_vec)[table(scina_vec) != 0]]) + 
     theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), legend.position="none")
 
  layer_garnett <- colData(dlpfcSCINA)[paste0("garnett", id, "top", top)][[1]]
  if (ID %in% name_ID12[5:8]) {dlpfcC <- readRDS(paste0(dir.exp, ID, "CorrectSpecify_Garnett.rds")) 
        layer_garnett <- colData(dlpfcC)[paste0("garnett", id, "top", top)][[1]]}
  layer_garnett[layer_garnett == "unknown"] <- 'Unknown'
  garnett_vec <- factor(layer_garnett, level=c("Layer1","Layer2","Layer3","Layer4","Layer5","Layer6","WM","Unknown"), order=T)
  p7 <- clusterPlot(dlpfcSCINA, label= garnett_vec, palette=NULL, size=0.05) +
    labs(title= "Garnett") + scale_fill_manual(values=color[levels(garnett_vec)[table(garnett_vec) != 0]]) + 
     theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), legend.position="none")
 

  layer_cellassign <- colData(dlpfcSCINA)[paste0("layer_cellassign", id,"top", top)][[1]] 
  cellassign_vec <- factor(layer_cellassign, level=c("Layer1","Layer2","Layer3","Layer4","Layer5","Layer6","WM","Unknown"), order=T)
  p4 <- clusterPlot(dlpfcSCINA, label= cellassign_vec, palette=NULL, size=0.05) +
     labs(title= "CellAssign") +  scale_fill_manual(values=color[levels(cellassign_vec)[table(cellassign_vec) != 0]]) + 
      theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), legend.position="none")
   
 # ----------------  
  # prepared for PAGA
  dlpfc <- readRDS(paste0("./DLPFC/unsupervised/", ID, ".rds"))   
  # bayesspace  
  ## low dimensional representation
  pc <- reducedDims(dlpfcSCINA)$PCA  
  write.table(pc, paste0("./DLPFC/unsupervised/", ID, "_byesSpace_Ez.txt"), row.names =FALSE, col.names =FALSE)   

  # dr-sc
  ## low dimensional representation
  output <- readRDS(paste0("./DLPFC/unsupervised/", ID, "_DR_SC.rds"))
  hZ <- output$hZ
  write.table(hZ, paste0("./DLPFC/unsupervised/", ID, "_drsc_Ez.txt"), row.names =FALSE, col.names =FALSE)   
 
# ----------------

  p <- p0 + p1 + p2 + p3 + p7 + p4 + plot_layout(byrow = T, nrow=2, ncol=3)

  #p <- p0 + p1 + p2 + p3 + p7 + p4 + p5 + p6 + plot_layout(byrow = T, nrow=2, ncol=4)
  pdf(file = paste0("Brain_heatmap", ID, ".pdf"))  
    p
  dev.off()

# plot each cluster seperately
 

 
########################################################################################################################
### 3.  
library(ggplot2) 
library(Seurat)
library(ggsci)
library(patchwork)
#ari <- matrix(0, ncol = 3, nrow =12)


iter <- 10# 151510 sample
color <- pal_d3("category10")(8)
names(color) <- c("Layer1","Layer2","Layer3","Layer4","Layer5","Layer6","WM","Unknown")

name_ID12 <- as.character(c(151507, 151508, 151509, 151510, 
                            151669, 151670, 151671, 151672, 
                            151673, 151674, 151675, 151676))


ID <- name_ID12[iter] 
dir.exp <- "./Real_data_results/dataFiles/DLPFC/Brain12cross/"
dir.spark <- "./Real_data_analysis/DLPFC/Datasets/brain12_spark/"
marker_dir <- "./Real_data_analysis/DLPFC/Datasets/markerfiles/"
marker_file <- paste0("151507", "_DEGmarkerTop_", 5, ".rds")
markers <- readRDS(paste0(marker_dir, marker_file))
if (ID %in% name_ID12[5:8]) {
    markers$Layer1 = NULL
    markers$Layer2 = NULL
}  
fit_s <- readRDS(paste0(dir.exp, ID, "_icmem.rds"))
dlpfc <- readRDS(paste0(dir.exp, ID, "Manual.rds"))

library(Seurat)
seu <- as.Seurat(dlpfc)
# -----------------------------------
#pca + BayesSPace
layer_bs <- read.table(paste0("./DLPFC/unsupervised/", "dlpfc_", ID, "_byesSpace_celltype.txt"))
bs_vec <- factor(layer_bs$V1, 
                       level=c("Layer1","Layer2","Layer3","Layer4","Layer5","Layer6","WM","Unknown"), order=T)
Idents(seu) <-  bs_vec
pc <- reducedDims(dlpfc)$PCA
seu@reductions$"pca" <- CreateDimReducObject(embeddings = pc, key='PC_', assay=DefaultAssay(seu))
seu <- RunTSNE(seu, reduction="pca", dims=1:15, verbose = F, check_duplicates = FALSE)
ncol(seu)
length(unique(colnames(seu)))
p3 <- TSNEPlot(seu)
tSNE_pca <- seu@reductions$"tsne"
seu@reductions$"tSNE_pca" <- tSNE_pca
p3 <- DimPlot(seu, cells=which(layer_bs != "Unknown"), reduction = 'tSNE_pca', pt.size = 1, cols= color[levels(bs_vec)[table(bs_vec) != 0]]) + 
        ggtitle("PCA") + theme(plot.title = element_text(face = "bold", size = rel(2), hjust = 0.5),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.title = element_text(size=25),
                axis.line = element_line(colour="black"),
                legend.text = element_text(size=25),
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                legend.title = element_blank())

# -----------------------------------
# spatialAno
sp_vec <- factor(fit_s$anno, 
                       level=c("Layer1","Layer2","Layer3","Layer4","Layer5","Layer6","WM","Unknown"), order=T)
Idents(seu) <- sp_vec
#-

#-
hZ <- fit_s$Ez_u
row.names(hZ) <- colnames(seu)
colnames(hZ) <- paste0('spanno_', 1:ncol(hZ))

seu@reductions$"spanno" <- CreateDimReducObject(embeddings = hZ, key='spanno_', assay=DefaultAssay(seu))
seu <- RunTSNE(seu, reduction="spanno", dims=1:15, verbose = F, check_duplicates = FALSE)
ncol(seu)
length(unique(colnames(seu)))
p1 <- TSNEPlot(seu)
tSNE_spanno <- seu@reductions$"tsne"
seu@reductions$"tSNE_spanno" <- tSNE_spanno
p1 <- DimPlot(seu, cells=which(fit_s$anno != "Unknown"), reduction = 'tSNE_spanno', pt.size = 1, cols= color[levels(sp_vec)[table(sp_vec) != 0]]) + 
        ggtitle("SpatialAnno") + theme(plot.title = element_text(face = "bold", size = rel(2), hjust = 0.5),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.title = element_text(size=25),
                axis.line = element_line(colour="black"),
                legend.text = element_text(size=25),
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                legend.title = element_blank())
# -----------------------------------
# drsc
layer_drsc <- read.table(paste0("./Real_data_results/dataFiles/DLPFC/unsupervised/", "dlpfc_", ID, "_drsc_celltype.txt"))
Ez <-  data.matrix(read.table(paste0("./Real_data_results/dataFiles/DLPFC/unsupervised/dlpfc_", ID, "_drsc_Ez.txt")))
drsc_vec <- factor(layer_drsc$V1, 
                       level=c("Layer1","Layer2","Layer3","Layer4","Layer5","Layer6","WM","Unknown"), order=T)
Idents(seu) <- drsc_vec 

row.names(Ez) <- colnames(seu)
colnames(Ez) <- paste0('drsc_', 1:ncol(Ez))

seu@reductions$"drsc" <- CreateDimReducObject(embeddings = Ez, key='drsc_', assay=DefaultAssay(seu))
seu <- RunTSNE(seu, reduction="drsc", dims=1:15, verbose = F, check_duplicates = FALSE)
ncol(seu)
length(unique(colnames(seu)))

p2 <- TSNEPlot(seu)
tSNE_drsc <- seu@reductions$"tsne"
seu@reductions$"tSNE_drsc" <- tSNE_drsc
p2 <- DimPlot(seu, cells=which(layer_drsc != "Unknown"), reduction = 'tSNE_drsc', pt.size = 1, cols= color[levels(drsc_vec)[table(drsc_vec) != 0]])+ 
        ggtitle("DR-SC") + theme(plot.title = element_text(face = "bold", size = rel(2), hjust = 0.5),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.title = element_text(size=25),
                axis.line = element_line(colour="black"),
                legend.text = element_text(size=25),
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                legend.title = element_blank())
#p2 
#dev.off()
# -----------------------------------
p <- p1 + p3 + p2 + plot_layout(byrow = T, nrow=1, ncol=3) 
 
ggsave(file=paste0('Brain_', ID, '_TSNE.pdf'), plot = p,
       width = 24, height =8, units = "in", dpi = 1000)
# -----------------------------------
# -----------------------------------
########################### ARI in Figure 2d
y <- as.character(dlpfc$layer_guess_reordered)
y[is.na(y)] <- 'Unknown'

seu <- FindNeighbors(seu, reduction = "spanno", dims = 1:15)
seu <- FindClusters(seu)  
ari[iter, 1] <- mclust::adjustedRandIndex(seu$seurat_clusters, y)

seu <- FindNeighbors(seu, reduction = "pca", dims = 1:15)
seu <- FindClusters(seu)  
ari[iter, 2] <- mclust::adjustedRandIndex(seu$seurat_clusters, y)

seu <- FindNeighbors(seu, reduction = "drsc", dims = 1:15)
seu <- FindClusters(seu)  
ari[iter, 3] <- mclust::adjustedRandIndex(seu$seurat_clusters, y)


# RGB plot in Figure 2d, Supplementary Figure 2c-13c
library(SpatialPCA)
set.seed(1234)
xy_coords <- data.frame(
        x_coord = seu@meta.data[, c("imagecol")], 
        y_coord = -seu@meta.data[, c("imagerow")]
    )
p1 = plot_RGB_tSNE(xy_coords, t(hZ), pointsize=3)[[2]] + theme(plot.title = element_blank())
p2 = plot_RGB_tSNE(xy_coords, t(pc), pointsize=3)[[2]] + theme(plot.title = element_blank())
p3 = plot_RGB_tSNE(xy_coords, t(Ez), pointsize=3)[[2]] + theme(plot.title = element_blank())
#pdf("LIBD_RGB_tSNE_1x.pdf",width=12, height=5)
pp <- p1  + p3  + p2  + plot_layout(byrow = T, nrow=1, ncol=3)  
ggsave(file=paste0('Brain_', ID, 'RGB_TSNE.pdf'), plot = pp, title="",
       width = 24, height =8)
if(FALSE) { 
  p1 = plot_RGB_UMAP(xy_coords, t(hZ), pointsize=3)[[2]] + theme(plot.title = element_blank())
  p2 = plot_RGB_UMAP(xy_coords, t(pc), pointsize=3)[[2]] + theme(plot.title = element_blank())
  p3 = plot_RGB_UMAP(xy_coords, t(Ez), pointsize=3)[[2]] + theme(plot.title = element_blank())
  pp <- p1  + p3  + p2  + plot_layout(byrow = T, nrow=1, ncol=3)  
  ggsave(file=paste0('Brain_', ID, 'RGB_UMAP.pdf'), plot = pp, title="",
       width = 24, height =8)
}

#
saveRDS(ari, "Brain12embedingARI.rds")
########################################################################################################################
# PAGA graph in Figure 2e, Supplementary Figure 2c-13c
### 4. run paga in python (dlpfc.py, dlpfc-drsc.py, dlpfc-bayesSpace.py) 



########################################################################################################################
### 5. # RGB plot in Figure 2c, Supplementary Figure 2a-13a
library(BayesSpace)
library(ggplot2)
library(dplyr)
library(ggsci)
library("org.Hs.eg.db")
library(patchwork)
K_id <- c(rep(7, 4), rep(5, 4), rep(7, 4))


iter <- 9# 151510 sample
 
  name_ID12 <- as.character(c(151507, 151508, 151509, 151510, 
                            151669, 151670, 151671, 151672, 
                            151673, 151674, 151675, 151676))
  ID <- name_ID12[iter] 
  color <- pal_d3("category10")(8)
  names(color) <- c("Layer1","Layer2","Layer3","Layer4","Layer5","Layer6","WM","Unknown")
  dir.exp <- "./Real_data_results/dataFiles/DLPFC/Brain12cross/"
  marker_dir <- "./Real_data_analysis/DLPFC/Datasets/markerfiles/"
  id <- "151507"
  top <- 5
  K <- K_id[iter] 
  dlpfcSCINA <- readRDS(paste0(dir.exp, ID, "SCINA_Garnett.rds") )
   
  y <- as.character(dlpfcSCINA$layer_guess_reordered)
  y[is.na(y)] <- 'Unknown'
  y_vec <- factor(y, level=c("Layer1","Layer2","Layer3","Layer4","Layer5","Layer6","WM","Unknown"), order=T)


  layer_sp <- colData(dlpfcSCINA)[paste0("icmemlfc", id, "top", top)][[1]]
  if(K == 7) Layers <- c("Layer1","Layer2","Layer3","Layer4","Layer5","Layer6","WM")
  if(K == 5) Layers <- c("Layer3","Layer4","Layer5","Layer6","WM")
pl <- list()
for(k in 1:K) {
  layer_k <- Layers[k]
  Lay <- layer_sp
  Lay[layer_sp !=  layer_k] <- "Unknown"
  pl[[k]] <-  clusterPlot(dlpfcSCINA, label=Lay, palette=NULL, size=0.05) +
    labs(title=layer_k) + scale_fill_manual(values=color[levels(y_vec)[table(y_vec) != 0]]) + 
     theme(plot.title = element_text(size = rel(0.8), hjust = 0.5), legend.position="none")
}

# feature plot
marker_file <- paste0(id, "_DEGmarkerTop_", top, ".rds")
markers <- readRDS(paste0(marker_dir, marker_file))
if (ID %in% name_ID12[5:8]) {
    markers$Layer1 = NULL
    markers$Layer2 = NULL
}  

#mapIds(org.Hs.eg.db, keys = markers[[7]], keytype = "ENSEMBL", column="SYMBOL")
firstup <- function(X) {
  x <- tolower(X)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
#L1 4 : AQP4
#L2 1 : HPCAL1 
#L3 1 : CARTPT
#L4 1 : NEFH
#L5 1 : PCP4
#L6 1 : KRT17
#WM 3 : MOBP
if(K==5) f_idx <- c(1,1,1,1,3) else f_idx <- c(4, 1, 1, 1, 1, 1, 3)
fl <- list()
for(k in 1:K) {
  gene <- markers[[Layers[k]]][f_idx[k]]
  gname <- firstup(mapIds(org.Hs.eg.db, keys = gene, keytype = "ENSEMBL", column="SYMBOL"))
fl[[k]] <-  featurePlot(dlpfcSCINA, feature=gene, color=NA) +  theme(plot.title = element_text(size = rel(0.8), hjust = 0.5, face = 'italic'), legend.position="none") +
   labs(title=gname)  
}
if(iter == 9){ # plot legend
  k <- 2
  gene <- markers[[Layers[k]]][f_idx[k]]
  fl_2 <-  featurePlot(dlpfcSCINA, feature=gene, color=NA) +  theme(plot.title = element_text(size = rel(0.8), hjust = 0.5, face = 'italic'), legend.position="top") +
  pdf(file = paste0("Brain_GE_lengend", ID, ".pdf"))  
  legend <- cowplot::get_legend(fl_2)  
  grid.newpage()
  grid.draw(legend)
  dev.off()
}

if(K==7) pf <-  fl[[1]] + fl[[2]] +  fl[[3]] + fl[[4]] + fl[[5]] + fl[[6]] + fl[[7]] + pl[[1]] + pl[[2]] +  pl[[3]] + pl[[4]] + pl[[5]] + pl[[6]] + pl[[7]]+ plot_layout(byrow = T, nrow=2)
  if(K==5) pf <- fl[[1]] + fl[[2]] +  fl[[3]] + fl[[4]] + fl[[5]] + pl[[1]] + pl[[2]] +  pl[[3]] + pl[[4]] + pl[[5]] + plot_layout(byrow = T, nrow=2)

pdf(file = paste0("Brain_Layer_marker", ID, ".pdf"))  
    pf
dev.off()

