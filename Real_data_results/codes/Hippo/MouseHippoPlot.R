setwd("./Real_data_results/dataFiles/Hippo/")
library(ggplot2)
library(patchwork)
library(dplyr)
library(grid)
library(Seurat)
embro <- 1 # 26177cells:36h, 43994 cells:48h 
iter <- 1
dataSource <- 1
mk_n <-  c(5, 10)
 

slide.seq <- readRDS(paste0("hippoAll_anno_", "satija", mk_n[iter], "_", embro, ".rds"))
types <- sort(unique(slide.seq@meta.data[, "icmemlfc"]))

corlors = c("#8B2029","#08F5EE","#FFFC3D","#7E6CF4","#51AF50",
            "#FC6143","#FC4326","#3F967E","#8A4D6B","#C1E1FC",
            "#03FF7B","#9A4CA7","#CC6454","#FF7900", "#FCC023", 
            "#0600FB","#FE81C8","#C78AAB","#C75E7D","#FF63B9")
cbp <- c(corlors, "#F19BFF", "#808080")

# ------------------------------------------------------------------------------
# TSNE plot 
fit_s <- readRDS(paste0("hippoAll_Manual_icmem_", embro, ".rds"))
Idents(slide.seq) <- factor(slide.seq@meta.data[, "icmemlfc"], levels=types)

seu <- slide.seq
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
pdf("hippo_tSNE.pdf")
p1 <- Seurat::DimPlot(seu, cells=which(slide.seq@meta.data[, "icmemlfc"] != "Unknown"), reduction = 'tSNE_spanno', pt.size = 1, cols= cbp) + 
        ggtitle("spatialAnno") + 
           theme(plot.title = element_text(face = "bold", size = rel(2), hjust = 0.5),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.title = element_text(size=25),
                axis.line = element_line(colour="black"),
                legend.text = element_text(size=25),
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                 legend.position = "none",
                legend.title = element_blank())
print(p1)
dev.off()


# pca
seu <- SCTransform(seu, assay = "RNA", return.only.var.genes = FALSE, verbose = FALSE)
seu <- RunPCA(seu)
seu <- RunTSNE(seu, reduction="pca", dims=1:15, verbose = F, check_duplicates = FALSE)
pc <- seu@reductions$pca@cell.embeddings

# drsc
drsc <- readRDS(paste0("hippo_", embro, "_DR_SC.rds"))
Ez <- drsc$hZ
row.names(Ez) <- colnames(seu)
colnames(Ez) <- paste0('drsc_', 1:ncol(Ez))
seu@reductions$"drsc" <- CreateDimReducObject(embeddings = Ez, key='drsc_', assay=DefaultAssay(seu))
seu <- RunTSNE(seu, reduction="drsc", dims=1:15, verbose = F, check_duplicates = FALSE)



rotate_angles <- c(60, 93)[embro]
library(scales)
coordinate_rotate <- function(pos, theta=0){# counter-clock rotation
  pos_new <- pos
  pos_new[,1] <- pos[,1]*cos(theta) - pos[,2]*sin(theta)
  pos_new[,2] <- pos[,1]*sin(theta) + pos[,2]*cos(theta)
  return(pos_new)
}  

pos <- cbind(slide.seq$row, slide.seq$col)
pos_new <- coordinate_rotate(pos, theta=rotate_angles)
#loc1 <- pos_new[, 1]; loc2 <- pos_new[, 2]
#loc1 <- slide.seq$row; loc2 <- slide.seq$col

coord.df = data.frame(x=pos_new[, 1], y=pos_new[, 2], stringsAsFactors=FALSE) # (stringsAsFactors only if also have a separate barcodes column)
rownames(coord.df) = attr(slide.seq@active.ident, "names")
slide.seq@images$image =  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = coord.df
  )

######################### Supplementary Figure 21a or 2b########################
# RGB plot
library(SpatialPCA)
set.seed(1234)
xy_coords <- data.frame(
        x_coord = pos_new[, 2],#seu@meta.data[, c("col")], 
        y_coord = -pos_new[, 1]#-seu@meta.data[, c("row")]
    )
p1 = plot_RGB_tSNE(xy_coords, t(hZ), pointsize=0.5)[[2]] + theme(plot.title = element_blank())
p2 = plot_RGB_tSNE(xy_coords, t(pc), pointsize=0.5)[[2]] + theme(plot.title = element_blank())
p3 = plot_RGB_tSNE(xy_coords, t(Ez), pointsize=0.5)[[2]] + theme(plot.title = element_blank())
pp <- p1  + p3  + p2  + plot_layout(byrow = T, nrow=1, ncol=3)  
ggsave(file=paste0('Hippo_', embro, 'RGB_TSNE.pdf'), plot = pp, title="",
       width = 24, height =8)

if(FALSE) { 
  p1 = plot_RGB_UMAP(xy_coords, t(hZ), pointsize=3)[[2]] + theme(plot.title = element_blank())
  p2 = plot_RGB_UMAP(xy_coords, t(pc), pointsize=3)[[2]] + theme(plot.title = element_blank())
  p3 = plot_RGB_UMAP(xy_coords, t(Ez), pointsize=3)[[2]] + theme(plot.title = element_blank())
  pp <- p1  + p3  + p2  + plot_layout(byrow = T, nrow=1, ncol=3)  
  ggsave(file=paste0('Hippo_', embro, 'RGB_UMAP.pdf'), plot = pp, title="",
       width = 24, height =8)
}

# ------------------------------------------------------------------------------
######################## Figure 4b ########################
 
K <- length(table(slide.seq@meta.data[, "icmemlfc"]))

#cbp <- c(hue_pal()(K-1), "#808080")


pos_gg <- coordinate_rotate(pos, theta=c(90, 60)[embro])
loc1 <- pos_gg[, 1]
loc2 <- pos_gg[, 2]

# ground truth
cluster <- factor(slide.seq@meta.data[, "icmemlfc"], levels=types)
 slide.seq@meta.data[1, "time_icmemlfc"]/3600 

      datt = data.frame(cluster, loc1, loc2); #datt <- datt[datt$cluster=='Oligodendrocyte', ]
      p0 = ggplot(datt, aes(x = loc1, y = loc2, color = cluster)) +
                geom_point( alpha =1, size=0.45) +
                scale_color_manual(values = cbp, drop=FALSE)
pdf(paste0("hippoAll_", ifelse(dataSource==1, "satija", "zhou"), mk_n[iter], "_", embro, "_icmem.pdf"))
p <- p0 +  theme_void()+
                theme(plot.title = element_blank(),
                     plot.margin = margin(-0.1, -0.1, -0.1, -0.1, unit = "pt"),
              legend.position = "none"
              )
     print(p)
dev.off() 
               

library(grid)
# plot legend
pdf(file = "hippo_ground_truth_legend.pdf",width = 5, height = 4)
p <- p0 + theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5)))
legend <- cowplot::get_legend(p)  
grid.newpage()
grid.draw(legend)
dev.off()


cluster = factor(slide.seq@meta.data[, "layer_scina"], levels=types) 
 slide.seq@meta.data[1, "time_scina"]/3600 
pdf(paste0("hippoAll_", ifelse(dataSource==1, "satija", "zhou"), mk_n[iter], "_", embro, "_scina.pdf"))
      datt = data.frame(cluster, loc1, loc2); #datt <- datt[datt$cluster=='Oligodendrocyte', ]
      p = ggplot(datt, aes(x = loc1, y = loc2, color = cluster)) +
                geom_point( alpha =1, size=0.45) +
                scale_color_manual(values = cbp, drop=FALSE)+
                theme_void()+
                theme(plot.title = element_blank(),
                     plot.margin = margin(-0.1, -0.1, -0.1, -0.1, unit = "pt"),
              legend.position = "none"
              ) #+  scale_x_reverse()
    print(p)
dev.off() 

#rts <- readRDS(paste0("hippoAll_Manual_scSorter_", mk_n[iter], "_", embro, ".rds")) 
#slide.seq@meta.data["layer_scSorter"] <- rts$Pred_Type

cluster = factor(slide.seq@meta.data[, "layer_scSorter"], levels=types)
 slide.seq@meta.data[1, "time_scSorter"]/3600 
pdf(paste0("hippoAll_", ifelse(dataSource==1, "satija", "zhou"), mk_n[iter], "_", embro, "_scSorter.pdf"))
      datt = data.frame(cluster, loc1, loc2); #datt <- datt[datt$cluster=='Oligodendrocyte', ]
      p = ggplot(datt, aes(x = loc1, y = loc2, color = cluster)) +
                geom_point( alpha =1, size=0.45) +
                scale_color_manual(values = cbp, drop=FALSE)+
                theme_void()+
                theme(plot.title = element_blank(),
                     plot.margin = margin(-0.1, -0.1, -0.1, -0.1, unit = "pt"),
              legend.position = "none"
              ) #+  scale_x_reverse()
    print(p)
dev.off() 


cluster = factor(slide.seq@meta.data[, "layer_cellassign"], levels=types)
slide.seq@meta.data[1, "time_cellassign"]/3600 
pdf(paste0("hippoAll_", ifelse(dataSource==1, "satija", "zhou"), mk_n[iter], "_", embro, "_cellAssign.pdf"))
      datt = data.frame(cluster, loc1, loc2); #datt <- datt[datt$cluster=='Oligodendrocyte', ]
      p = ggplot(datt, aes(x = loc1, y = loc2, color = cluster)) +
                geom_point( alpha =1, size=0.45) +
                scale_color_manual(values = cbp, drop=FALSE)+
                theme_void()+
                theme(plot.title = element_blank(),
                     plot.margin = margin(-0.1, -0.1, -0.1, -0.1, unit = "pt"),
              legend.position = "none"
              ) #+  scale_x_reverse()
    print(p)
dev.off() 
 

# garnett
garnett <- readRDS(paste0("hippo", embro, "_anno_garnett.rds"))
cluster <- factor(garnett, levels=types)
pdf(paste0("hippoAll_", ifelse(dataSource==1, "satija", "zhou"), mk_n[iter], "_", embro, "_garnett.pdf"))
      datt = data.frame(cluster, loc1, loc2); #datt <- datt[datt$cluster=='Oligodendrocyte', ]
      p = ggplot(datt, aes(x = loc1, y = loc2, color = cluster)) +
                geom_point( alpha =1, size=0.45) +
                scale_color_manual(values = cbp, drop=FALSE)+
                theme_void()+
                theme(plot.title = element_blank(),
              legend.position = "none"
              ) #+  scale_x_reverse()
    print(p)
dev.off() 

# clustering results of unsupervised methods.
cbp_spatialpca = c( "#FD7446" ,"#709AE1", "#31A354","#9EDAE5",
 "#DE9ED6" ,"#BCBD22", "#CE6DBD" ,"#DADAEB" ,
 "yellow", "#FF9896","#91D1C2", "#C7E9C0" ,
 "#00AFBB", "#E7B800", "#FC4E07", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7",
 "#66C2A5", "lightyellow2", #"cornflowerblue" ,"#E78AC3", "skyblue1", 
 "#6B6ECF", "#7B4173" )

pdf(paste0("hippo_", embro, "_bsCluster.pdf"))
cluster <- factor(readRDS(paste0("hippo_", embro, "_Unsupervised_cluster_bs.rds")), levels=1:25)
       datt = data.frame(cluster, loc1, loc2);  
      p = ggplot(datt, aes(x = loc1, y = loc2, color = cluster)) +
                geom_point(alpha =1, size=0.45 ) + 
                theme_void()+
                scale_color_manual(values = cbp_spatialpca)+
                theme(plot.title = element_blank(),
                 text = element_text(size = 35),
                  legend.key.size= unit(0.6, "cm"),
                  legend.title = element_blank(),
                 legend.position = "bottom"
              )  + guides(color=guide_legend(override.aes = list(size = 5)))  
    print(p)
dev.off() 

pdf(paste0("hippo_", embro, "_drscCluster.pdf"))
cluster <- factor(readRDS(paste0("hippo_", embro, "_DR_SC.rds"))$cluster, levels=1:25)
       datt = data.frame(cluster, loc1, loc2);  
      p = ggplot(datt, aes(x = loc1, y = loc2, color = cluster)) +
                geom_point(alpha =1, size=0.45 ) + 
                theme_void()+
                scale_color_manual(values = cbp_spatialpca)+
                theme(plot.title = element_blank(),
                 text = element_text(size = 35),
                 legend.key.size= unit(0.6, "cm"),
                  legend.title = element_blank(),
                 legend.position = "bottom"
              ) + guides(color=guide_legend(override.aes = list(size = 5)))  
    print(p) 
dev.off()  

# -------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------
# heatmap of markers genes
######################### Figure 4c or 4g#########################
 
seu <- slide.seq

featurePlot <- function(seu, feature=NULL, cols=NULL, pt_size=0.45){
  datt <- data.frame(loc1, loc2)

  if(is.null(feature)) feature <- row.names(seu)[1]
  datt$Expression <- log(1+seu[["RNA"]]@counts[feature,])
  alpha <- rep(1, length( datt$Expression )); alpha[seu[["RNA"]]@counts[feature,]==0] = 0.1  

  g <- ggplot(datt, aes(x = loc1, y = loc2, color = Expression)) + geom_point(alpha =alpha, size=pt_size) +
     scale_colour_gradient(low ="#bebebe", high = "#DE2D26") +  
     theme_void()+ theme(plot.title = element_blank(), legend.position = "none")
  g
}

pdf("Wfs1.pdf")
SpatialFeaturePlot(slide.seq, feature="Wfs1", alpha = c(0.5, 2)) + 
 theme(plot.title = element_blank(), legend.position = "none", plot.margin = margin(0, 0, 0, 0, unit = "pt"))
dev.off()


#  markers$'CA1 Principal cells' <- c("WFS1", "FIBCD1", "ATP2B1", "ITPKA", "PPP3CA") 
# markers$'CA3 Principal cells' <- c("CHGB", "HS3ST4", "NPTXR", "CPNE4", "NEUROD6") 
#  markers$'Dentate Principal cells' <- c("C1QL2", "PPP3CA", "FAM163B", "OLFM1", "NCDN")
pdf(paste0(embro, "_CA1.pdf"))
p0 <- featurePlot(slide.seq, feature="Wfs1") #+ theme(legend.position = "right")
p0
dev.off()

pdf(file = paste0(embro,"_GE_legend.pdf"), width = 5, height = 4)
p <- p0 + theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        #legend.position="bottom",
        legend.title = element_blank(),
        axis.ticks = element_blank()) 
legend <- cowplot::get_legend(p)  
grid.newpage()
grid.draw(legend)
dev.off()


pdf(paste0(embro, "_CA3.pdf"))
featurePlot(slide.seq, feature="Cpne4") #+ theme(legend.position = "right")
dev.off()

pdf(paste0(embro, "_Dentate.pdf"))
featurePlot(slide.seq, feature="C1ql2") #+ theme(legend.position = "right")
dev.off()
# -------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------
# spatial distribution of cells
Idents(slide.seq) <- slide.seq@meta.data[, "icmemlfc"] 
pt_size=0.45
palette <- c("#0173b2", "#de8f05", "#029e73", "#d55e00", "#cc78bc",
             "#ca9161", "#fbafe4", "#949494", "#ece133", "#56b4e9")[c(1, 8)]
DimPlot <- function(seu, highligth, plot_title=NULL)    {
   datt <- data.frame(loc1, loc2)
    ll <- Idents(slide.seq) 
    levels(ll)[levels(ll) != highligth] <- "New"  
    datt$label <- factor(ll, levels=c(highligth, "New")) 
    #datt$label <-  factor((Idents(slide.seq)  == highligth)+1, levels=c("a", "b"))

    g <- ggplot(datt, aes(x = loc1, y = loc2, color = label)) + geom_point(alpha =1, size=pt_size) +
      theme_void() + theme(legend.position = "none") + scale_color_manual(values=palette)
    if(!is.null(plot_title)) g <- g + labs(title=plot_title) + theme(plot.title = element_text(face = "bold", size = rel(3), hjust = 0.5))
    if( is.null(plot_title)) g <- g + theme(plot.title =  element_blank())
    g 
}


pdf(paste0(embro, "_CA1_cell.pdf"))
DimPlot(slide.seq, highligth = "CA1 Principal cells")
dev.off()
 
pdf(paste0(embro, "_CA3_cell.pdf"))
DimPlot(slide.seq, highligth = "CA3 Principal cells")
dev.off()

pdf(paste0(embro, "_Dentate_cell.pdf"))
DimPlot(slide.seq, highligth = "Dentate Principal cells")
dev.off()


slide.seq@meta.data[, "icmemlfc"][slide.seq@meta.data[, "icmemlfc"]=="TH"] = "Hb neuron"
slide.seq@meta.data[, "layer_scina"][slide.seq@meta.data[, "layer_scina"]=="TH"] = "Hb neuron"
slide.seq@meta.data[, "layer_scSorter"][slide.seq@meta.data[, "layer_scSorter"]=="TH"] = "Hb neuron"
slide.seq@meta.data[, "layer_cellassign"][slide.seq@meta.data[, "layer_cellassign"]=="TH"] = "Hb neuron"
slide.seq@meta.data$garnett = garnett
slide.seq@meta.data[, "garnett"][slide.seq@meta.data[, "garnett"]=="TH"] = "Hb neuron"
################################### Figure 4e ##########################################
#######  technique difference between two hippo1 and hippo2
seu1 <- readRDS(paste0("hippoAll_anno_", ifelse(dataSource==1, "satija", "zhou"), mk_n[iter], "_", 1, ".rds"))
seu2 <- readRDS(paste0("hippoAll_anno_", ifelse(dataSource==1, "satija", "zhou"), mk_n[iter], "_", 2, ".rds"))

#feature <- intersect(rownames(seu1), rownames(seu2))
#V1 <- log(rowSums(seu1[["RNA"]]@counts[feature,]))
#V2 <- log(rowSums(seu2[["RNA"]]@counts[feature,]))
#umi <- data.frame(V1=V1, V2=V2) 

seuList <- readRDS(file='slideV2_mouseHip2_seu.RDS')
V1 <- data.frame(pf = "Slide-seq",   Counts=log10(seuList[[1]]@meta.data$nFeature_RNA))
V2 <- data.frame(pf = "Slide-seqV2", Counts=log10(seuList[[2]]@meta.data$nFeature_RNA))
umi <- rbind(V1=V1, V2=V2) 
umi$pf <- as.factor(umi$pf)

 
a <- ggplot(umi, aes(x = Counts, color = pf, fill=pf))  
pdf("slide_seq.pdf")
gg <-   a + geom_histogram(bins = 120, alpha = 0.3, position = "identity") +
        #a +     geom_boxplot() 
            scale_fill_manual(values = c("#E7B800", "#00AFBB")) +
            scale_color_manual(values = c( "#E7B800", "#00AFBB")) + 
            ylab("Number of beads") +
            xlab("Log10 Number of UMIs") + 
            theme(plot.title = element_blank(),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 text = element_text(size = 35),
                 axis.title.x = element_blank(),
                 legend.key.size= unit(0.6, "cm"),
                 legend.title = element_blank(),
                 legend.position = "bottom"
            )  
print(gg)
dev.off()
  
##################### extract legend for Figure 4b############################
pi1 <- table(factor(seu1@meta.data[, method], levels=types))/ncol(seu1)
pi2 <- table(factor(seu2@meta.data[, method], levels=types))/ncol(seu2)
dfAll <- data.frame(pi1, pi2);dfAll$Var1.1 <- NULL 
colnames(dfAll) <- c("types", "V1", "V2")
p <- ggplot(dfAll, aes(V1, V2, color = types)) +  scale_color_manual(values = cbp)  + 
     geom_point(size=2) + 
     theme( #plot.title = element_blank(),
            #text = element_text(size=40),
            #strip.text = element_text(face="bold"), 
            #plot.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "pt"),
            legend.key.size= unit(0.4, "cm"),
            #legend.title = element_blank(),
            legend.position = "right") + theme_void() + labs(col="") +
            guides(color=guide_legend(override.aes = list(size = 5))) 
  

leg <- get_legend(p)

# Convert to a ggplot and print
pdf("hippo_legend.pdf")
grid.newpage()                              
grid.draw(leg) 
dev.off()


############################## Figure 4f ################################
####### composition difference between two hippo1 and hippo2
dfAll$cbp <- cbp
df <- dfAll[-22, ]
saveRDS(df, "prop.rds")
#  setwd("~/Cloud/Research/SSL/code/data/Slide-seqV2/") # on local laptop
df <- readRDS("prop.rds")
cor(df$V1, df$V2)

fit <- summary(lm(V2~V1, data=df))
R2 <- round(fit$r.squared, 2)
slope <- round(fit$ coefficients[2,1], 2)
pv <- formatC(pf(fit$fstatistic[1],fit$fstatistic[2],fit$fstatistic[3],lower.tail=FALSE), format = "e", digits = 2)
pdf("Slide-seqFreq.pdf")
p <- ggplot(df,aes(V1, V2)) +
  geom_point(color=df$cbp, size=5) +
  geom_smooth(method='lm', color=1, lwd=2) + 
  geom_abline(slope = 1, intercept=0, lty=2, lwd=2) +
   labs(x='Slide-seq', y='Slide-seqV2', title=paste0('slope=', slope, " (", pv, ")")) 
p +   theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
            text = element_text(size=25),
            panel.background = element_blank(),
            plot.background = element_rect(colour = NA),
            #axis.title = element_blank(),
            #panel.grid.minor = element_blank(),
            axis.line = element_line(colour="black"),
            #plot.margin=unit(c(5,7,5,0),"mm"),
            strip.text = element_text(face="bold"),
            legend.position = "none"
            )  
dev.off()  


 