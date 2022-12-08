################## Supplementary Figure 18-20, 22-24

setwd("./Real_data_results/dataFiles/Hippo/")
library(ggplot2)
library(patchwork)
library(dplyr)
library(grid)
library(Seurat)
embro <- 2 # 26177cells:36h, 43994 cells:48h 
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

 
# ------------------------------------------------------------------------------
# heatmap  
 
 

pos_gg <- coordinate_rotate(pos, theta=c(90, 60)[embro])
loc1 <- pos_gg[, 1]
loc2 <- pos_gg[, 2]
  
# garnett
garnett <- readRDS(paste0("hippo", embro, "_anno_garnett.rds"))
 
# -------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------
# heatmap of markers genes

  
# -------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------
# spatial distribution of cells
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
 
slide.seq@meta.data[, "icmemlfc"][slide.seq@meta.data[, "icmemlfc"]=="TH"] = "Hb neuron"
slide.seq@meta.data[, "layer_scina"][slide.seq@meta.data[, "layer_scina"]=="TH"] = "Hb neuron"
slide.seq@meta.data[, "layer_scSorter"][slide.seq@meta.data[, "layer_scSorter"]=="TH"] = "Hb neuron"
slide.seq@meta.data[, "layer_cellassign"][slide.seq@meta.data[, "layer_cellassign"]=="TH"] = "Hb neuron"
slide.seq@meta.data$garnett = garnett
slide.seq@meta.data[, "garnett"][slide.seq@meta.data[, "garnett"]=="TH"] = "Hb neuron"


######### plot different methods together
celltype <- sort(unique(slide.seq@meta.data[, "icmemlfc"]))
meth <- c("icmemlfc", "layer_scina", "layer_scSorter", "garnett", "layer_cellassign")
Mname <- c("SpatialAnno", "SCINA", "scSorter", "Garnett", "CellAssign")
M <- length(meth)
plist <- list()# 1:8, 9:16, 17:22
for(k in 1:8){  
    for (m in 1:M) {
        Idents(slide.seq) <- slide.seq@meta.data[, meth[m]] 
        pg <- DimPlot(slide.seq, highligth = celltype[[k]], plot_title=celltype[[k]])

        if (k==1) pg <- pg + labs(title=Mname[m])
        if (k!=1) pg <- pg + labs(title=element_blank()) 
        if (m==1) pg <- pg +
          scale_x_continuous(expand = expansion(mult = c(0.05, 0.01))) +
          annotate("text", min(loc1)-300, median(loc2), label = celltype[[k]], angle='90', size = 11)  # 10 mm
        plist[[(k-1)*5+m]] <- pg
    }
} 
library(gridExtra) # also loads grid
png(filename = paste0(embro, "_spaDistr_1.png"), width = 800, height =1280, units = "mm", res=200, type = "cairo")
grid.arrange(grobs=plist, nrow = 8)
dev.off()

plist <- list()# 1:8, 9:16, 17:22
for(k in 1:8){  
    for (m in 1:M) {
        Idents(slide.seq) <- slide.seq@meta.data[, meth[m]] 
        pg <- DimPlot(slide.seq, highligth = celltype[[k+8]], plot_title=celltype[[k+8]])

        if (k==1) pg <- pg + labs(title=Mname[m])
        if (k!=1) pg <- pg + labs(title=element_blank()) 
        if (m==1) pg <- pg +
          scale_x_continuous(expand = expansion(mult = c(0.05, 0.01))) +
          annotate("text", min(loc1)-300, median(loc2), label = celltype[[k+8]], angle='90', size = 11)  # 10 mm
        plist[[(k-1)*5+m]] <- pg
    }
} 
library(gridExtra) # also loads grid
png(filename = paste0(embro, "_spaDistr_2.png"), width = 800, height =1280, units = "mm", res=200, type = "cairo")
grid.arrange(grobs=plist, nrow = 8)
dev.off()


plist <- list()# 1:8, 9:16, 17:22
for(k in 1:6){  
    for (m in 1:M) {
        Idents(slide.seq) <- slide.seq@meta.data[, meth[m]] 
        pg <- DimPlot(slide.seq, highligth = celltype[[k+16]], plot_title=celltype[[k+16]])

        if (k==1) pg <- pg + labs(title=Mname[m])
        if (k!=1) pg <- pg + labs(title=element_blank())
        if (m==1) pg <- pg +
          scale_x_continuous(expand = expansion(mult = c(0.05, 0.01))) + 
          annotate("text", min(loc1)-300, median(loc2), label = celltype[[k+16]], angle='90', size = 11)  # 10 mm
        plist[[(k-1)*5+m]] <- pg
    }
} 
library(gridExtra) # also loads grid
png(filename = paste0(embro, "_spaDistr_3.png"), width = 800, height =960, units = "mm", res=200, type = "cairo")
grid.arrange(grobs=plist, nrow = 6)
dev.off()
 
  