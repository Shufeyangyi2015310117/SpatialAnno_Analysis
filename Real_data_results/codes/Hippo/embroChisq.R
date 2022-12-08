########################## Figure 4d #########################

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
set.seed(1234)

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


coord.df = data.frame(x=pos_new[, 1], y=pos_new[, 2], stringsAsFactors=FALSE) # (stringsAsFactors only if also have a separate barcodes column)
rownames(coord.df) = attr(slide.seq@active.ident, "names")
slide.seq@images$image =  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = coord.df
  )


garnett <- readRDS(paste0("hippo", embro, "_anno_garnett.rds"))
slide.seq@meta.data$garnett = garnett
slide.seq@meta.data[, "garnett"][slide.seq@meta.data[, "garnett"]=="TH"] = "Hb neuron"


# loc: position
# x: a matrix to be smoothed, with its neighbors
smoothX <- function(Adj_sp, x)
{ 
    #x <- dat_hl[, -(1:3)]
    n <- nrow(x)
    for (i in 1:n) {
        #x0 <- x[i,]
        i_nb <- which(Adj_sp[i, ] !=0)
        if (length(i_nb) >= 2)  {
            x_nb <- x[i_nb, ]
            x[i,] <- colMeans(x_nb)
        }
    }
    x
}
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
smoothY <- function(Adj_sp, y)
{
    #y <- dat_hl$y
    n <- length(y)
    for (i in 1:n) {
        y0 <- y[i]
        i_nb <- which(Adj_sp[i, ] !=0)
        if (length(i_nb) >= 2)  {
            y_mode <- Mode(y[i_nb])
            if(y0 != y_mode) y[i] <- y_mode
        }
    }
    y
}


####### 
highlights = c("CA1 Principal cells", "CA3 Principal cells", "Dentate Principal cells")

methods <- c("icmemlfc", "layer_scina", "layer_scSorter", "garnett", "layer_cellassign") 
j <- 1
method <- methods[j]
seu <- slide.seq
Y <- cbind(seu@meta.data[, methods])
hipR <- (Y$icmemlfc%in%highlights) + (Y$layer_scina%in%highlights)  + (Y$layer_scSorter%in%highlights) + (Y$layer_cellassign%in%highlights) + (Y$garnett%in%highlights)
beads_hl <- rownames(Y)[hipR != 0]
 markers <- list()
  markers$'CA1 Principal cells' <- c("WFS1", "FIBCD1", "ATP2B1", "ITPKA", "PPP3CA")#[1] 
  markers$'CA3 Principal cells' <- c("CHGB", "HS3ST4", "NPTXR", "CPNE4", "NEUROD6")#[4]
  markers$'Dentate Principal cells' <- c("C1QL2", "PPP3CA", "FAM163B", "OLFM1", "NCDN")#[1]
R <- matrix(0, nrow = 5, ncol = 3); rownames(R) <- methods
Fs <- matrix(0, nrow = 5, ncol = 3)
Chi <- matrix(0, nrow = 5, ncol = 3)
Odd <- matrix(0, nrow = 5, ncol = 3)
Ncell <- matrix(0, nrow = 5, ncol = 3)
Nb <- matrix(0, nrow = 5, ncol = 3)
library(SpatialPack)
Adj <- SC.MEB:::getneighborhood_fast(cbind(seu$row, seu$col), cutoff = 30)
mg <- c(1, 4, 1)
#library(pROC)
for(k in 1:3) {    
    feature <- markers[[k]]
    X <- t(as.matrix(seu[["RNA"]]@counts[match(feature, toupper(rownames(seu))), ]))
    xsmooth <- smoothX(Adj, X)

    for(m in 1:5)   {
        highlight <- highlights[k]
        y <- 1 * (Y[, m] == highlight)
        ysmooth <- smoothY(Adj, y)

        Ncell[m, k] <- sum(y)
        dat <- data.frame(coords1 = seu$row, coords2 = seu$col, y=ysmooth, X=xsmooth)
        dat_hl <- dat
        coords <- cbind(dat_hl$coords1, dat_hl$coords2)
        ctj <- table(1*(dat_hl[,3+mg[k]]!=0), dat_hl$y)
        Chi[m, k] <- chisq.test(1*(dat_hl[,3+mg[k]]!=0), dat_hl$y)$statistic #estimate
        Odd[m, k] <- fisher.test(ctj)$estimate
        mcnemar.test(ctj)
        dd = dat_hl[dat_hl$y==1, (1:2)]
        Adj_sp <- SC.MEB:::getneighborhood_fast(cbind(dd[,1], dd[,2]), cutoff = 30)
        Nb[m, k] <- mean(rowSums(Adj_sp))
    }
}

rownames(Chi) <- c("SpatialAnno", "SCINA", "scSorter", "Garnett", "CellAssign")
colnames(Chi) <- c("CA1", "CA3", "DG")

rownames(Odd) <- c("SpatialAnno", "SCINA", "scSorter", "Garnett", "CellAssign")
colnames(Odd) <- c("CA1", "CA3", "DG")

saveRDS(list(Chi=Chi, Odd=Odd), paste0("embro", "Chi", embro, ".rds"))

