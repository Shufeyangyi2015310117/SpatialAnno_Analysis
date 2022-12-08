
##########################################
##########################################
# prepare input for Brain12plot.R
##########################################  
##########################################
# 1-------------------------------------------------------------------------
# results from manually selected markers.
library(SingleCellExperiment) 
library(psych)

name_ID12 <- as.character(c(151507, 151508, 151509, 151510, 151669, 151670,
                              151671, 151672, 151673, 151674, 151675, 151676))
samp <- rep(1:3, each=4)
 

dir.exp <- "./Real_data_results/dataFiles/DLPFC/Brain12cross/"

# computational times
tim <- matrix(0, nrow=12, ncol=4)
for (iter in (1:12)) {
    ID <- name_ID12[iter] 
    dlpfc <- readRDS(paste0(dir.exp, ID, "SCINA_Garnett.rds")) 
    top <- 5
    tim[iter, 1] <- colData(dlpfc)[paste0("time_icmemlfc", id, "top", top)][1,]
    tim[iter, 2] <- colData(dlpfc)[paste0("layer_scSorter", id, "top", top)][1,]
    tim[iter, 3] <- colData(dlpfc)[paste0("time_scina", id, "top", top)][1,]
    tim[iter, 4] <- colData(dlpfc)[paste0("time_cellassign", id, "top", top)][1,]           
}

write.table(round(tim, 3), "time.txt",sep="\t")

# -------------------------------------------------------------
# 1.1 over specified cell types for the 2nd sample.
kp <- matrix(0, nrow = 12, ncol = 5)
f1 <- matrix(0, nrow = 12, ncol = 5)
ac <- matrix(0, nrow = 12, ncol = 5)

for (iter in (5:8)) {
    ID <- name_ID12[iter] 
    dlpfc <- readRDS(paste0(dir.exp, ID, "Manual.rds")) 
    y <- as.character(dlpfc$layer_guess_reordered)
    y[is.na(y)] <- 'Unknown'
    layer_sp <- colData(dlpfc)["icmemlfc"][[1]]
    res_sp <- evaluate(y, layer_sp, Indices = NULL)
    kp[iter, 1] <- cohen.kappa(x=cbind(y, layer_sp))$kappa
    f1[iter, 1] <- mean(res_sp$F1)
    ac[iter, 1] <- mean(y == layer_sp)  

    layer_sc <- colData(dlpfc)["layer_scSorter"][[1]] 
    res_sc <- evaluate(y, layer_sc, Indices = NULL)
    kp[iter, 2] <- cohen.kappa(x=cbind(y, layer_sc))$kappa
    f1[iter, 2] <- mean(res_sc$F1)
    ac[iter, 2] <- mean(y == layer_sc)

    layer_scina <- colData(dlpfc)["layer_scina"][[1]] 
    layer_scina[layer_scina == "unknown"] <- 'Unknown'
    res_scina <- evaluate(y, layer_scina, Indices = NULL)
    kp[iter, 3] <- cohen.kappa(x=cbind(y, layer_scina))$kappa
    f1[iter, 3] <- mean(res_scina$F1)
    ac[iter, 3] <- mean(y == layer_scina)

    layer_cellassign <- colData(dlpfc)["layer_cellassign"][[1]] 
    res_cellassign <- evaluate(y, layer_cellassign, Indices = NULL)
    kp[iter, 4] <- cohen.kappa(x=cbind(y, layer_cellassign))$kappa
    f1[iter, 4] <- mean(res_cellassign$F1)
    ac[iter, 4] <- mean(y == layer_cellassign)

    # at beginning, garnett was run with overspecified.
    dir.exp1 <- "./Real_data_results/dataFiles/DLPFC/Brain12cross/"
    dlpfcSCINA <- readRDS(paste0(dir.exp1, ID, "SCINA_Garnett.rds") )
    layer_garnett <- colData(dlpfcSCINA)[paste0("garnett", "151507","top", 5)][[1]] 
    res_garnett <- evaluate(y, layer_garnett, Indices = NULL)
    kp[iter, 5] <- cohen.kappa(x=cbind(y, layer_garnett))$kappa
    f1[iter, 5] <- mean(res_garnett$F1)
    ac[iter, 5] <- mean(y == layer_garnett)
}
colnames(kp) <- c("SpatialAnno", "scSorter", "SCINA", "CellAssign", "Garnett")
colnames(ac) <- c("SpatialAnno", "scSorter", "SCINA", "CellAssign", "Garnett")
colnames(f1) <- c("SpatialAnno", "scSorter", "SCINA", "CellAssign", "Garnett")
rownames(kp) <- name_ID12; rownames(ac) <- name_ID12; rownames(f1) <- name_ID12; 

saveRDS(list(kp=kp[5:8,], ac=ac[5:8,], f1=f1[5:8,]), file="dlpfc_sample2_OverSpecify.rds")

# 1.2 correctly specified cell types for the 2nd sample.
kpC <- matrix(0, nrow = 12, ncol = 5)
f1C <- matrix(0, nrow = 12, ncol = 5)
acC <- matrix(0, nrow = 12, ncol = 5)
id <- '151507'
top <- 5
for (iter in (5:8)) {
    ID <- name_ID12[iter] 
   
    dir.exp <- "./Real_data_results/dataFiles/DLPFC/Brain12cross/"
    dlpfc <- readRDS(paste0(dir.exp, ID, "SCINA_Garnett.rds") )

    y <- as.character(dlpfc$layer_guess_reordered)
    y[is.na(y)] <- 'Unknown'
    layer_sp <- colData(dlpfc)[paste0("icmemlfc", id, "top", top)][[1]]
    res_sp <- evaluate(y, layer_sp, Indices = NULL)
    kpC[iter, 1] <- cohen.kappa(x=cbind(y, layer_sp))$kappa
    f1C[iter, 1] <- mean(res_sp$F1)
    acC[iter, 1] <- mean(y == layer_sp)  

    layer_sc <- colData(dlpfc)[paste0("layer_scSorter", id,"top", top)][[1]] 
    res_sc <- evaluate(y, layer_sc, Indices = NULL)
    kpC[iter, 2] <- cohen.kappa(x=cbind(y, layer_sc))$kappa
    f1C[iter, 2] <- mean(res_sc$F1)
    acC[iter, 2] <- mean(y == layer_sc)

    layer_scina <- colData(dlpfc)[paste0("layer_scina", id,"top", top)][[1]] 
    layer_scina[layer_scina == "unknown"] <- 'Unknown'
    res_scina <- evaluate(y, layer_scina, Indices = NULL)
    kpC[iter, 3] <- cohen.kappa(x=cbind(y, layer_scina))$kappa
    f1C[iter, 3] <- mean(res_scina$F1)
    acC[iter, 3] <- mean(y == layer_scina)

    layer_cellassign <- colData(dlpfc)[paste0("layer_cellassign", id,"top", top)][[1]] 
    res_cellassign <- evaluate(y, layer_cellassign, Indices = NULL)
    kpC[iter, 4] <- cohen.kappa(x=cbind(y, layer_cellassign))$kappa
    f1C[iter, 4] <- mean(res_cellassign$F1)
    acC[iter, 4] <- mean(y == layer_cellassign)

    # at beginning, garnett was run with overspecified.
    dlpfcC <- readRDS(paste0(dir.exp, ID, "CorrectSpecify_Garnett.rds")) 
    layer_garnett <- colData(dlpfcC)[paste0("garnett", id, "top", top)][[1]]
    res_garnett <- evaluate(y, layer_garnett, Indices = NULL)
    kpC[iter, 5] <- cohen.kappa(x=cbind(y, layer_garnett))$kappa
    f1C[iter, 5] <- mean(res_garnett$F1)
    acC[iter, 5] <- mean(y == layer_garnett)
}
colnames(kpC) <- c("SpatialAnno", "scSorter", "SCINA", "CellAssign", "Garnett")
colnames(acC) <- c("SpatialAnno", "scSorter", "SCINA", "CellAssign", "Garnett")
colnames(f1C) <- c("SpatialAnno", "scSorter", "SCINA", "CellAssign", "Garnett")
rownames(kpC) <- name_ID12; rownames(acC) <- name_ID12; rownames(f1C) <- name_ID12; 
saveRDS(list(kpC=kpC[5:8,], acC=acC[5:8,], f1C=f1C[5:8,]), file="dlpfc_sample2_CorrectSpecify.rds")


# 2-------------------------------------------------------------------------
# semi-supervised resluts based on markers from DeG of other Visium slices.
setwd("./Real_data_results/dataFiles/DLPFC/Brain12cross")
library(SingleCellExperiment) 
library(psych)

name_ID12 <- as.character(c(151507, 151508, 151509, 151510, 
                            151669, 151670, 151671, 151672, 
                            151673, 151674, 151675, 151676))
samp <- rep(1:3, each=4)
 
#topMarker <- c(10, 5, 5, 5, 
                #5, 10, 10, 5, 
                #5, 5, 5, 5)
dir.exp <- "./Real_data_results/dataFiles/DLPFC/Brain12cross/"
p <- 24
kp_sp <- matrix(0, nrow = 12, ncol = p)
f1_sp <- matrix(0, nrow = 12, ncol = p)
ac_sp <- matrix(0, nrow = 12, ncol = p)

kp_icmemlfc0 <- matrix(0, nrow = 12, ncol = p)
f1_icmemlfc0 <- matrix(0, nrow = 12, ncol = p)
ac_icmemlfc0 <- matrix(0, nrow = 12, ncol = p)

kp_icmemlfc2 <- matrix(0, nrow = 12, ncol = p)
f1_icmemlfc2 <- matrix(0, nrow = 12, ncol = p)
ac_icmemlfc2 <- matrix(0, nrow = 12, ncol = p)


kp_icmemlfc5 <- matrix(0, nrow = 12, ncol = p)
f1_icmemlfc5 <- matrix(0, nrow = 12, ncol = p)
ac_icmemlfc5 <- matrix(0, nrow = 12, ncol = p)

kp_sc <- matrix(0, nrow = 12, ncol = p)
f1_sc <- matrix(0, nrow = 12, ncol = p)
ac_sc <- matrix(0, nrow = 12, ncol = p)
kp_scina <- matrix(0, nrow = 12, ncol = p)
f1_scina <- matrix(0, nrow = 12, ncol = p)
ac_scina <- matrix(0, nrow = 12, ncol = p)
kp_cellassign <- matrix(0, nrow = 12, ncol = p)
f1_cellassign <- matrix(0, nrow = 12, ncol = p)
ac_cellassign <- matrix(0, nrow = 12, ncol = p)
kp_garnett <- matrix(0, nrow = 12, ncol = p)
f1_garnett <- matrix(0, nrow = 12, ncol = p)
ac_garnett <- matrix(0, nrow = 12, ncol = p)


N_marker <- c(5, 10, 15)
for (iter in (1:12)) {
    ID <- name_ID12[iter] 
    #dlpfc <- readRDS(paste0(dir.exp, ID, "CellAssign.rds") )
    dlpfcSCINA <- readRDS(paste0(dir.exp, ID, "SCINA_Garnett.rds") )
    

    y <- as.character(dlpfcSCINA$layer_guess_reordered)
    y[is.na(y)] <- 'Unknown'
    i =  1
    for (top in N_marker) {
        #for (id in name_ID12[samp==samp[iter]]) {
        for (id in name_ID12[samp != 2])  {
            layer_sp <- colData(dlpfcSCINA)[paste0("icmemlfc", id, "top", top)][[1]]
            res_sp <- evaluate(y, layer_sp, Indices = NULL)
            kp_sp[iter, i] <- cohen.kappa(x=cbind(y, layer_sp))$kappa
            f1_sp[iter, i] <- mean(res_sp$F1)
            ac_sp[iter, i] <- mean(y == layer_sp)

            layer_icmemlfc <- colData(dlpfcSCINA)[paste0("icmemlfc0", id, "top", top)][[1]]
            res_icmemlfc <- evaluate(y, layer_icmemlfc, Indices = NULL)
            kp_icmemlfc0[iter, i] <- cohen.kappa(x=cbind(y, layer_icmemlfc))$kappa
            f1_icmemlfc0[iter, i] <- mean(res_icmemlfc$F1)
            ac_icmemlfc0[iter, i] <- mean(y == layer_icmemlfc)

            layer_icmemlfc <- colData(dlpfcSCINA)[paste0("icmemlfc02", id, "top", top)][[1]]
            res_icmemlfc <- evaluate(y, layer_icmemlfc, Indices = NULL)
            kp_icmemlfc2[iter, i] <- cohen.kappa(x=cbind(y, layer_icmemlfc))$kappa
            f1_icmemlfc2[iter, i] <- mean(res_icmemlfc$F1)
            ac_icmemlfc2[iter, i] <- mean(y == layer_icmemlfc)

            layer_icmemlfc <- colData(dlpfcSCINA)[paste0("icmemlfc05", id, "top", top)][[1]]
            res_icmemlfc <- evaluate(y, layer_icmemlfc, Indices = NULL)
            kp_icmemlfc5[iter, i] <- cohen.kappa(x=cbind(y, layer_icmemlfc))$kappa
            f1_icmemlfc5[iter, i] <- mean(res_icmemlfc$F1)
            ac_icmemlfc5[iter, i] <- mean(y == layer_icmemlfc)

            layer_sc <- colData(dlpfcSCINA)[paste0("layer_scSorter", id,"top", top)][[1]] 
            res_sc <- evaluate(y, layer_sc, Indices = NULL)
            kp_sc[iter, i] <- cohen.kappa(x=cbind(y, layer_sc))$kappa
            f1_sc[iter, i] <- mean(res_sc$F1)
            ac_sc[iter, i] <- mean(y == layer_sc)

            layer_scina <- colData(dlpfcSCINA)[paste0("layer_scina", id,"top", top)][[1]] 
            layer_scina[layer_scina == "unknown"] <- 'Unknown'
            res_scina <- evaluate(y, layer_scina, Indices = NULL)
            kp_scina[iter, i] <- cohen.kappa(x=cbind(y, layer_scina))$kappa
            f1_scina[iter, i] <- mean(res_scina$F1)
            ac_scina[iter, i] <- mean(y == layer_scina)

            layer_cellassign <- colData(dlpfcSCINA)[paste0("layer_cellassign", id,"top", top)][[1]] 
            res_cellassign <- evaluate(y, layer_cellassign, Indices = NULL)
            kp_cellassign[iter, i] <- cohen.kappa(x=cbind(y, layer_cellassign))$kappa
            f1_cellassign[iter, i] <- mean(res_cellassign$F1)
            ac_cellassign[iter, i] <- mean(y == layer_cellassign)

            layer_garnett <- colData(dlpfcSCINA)[paste0("garnett", id,"top", top)][[1]] 
            if (iter %in% (5:8))     {
                dlpfcC <- readRDS(paste0(dir.exp, ID, "CorrectSpecify_Garnett.rds")) 
                layer_garnett <- colData(dlpfcC)[paste0("garnett", id,"top", top)][[1]] 
            }
            res_garnett <- evaluate(y, layer_garnett, Indices = NULL)
            kp_garnett[iter, i] <- cohen.kappa(x=cbind(y, layer_garnett))$kappa
            f1_garnett[iter, i] <- mean(res_garnett$F1)
            ac_garnett[iter, i] <- mean(y == layer_garnett)

            i = i + 1
        }
    }    
}

colnames(kp_sp) <- paste0(name_ID12[samp != 2], "top", rep(N_marker, each=8))
colnames(f1_sp) <- paste0(name_ID12[samp != 2], "top", rep(N_marker, each=8))
colnames(ac_sp) <- paste0(name_ID12[samp != 2], "top", rep(N_marker, each=8))


res <- list(sp = list(kp_sp = kp_sp, f1_sp = f1_sp, ac_sp=ac_sp), 
            icmemlfc0 = list(kp_icmemlfc0 = kp_icmemlfc0, f1_icmemlfc0 = f1_icmemlfc0, ac_icmemlfc0=ac_icmemlfc0), 
            icmemlfc2 = list(kp_icmemlfc2 = kp_icmemlfc2, f1_icmemlfc2 = f1_icmemlfc2, ac_icmemlfc2=ac_icmemlfc2), 
            icmemlfc5 = list(kp_icmemlfc5 = kp_icmemlfc5, f1_icmemlfc5 = f1_icmemlfc5, ac_icmemlfc5=ac_icmemlfc5),
            sc = list(kp_sc = kp_sc, f1_sc = f1_sc, ac_sc=ac_sc), 
            scina = list(kp_scina = kp_scina, f1_scina = f1_scina, ac_scina=ac_scina), 
            cellassign = list(kp_cellassign = kp_cellassign, f1_cellassign = f1_cellassign, ac_cellassign=ac_cellassign),
            garnett = list(kp_garnett = kp_garnett, f1_garnett = f1_garnett, ac_garnett=ac_garnett)
            ) 
saveRDS(res, "brain12_semisupervised.rds")


  
 # -------------------------------------------------------------------------
# produce input for PAGA
setwd("./Real_data_results/dataFiles/DLPFC/Brain12cross/")
library(SingleCellExperiment) 
 
name_ID12 <- as.character(c(151507, 151508, 151509, 151510, 
                            151669, 151670, 151671, 151672, 
                            151673, 151674, 151675, 151676))
  

dir.exp <- "./Real_data_results/dataFiles/DLPFC/Brain12cross/"
dir.spark <- "./Real_data_analysis/DLPFC/Datasets/brain12_spark/"
 
N_marker <- c(5, 10, 15)

for (iter in (1:12)) {
    ID <- name_ID12[iter] 
    dlpfcSCINA <- readRDS(paste0(dir.exp, ID, "SCINA_Garnett.rds") )
    #dlpfc <- readRDS(paste0(dir.exp, ID, ".rds") )
    load(paste0(dir.spark,"brain_", ID,"_spark.Rdata") )
    set.seed(101)
    adjPval <- PvalDF[,2]
    names(adjPval) <- row.names(PvalDF)
    sort_adjPval <- sort(adjPval) 
    num_cut <- 2000
    sp_sig_genes <- names(sort_adjPval)[1:num_cut]
  
    logCount <- assay(dlpfcSCINA, "logcounts")
    sp_logCount <- logCount[sp_sig_genes, ]
  
    X <- as.matrix(t(sp_logCount)) # obtain data

    write.table(X, paste0("dlpfc_", ID, "_count.txt"))

    ## annoated layer 
    top <- 5
    id <- name_ID12[1] 
    layer_sp <- colData(dlpfcSCINA)[paste0("icmemlfc", id, "top", top)][[1]]
    write.table(layer_sp, paste0("dlpfc_", ID, "_celltype.txt"), row.names =FALSE, col.names =FALSE)   

    ## low dimensional representation
    fit_s <- readRDS(paste0("./Real_data_results/dataFiles/DLPFC/unsupervised/", ID, "_icmem.rds"))
    hZ <- fit_s$Ez_u
    write.table(hZ, paste0("dlpfc_", ID, "_Ez.txt"), row.names =FALSE, col.names =FALSE)   
}

