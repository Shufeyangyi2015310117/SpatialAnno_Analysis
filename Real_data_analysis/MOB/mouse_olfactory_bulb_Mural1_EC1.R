set.seed(99)
setwd("/home/yangyi/Xingjie/MOB2")
load("Rep12_MOB_count_matrix-1.RData")
library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

pos = cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(MOB_raw),split="x"),"[",1)),y=as.numeric(sapply(strsplit(colnames(MOB_raw),split="x"),"[",2)))
X <- t(as.matrix(MOB_raw))

MOB <- CreateSeuratObject(counts = MOB_raw, project = "MOB")

MOB = SCTransform(MOB)

X = t(MOB@assays$SCT@scale.data)


DEG = readRDS("/home/yangyi/Xingjie/MOB2/markerList_MOB_zhou.rds")

cutoff <- c(1.5, 1)
#for(j in 1:2) { try different makers
j = 1
K_Lay <- length(unique(DEG$cluster))
markers <- vector("list", K_Lay)

## intersection of colnames(X) and DEG$gene
int_gene = intersect(colnames(X), DEG$gene)
idx = match(colnames(X), DEG$gene)


for (i in 1:K_Lay){
  if (i==1){
    DEG_i = dplyr::filter(DEG, cluster == as.character(unique(DEG$cluster)[i]))
    idx = match(int_gene, DEG_i$gene)
    DEG2_i = DEG_i[sort(idx[is.na(idx)!=T]),]
    DEG2 = DEG2_i
  }else{
    DEG_i = dplyr::filter(DEG, cluster == as.character(unique(DEG$cluster)[i]))
    idx = match(int_gene, DEG_i$gene)
    DEG2_i = DEG_i[sort(idx[is.na(idx)!=T]),]
    DEG2 = rbind(DEG2, DEG2_i)
  }
}
dim(DEG2)
library(dplyr)
## DEG2 is the DEG that are included in the 2000 highly variable genes


for (k in 1:K_Lay){
  
  signature_mat <- dplyr::filter(DEG2, cluster == as.character(unique(DEG2$cluster)[k]) & 
                                   avg_log2FC > 0 & p_val_adj < 0.05) %>%
    top_n(n = 4, wt = avg_log2FC)
  signature = signature_mat$gene
  print(length(signature))
  
  
  if (length(signature)>0){
    idx = match(signature, colnames(X))
    names(markers)[k] = as.character(unique(DEG2$cluster)[k])
    markers[[k]] = signature[!is.na(idx)]
  }else{
    print("else")
    
    signature_mat <- dplyr::filter(DEG2, cluster == as.character(unique(DEG2$cluster)[k]) & 
                                     avg_log2FC > 0) %>%
      top_n(n = 4, wt = avg_log2FC)
    signature = signature_mat$gene
    print(length(signature))
    
    idx = match(signature, colnames(X))
    names(markers)[k] = as.character(unique(DEG2$cluster)[k])
    markers[[k]] = signature[!is.na(idx)]
  }
}

markers[[1]][4] = "Penk"
markers[[2]][4] = "Th"
idx = match(c("Penk", "Th"), colnames(X))
markers2 = vector("list", K_Lay+2)
names(markers2)[1:5] = names(markers)
names(markers2)[6] = "Mural1"
names(markers2)[7] = "EC1"
markers2[1:5] = markers
markers2[[6]] = c("Cald1","Slco1a4", "Ly6c1","Igfbp7")
markers2[[7]] = c("Ly6c1","Slco1a4", "Ly6a","Cldn5")
markers = markers2

save(markers, file = "Rep12_MOB_4markers_replace_with_Th_Penk_Mural1_EC1.RData")


#------------------------------------------------------------------------------------------------------------
#setwd("~/Cloud/Research/SSL/code/")
setwd("/home/yangyi/Xingjie/run")
source("function.R")
source("function_f1score.R")



rho <- marker_list_to_mat(markers, TRUE)
#rho <- rho[, which(colSums(rho) !=0 )]

err <- matrix(0, 1, 10)
acc <- matrix(0, 1, 10)
f1s <- matrix(0, 1, 10)
time <- matrix(0, 1, 10)
kappa <- matrix(0, 1, 10)


Adj_sp <- getneighborhood_fast(as.matrix(pos), cutoff = 1.2)
fit_s2 <- SpatialAnno(X = X, Adj_sp = Adj_sp, marker = marker, initial = "SCINA")  


library(scSorter)
dat <- t(X)
tic <- proc.time()
rts <- scSorter(dat, anno)
toc <- proc.time()

setwd("/home/yangyi/Xingjie/MOB2")
save(err, acc, f1s, kappa, time, results, fit_s2, rho, pos, rts, file = "Rep12_MOB_4markers_SCINA_initial_icmem_lfc_replace_with_Penk_Th_Mural1_EC1.RData")






