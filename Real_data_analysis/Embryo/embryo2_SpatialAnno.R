setwd("/home/yangyi/Xingjie")
library(Seurat)
library(SpatialAnno)
Emall <- readRDS('counts.Rds')
meta = readRDS("metadata.Rds")
head(meta)
idx1 = which(meta$embryo == "embryo2")
idx = which(meta$embryo == "embryo1")

### Take out sample 3
y <- meta$celltype_mapped_refined[idx]
X <- t(as.matrix(Emall[,idx1])) 
pos <- cbind(meta$x_global[idx1], meta$y_global[idx1])


library(SingleCellExperiment)
library(BayesSpace)
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


cutoff <- c(1.5, 1)
#for(j in 1:2) { try different makers
j = 1
K_Lay <- 24
markers <- vector("list", K_Lay)
freq_Lay <- table(y)

DEG = read.csv("embryo3_DEGs_top30_annotatedLabels.csv")


library(dplyr)

for (k in 1:24){
  
  signature_mat <- dplyr::filter(DEG, cluster == as.character(unique(y)[k]) & 
                                   avg_log2FC > 1 & p_val_adj < 1e-5) %>%
    top_n(n = 8, wt = avg_log2FC)
  signature = signature_mat$gene
  print(length(signature))
  
  
  if (length(signature)>0){
    idx = match(signature, colnames(X))
    names(markers)[k] = as.character(unique(y)[k])
    markers[[k]] = signature[!is.na(idx)]
  }else{
    
    signature_mat <- dplyr::filter(DEG, cluster == as.character(unique(y)[k]) & 
                                     avg_log2FC > 1) %>%
      top_n(n = 4, wt = avg_log2FC)
    signature = signature_mat$gene
    print(length(signature))
    
    idx = match(signature, colnames(X))
    names(markers)[k] = as.character(unique(y)[k])
    markers[[k]] = signature[!is.na(idx)]
  }
}

y <- meta$celltype_mapped_refined[idx1]

#------------------------------------------------------------------------------------------------------------
#setwd("~/Cloud/Research/SSL/code/")
setwd("/home/yangyi/Xingjie/run")
source("function.R")
source("function_f1score.R")

markers2 = markers[-c(1,9,18)]
markers2[[4]]
markers2[[18]]
markers3 = markers2
load("/home/yangyi/Xingjie/run/embryo1_markers2.RData")
markers3[[18]] = markers2[[18]]
markers3[[4]] = markers2[[4]]
markers2 = markers3
markers2[[18]]
markers2[[4]]

rho <- marker_list_to_mat(markers2, TRUE)
#rho <- rho[, which(colSums(rho) !=0 )]
#save(rho, file = "embryo2_rho_8marker_presomitic_8m_Gut.RData")

err <- matrix(0, 1, 10)
acc <- matrix(0, 1, 10)
f1s <- matrix(0, 1, 10)
time <- matrix(0, 1, 10)
kappa <- matrix(0, 1, 10)


marker = markers2
Adj_sp <- getneighborhood_fast(as.matrix(pos), cutoff = 1.2)
fit_s2 <- SpatialAnno(X = X, Adj_sp = Adj_sp, marker = marker, initial = "scSorter")  

save(err, acc, f1s,kappa, time, fit_s2, rho,  pos, y, file = "embryo2_allcelltype_unknown_8markers_scSorter_initial_icmem_lfc=_8marker_presomitic_8m_Gut_adjusted.RData")






