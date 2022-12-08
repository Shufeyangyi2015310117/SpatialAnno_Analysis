 
setwd("./Real_data_analysis/Hippo/Datasets/")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
#slide.seq <- LoadData("ssHippo")
seulist <- readRDS("seulist_Hip2_minFeature15_preprocess.rds")

# 1 or 2
embro <- 1
slide.seq <- seulist[[embro]]
 
pos <- cbind(slide.seq$row, slide.seq$col)
library(Matrix)
Adj_sp <- getneighborhood_fast(as.matrix(pos), cutoff = 30)
summary(rowSums(Adj_sp))

slide.seq <- SCTransform(slide.seq, assay = "RNA", return.only.var.genes = FALSE, verbose = FALSE)
X <- t(slide.seq@assays$SCT@scale.data)

iter <- 1
  
dataSource <- 1
mk_n <-  c(5, 10)
markers <-  readRDS("./hippo_markers.rds")
 
K <- length(markers)
rho <- marker_list_to_mat(markers, TRUE)
    
df_all = as.data.frame(matrix(0,0,2))
colnames(df_all) = c("Marker", "Type")
for (i in 1:K){
    if (sum(rho[,i]) != 0)  {
        print(rownames(rho)[rho[,i]==1])
        df = as.data.frame(rownames(rho)[rho[,i]==1])
        colnames(df) = "Marker"
        df$Type = colnames(rho)[i]
        df_all = rbind(df_all, df)
    }
}
anno = df_all[, c(2,1)]
#---------------------------
 
  library(scSorter)
  tic <- proc.time()
  rts <- scSorter(t(X), anno, alpha = 0)
  toc <- proc.time()
  print(toc - tic)
  slide.seq@meta.data["layer_scSorter"] <- rts$Pred_Type
  slide.seq@meta.data["time_scSorter"] <- (toc - tic)[3]

    
#---------------------------
# prepared data, markers are in the same order of columns as scSorters.
      anno_processed = scSorter:::design_matrix_builder(anno, weight=2)
      dat <- scSorter:::data_preprocess(t(X), anno_processed)
      dat$designmat$Unknown = 0
      rho <- as.matrix(dat$designmat)
      m <- nrow(rho)
      X_m <- t(dat$dat[1:m, ])
      X_u <- t(dat$dat[-(1:m), ])
      n <- nrow(X_m)
      K <- ncol(rho)
      
if(FALSE){
      Sys.setenv(RETICULATE_PYTHON = "./miniconda2/envs/py37/bin/python")
      library(cellassign)
      Count <- slide.seq@assays$RNA@counts  
      mat <- as.matrix(Count)
      s <- colSums(mat)/mean(colSums(mat))
      rownames(mat) <- toupper(rownames(mat))
      share_gene <- intersect(rownames(rho), rownames(mat))
      mat2 <- mat[share_gene, ]
      spot_idx <- which(colSums(mat2) != 0)
      mat3 <- mat2[, spot_idx]
      tic <- proc.time()
      fit <- cellassign(exprs_obj = t(mat3),
                          marker_gene_info = rho,
                          s = s[spot_idx], min_delta = 0.2,  
                          learning_rate = 1e-2, shrinkage = TRUE, verbose = FALSE)
      toc <- proc.time()
      y_hat <- rep("Unknown", n)
      y_hat[spot_idx] <- fit$cell_type
      slide.seq@meta.data["layer_cellassign"] <- y_hat
      slide.seq@meta.data["time_cellassign"] <- (toc - tic)[3]
}

      library(SCINA)
      tic <- proc.time()
      results = SCINA(dat$dat, markers, max_iter = 100, convergence_n = 10, rm_overlap = 0)
      toc <- proc.time()
      slide.seq@meta.data["layer_scina"] <- results$cell_labels
      slide.seq@meta.data["time_scina"] <- (toc - tic)[3]

      output <- "~/Real_data_results/dataFiles/Hippo/"
 
##############################################################################
      toc <- proc.time()
      fit_s <- SpatialAnno(X = X, Adj_sp = Adj_sp, marker = markers, initial = "SCINA")
      toc <- proc.time()
      slide.seq@meta.data["icmemlfc"] <- colnames(rho)[fit_s$type]
      slide.seq@meta.data["time_icmemlfc"] <- (toc - tic)[3]
      print((toc - tic)[3])
     
      saveRDS(fit_s, paste0(output, "hippoAll_Manual_icmem_", embro, ".rds")) 
      saveRDS(slide.seq, paste0(output, "hippoAll_anno_", "satija", mk_n[iter], "_", embro, ".rds")) 
 



    