setwd("~/Xingjie/simulation")
set.seed(1234)
library(BayesSpace)
library(DR.SC)
library(Seurat)
library(SCINA)
library(psych)
library(SingleCellExperiment)
library(Rcpp)
sourceCpp("~/Xingjie/run/em_jigsaw.cpp")
sourceCpp("~/Xingjie/run/wpca.cpp")
#library(garnett)
library(scuttle)
#library(cellassign)
library(scSorter)
source("~/Xingjie/run/function_f1score.R")


err <- matrix(0, 50, 10)
acc <- matrix(0, 50, 10)
f1s <- matrix(0, 50, 10)
time <- matrix(0, 50, 10)
kappa <- matrix(0, 50, 10)
ARI <- matrix(0, 50, 10)

## read position and annotation label from real data dlpfc
dlpfc <- getRDS("2020_maynard_prefrontal-cortex", "151673")
pos = colData(dlpfc)[,c("row", "col")]

for (iter in 1:50){
  print(iter)
  set.seed(NULL)
  
  n <- dim(dlpfc)[2] # number of spots
  p <- 100 ## number of non-markers
  y <- as.numeric(colData(dlpfc)[,c("layer_guess_reordered")])
  y[is.na(y)] = 8
  y2 <- paste0("ct", y)
  y2[y2=="ct8"] = "Unknown"
  K <- length(table(y)) -1 ## number of clusters
  num_mk_per_ct = 5 ## number of markers per cell type
  m <- num_mk_per_ct  * K ## number of markers
  

  generate_count <- function(J = 100, de_facLoc = 0, de_facScale = 1){
    library(BayesSpace)
    library(SingleCellExperiment)
    library(splatter)
    setwd("~/Xingjie/simulation")
    ## read position and annotation label from real data dlpfc
    dlpfc <- getRDS("2020_maynard_prefrontal-cortex", "151673")
    pos = colData(dlpfc)[,c("row", "col")]
    y <- as.numeric(colData(dlpfc)[,c("layer_guess_reordered")])
    y[is.na(y)] = 8
    
    dec <- scran::modelGeneVar(dlpfc)
    top <- scran::getTopHVGs(dec, n = 2000)
    cnts = as.matrix(counts(dlpfc)[top,])
    init_params <- splatEstimate(cnts)
    
    num_mk_per_ct = 5
    batch_facLoc = 0
    C = 8
    I = NULL
    N = 7000
    L = 1
    
    de_prop = rep(0.5,8)
    sim_seed = as.integer(Sys.time())
    debug = FALSE
    
    # 1.simulate count data
    noBatch <- ifelse(batch_facLoc == 0, TRUE, FALSE)
    group_prob <- as.vector(table(y)/length(y))
    
    
    params <- setParams(
      init_params, 
      batchCells = rep(N, L), # 3N here represents a large number such that
      # we have sufficient cells of each type to be 
      # allocated to the spatial transcriptomics data
      batch.rmEffect = noBatch,
      batch.facLoc = batch_facLoc,
      nGenes = J,
      group.prob = group_prob,
      out.prob = 0,
      de.prob = de_prop,
      de.facLoc = de_facLoc,
      de.facScale = de_facScale,
      seed = sim_seed)
    
    sim_groups <- splatSimulate(
      params = params, 
      method = "groups", 
      verbose = FALSE)
    
    library(Seurat)
    library(scuttle)
    sim_groups <- logNormCounts(sim_groups)
    seu = as.Seurat(sim_groups)
    seu = SCTransform(seu, assay = "originalexp")
    Idents(seu) = seu@meta.data$Group
    all.markers = FindAllMarkers(seu, assay = "SCT", logfc.threshold = 0.1)
    
    library(dplyr)
    all.markers %>%
      group_by(cluster) %>%
      top_n(n = num_mk_per_ct, wt = avg_log2FC) -> top5
    
    out = list(top5 = top5, sim_groups = sim_groups, all.markers = all.markers)
    return(out)
    
  }
  
  ## generate the count
  res1 = generate_count(J = 100, de_facLoc = c(1,2,3,4,5,6,7,8)*0.1, de_facScale = c(0.1,0.1,1,0.1,1,1,0.10,0.1))
  res2 = generate_count(J = 2000, de_facLoc = c(1,2,3,4,5,6,7,8)*0.1, de_facScale = c(0.1,0.1,1,0.1,1,1,0.10,0.1)*0.4)
  
  
  res1$all.markers %>%
    group_by(cluster) %>%
    top_n(n = 15, wt = avg_log2FC) -> top15
  adjusted_top5_gene = top15$gene[c(11:15,26:30,31:35,56:60,61:65,86:90,101:105)]
  ##################################1:5,16:20,31:35,46:50,61:65,76:80,91:95
  ### c(11:15,26:30,31:35,56:60,61:65,86:90,101:105)
  
  
  ## reorder the sample to  be the same as y
  Groupy = paste0("Group", y)
  num_each_celltype = table(y)
  idx1 = rep(0, length(y))
  for (i in 1:8){
    idx1[y==i] = which(colData(res1$sim_groups)[,3] == paste0("Group",i))[1:num_each_celltype[i]] 
  }
  
  idx2 = rep(0, length(y))
  for (i in 1:8){
    idx2[y==i] = which(colData(res2$sim_groups)[,3] == paste0("Group",i))[1:num_each_celltype[i]] 
  }
  
  
  ## ordered samples
  X1 = counts(res1$sim_groups)[unique(adjusted_top5_gene),idx1]
  rownames(X1) = gsub("Gene","mk", rownames(X1))
  X2 = counts(res2$sim_groups)[,idx2]
  colsum = apply(X2,1,sum)
  if (min(colsum) == 0){
    idx_not_equal_zero = which(colsum!=0)
    X2 = X2[idx_not_equal_zero,]
  }
  X2 = X2[1:1000,]
  
  print(all(table(colData(res1$sim_groups)[idx1,3]) == table(colData(res2$sim_groups)[idx2,3])))
  
  ## generate rho
  rho <- matrix(0, dim(X1)[1], K+1)
  rownames(rho) <- gsub("Gene","mk", rownames(X1))
  colnames(rho)[1:(K+1)] <- paste0("ct", 1:(K+1))
  colnames(rho)[K+1] <- "Unknown"
  res1_top5_gene = gsub("Gene","mk", adjusted_top5_gene)
  for (k in 1:K) {
    rho[res1_top5_gene[((k-1)*num_mk_per_ct + 1):(k*num_mk_per_ct)], k] <- 1
  }
  
  ## define markers
  marker = list()
  for (k in 1:K){
    marker[[k]] = rownames(rho)[rho[,k]==1]
    names(marker)[k] = colnames(rho)[k]
  }
  
  combined_counts = rbind(X1, X2)
  fit_pca = prcomp(t(combined_counts))
  fit_y = fit_pca$x[,1:15]
  # -------------------------------------------------
  # make BayesSpace metadata used in BayesSpace
  counts <- combined_counts
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
  
  sce <- logNormCounts(sce)
  X = t(logcounts(sce))
  print(X[1:10,1:10])
  
    
  
  df_all = as.data.frame(matrix(0,0,2))
  colnames(df_all) = c("Marker", "Type")
  for (i in 1:(dim(rho)[2]-1)){
    print(rownames(rho)[rho[,i]==1])
    df = as.data.frame(rownames(rho)[rho[,i]==1])
    colnames(df) = "Marker"
    df$Type = colnames(rho)[i]
    df_all = rbind(df_all, df)
  }
  anno = df_all[,c(2,1)]
  
  
  # --------------------------------------------------------------------
  # Unknown is a swift controlling whether novel celltype is allowed.
  Unknown = TRUE
  anno_processed = scSorter:::design_matrix_builder(anno, weight=2)
  dat <- scSorter:::data_preprocess(t(X), anno_processed)
  if (Unknown == TRUE) dat$designmat$Unknown = 0
  rho <- as.matrix(dat$designmat)
  m <- nrow(rho)
  X_m <- as.matrix(t(dat$dat[1:m, ]))
  X_u <- as.matrix(t(dat$dat[-(1:m), ]))
  n <- nrow(X_m)
  K_unknown <- ncol(rho)
  
  
  dat2 = t(X)
  tic <- proc.time()
  results = SCINA(dat2, marker, max_iter = 100, convergence_n = 10, rm_overlap = 0, allow_unknown=TRUE)
  toc <- proc.time()
  time[iter, 1] = (toc - tic)[3]

  index1 = which(y2!="Unknown")
  index2 = which(results$cell_labels!="unknown")
  index3 = intersect(index1,index2)
  err[iter,1] = 1 - mean(results$cell_labels[index3] == as.character(y2)[index3])
  f1score = evaluate(as.character(y2)[index3], results$cell_labels[index3])
  acc[iter,1] = f1score$Acc
  f1s[iter,1] = mean(f1score$F1)
  kappa[iter,1] = cohen.kappa(x=cbind(as.character(y2)[index3], results$cell_labels[index3]))$kappa



  #### SpatialAnno ####################################################
  bet_min <- numeric(length(results$theta))
  bet_max <- numeric(length(results$theta))
  alpha_int <- numeric(m)
  bet_int <- matrix(0, nrow = m, ncol = K_unknown)
  for(i in 1:length(results$theta)){
    bet_max[i] <- max(results$theta[[i]]$mean[,1] - results$theta[[i]]$mean[,2])
    bet_min[i] <- min(results$theta[[i]]$mean[,1] - results$theta[[i]]$mean[,2])

    bet_int[rho[,i]!=0, i] <- results$theta[[i]]$mean[, 1] - results$theta[[i]]$mean[, 2]
    alpha_int[rho[,i]!=0] <- results$theta[[i]]$mean[, 2]
  }
  lfc <- median(bet_min)



  ################################################
  ## SCINA initial value
  y_hat <- results$cell_labels;
  y_hat[y_hat == "unknown"] = "Unknown"
  y_hat <- results$cell_labels;
  y_hat[y_hat == "unknown"] = "Unknown"
  if (length(unique(y_hat))<K_unknown){
    idx = match(names(marker), unique(y_hat))
    idx2 = which(is.na(idx)==TRUE)
    if (length(idx2) > 0){
      y_hat[1:length(names(marker)[idx2])]=names(marker)[idx2]
    }
    if (is.na(match("Unknown", unique(y_hat)))){
      y_hat[length(names(marker)[idx2])+1]="Unknown"
    }
  }
  y_int <- match(y_hat, colnames(rho))
  R_int <- matrix(0, nrow = n, ncol = K_unknown)
  for(i in 1:nrow(X_m)) R_int[i, match(y_hat, colnames(rho))[i]] <- 1
  Pi_u_int <- colMeans(R_int)
  sigma_int <- update_sigma(X_m, rho, R_int, alpha_int, bet_int)


  princ <- wpca(X_u, q=15, weighted=TRUE)
  Lam_u_int <- princ$Lam_vec
  W_u_int <- princ$loadings
  hZ <- princ$PCs
  n_c <- colSums(R_int)
  Mu_u_int <- t(sapply(1:K_unknown, function(k) 1/n_c[k] * colSums(R_int[, k] * hZ)))
  q <- ncol(Mu_u_int)
  Sgm_u_int <- init_Sgm(R_int, hZ, matrix(0, ncol=q, nrow=q), Mu_u_int, FALSE)
  # -----------------------------------------------------------------



  ## alternative cutoff
  library(SC.MEB)
  library(Matrix)
  # markov random fields
  Adj_sp <- SC.MEB:::find_neighbors2(sce, "Visium")
  summary(rowSums(Adj_sp))


  xi_int <- 1.5
  xi_grid=seq(0.1, 2.5, by=0.2)
  # sourceCpp("em_jigsaw.cpp")
  tic <- proc.time()
  fit_s2 <- icmem(X_m, X_u, Adj_sp, rho, lfc,
                  y_int, Pi_u_int*0, xi_int, xi_grid,
                  alpha_int, bet_int, sigma_int,
                  Mu_u_int, W_u_int, Sgm_u_int, Lam_u_int,
                  300, 10, 1e-6, TRUE,
                  FALSE, FALSE)
  toc <- proc.time()
  time[iter, 2] = (toc - tic)[3]
  index1 = which(y2!="Unknown")
  index2 = which(colnames(rho)[fit_s2$type]!="Unknown" )
  index3 = intersect(index1,index2)
  err[iter, 2] = 1 - mean(colnames(rho)[fit_s2$type][index3] == as.character(y2)[index3])
  f1score = evaluate(as.character(y2)[index3], colnames(rho)[fit_s2$type][index3])
  acc[iter, 2] = f1score$Acc
  f1s[iter, 2] = mean(f1score$F1)
  kappa[iter, 2] = cohen.kappa(x=cbind(as.character(y2)[index3], colnames(rho)[fit_s2$type][index3]))$kappa
  print(kappa[iter, 1:7])
  
  ######## scSorter ###############################################
  dat <- t(X)
  run_scSorter = function(dat, anno){
    rts <- scSorter(dat, anno)
    return(rts)
  }
  tic <- proc.time()
  rts = run_scSorter(dat, anno)
  toc <- proc.time()

  time[iter, 3] = (toc - tic)[3]
  index1 = which(y2!="Unknown")
  index2 = which(rts$Pred_Type!="Unknown")
  index3 = intersect(index1,index2)
  err[iter, 3] = 1 - mean(rts$Pred_Type[index3] == as.character(y2)[index3])
  f1score = evaluate(as.character(y2)[index3], rts$Pred_Type[index3])
  acc[iter, 3] = f1score$Acc
  f1s[iter, 3] = mean(f1score$F1)
  kappa[iter, 3] <- cohen.kappa(x=cbind(as.character(y2)[index3], rts$Pred_Type[index3]))$kappa
  

  save(acc, f1s, kappa, err, time, file = "~/Xingjie/simulation/simulation162_splatter_method1-2.RData")
  
  
  ##### Garnett #################################################
  mat <- counts(sce)
  fdata <- data.frame(gene_short_name = rownames(mat), num_cells_expressed = rowSums(mat != 0))
  pdata <- as.data.frame(matrix(1, dim(mat)[2], 2))
  rownames(pdata) = colnames(mat)
  # create a new CDS object
  pd <- new("AnnotatedDataFrame", data = pdata)
  fd <- new("AnnotatedDataFrame", data = fdata)
  pbmc_cds <- newCellDataSet(as(mat, "dgCMatrix"),
                             phenoData = pd,
                             featureData = fd)
  # generate size factors for normalization later
  pbmc_cds <- estimateSizeFactors(pbmc_cds)


  sink("/home/yangyi/Xingjie/run/simulation1_markers.txt")
  for (i in 1:length(marker)){
    cat(paste0(">",names(marker)[i],"\n"))
    cat(paste0("expressed: ", paste(marker[[i]][1], marker[[i]][2], marker[[i]][3],
                                    marker[[i]][4], marker[[i]][5], sep = ", ")))
    cat("\n")
  }
  sink()
  file.show("/home/yangyi/Xingjie/run/simulation1_markers.txt")

  tic <- proc.time()
  marker_file_path <- "/home/yangyi/Xingjie/run/simulation1_markers.txt"
  pbmc_classifier <- train_cell_classifier(cds = pbmc_cds,
                                           marker_file = marker_file_path,
                                           db="none",
                                           num_unknown = 50)
  pbmc_cds <- classify_cells(pbmc_cds, pbmc_classifier,
                             db = "none",
                             cluster_extend = TRUE)
  toc <- proc.time()

  time[iter, 5] = (toc - tic)[3]
  layer_garnett = pData(pbmc_cds)[,5]
  index1 = which(y2!="Unknown")
  index2 = which(layer_garnett!="Unknown")
  index3 = intersect(index1,index2)
  err[iter, 5] = 1 - mean(layer_garnett[index3] == as.character(y2)[index3])
  f1score = evaluate(as.character(y2)[index3], layer_garnett[index3])
  acc[iter, 5] = f1score$Acc
  f1s[iter, 5] = mean(f1score$F1)
  kappa[iter, 5] <- cohen.kappa(x=cbind(as.character(y2)[index3], layer_garnett[index3]))$kappa

  save(acc, f1s, kappa, err, file = "simulation1_method1-5.RData")
  
  
}

