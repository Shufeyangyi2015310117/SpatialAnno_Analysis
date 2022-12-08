
# run SpatialAnno, SCINA, scSorter, CellAssign

library(SingleCellExperiment) 
library(SpatialAnno) 

 
dlpfcCross <- function(iter) {
  dir.spark <- "/Real_data_analysis/DLPFC/Datasets/brain12_spark/"
  dir.exp <- "./Real_data_results/dataFiles/DLPFC/unsupervised/"
  ticAll <- Sys.time()
  name_ID12 <- as.character(c(151507, 151508, 151509, 151510, 
                              151669, 151670, 151671, 151672, 
                              151673, 151674, 151675, 151676))
  samp <- rep(1:3, each=4)
  num_cut <- 2000
  #trueK_set <- c(rep(7,4), rep(5,4), rep(7,4))
  #------------------------------------------
  # load and read data
  ID <- name_ID12[iter] 

  ## get the first 2000 genes in SVGs order
  dlpfc <- readRDS(paste0(dir.exp, ID, ".rds") )
  load(paste0(dir.spark,"brain_", ID,"_spark.Rdata") )
  set.seed(101)
  adjPval <- PvalDF[,2]
  names(adjPval) <- row.names(PvalDF)
  sort_adjPval <- sort(adjPval) 
  sp_sig_genes <- names(sort_adjPval)[1:num_cut]
  
  logCount <- assay(dlpfc, "logcounts")
  sp_logCount <- logCount[sp_sig_genes, ]
  
  X <- as.matrix(t(sp_logCount)) # obtain data
  pos <- cbind(dlpfc$row, dlpfc$col) 
  p <- ncol(X); n <- nrow(X)
  #  make BayesSpace metadata used in BayesSpace-------------------------------------------------
  counts <- t(X)
  rownames(counts) <- paste0("gene_", seq_len(p))
  colnames(counts) <- paste0("spot_", seq_len(n))
  
  ## Make array coordinates - filled rectangle
  cdata <- list()
  cdata$row <- pos[,1]
  cdata$col <- pos[,2]
  cdata <- as.data.frame(do.call(cbind, cdata))
  
  cdata$imagerow <- cdata$row
  cdata$imagecol <- cdata$col 
  ## Make SCE
  sce <- SingleCellExperiment(assays=list(logcounts=counts), colData=cdata)
  princ <- princomp(X)
  reducedDim(sce, "PCA") <- princ$scores[,1:50]
  # hq <- selectFacNumber(X)$q
  
  y <- as.character(dlpfc$layer_guess_reordered)
  y[is.na(y)] <- 'Unknown'

  K <- length(unique(y))


  # calculate the Adjoint matrix
  library(purrr)
  library(Matrix)
  Adj_sp <- find_neighbors2(sce, platform="Visium")
  #------------------------------------------
  #K <- trueK_set[iter]
  hq <- 15 # default as 15 factors
  set.seed(20132014)
  #----------------------------------------------------------prepare a marker list
  #----------------------------------------------------------
  marker_dir <- "./Real_data_analysis/DLPFC/Datasets/markerfiles/"
  for (top in c(5, 10, 15)) {
    for(id in name_ID12[samp != 2])  {
      marker_file <- paste0(id, "_DEGmarkerTop_", top, ".rds")
      markers <- readRDS(paste0(marker_dir, marker_file))
      if (ID %in% name_ID12[5:8]) {
        markers$Layer1 = NULL
        markers$Layer2 = NULL
      }  
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
      # scSorter
      library(scSorter)
      tic <- proc.time()
      rts <- scSorter(t(X), anno, alpha = 0)
      toc <- proc.time()
      colData(dlpfc)[paste0("layer_scSorter", id, "top", top)] <- rts$Pred_Type
      colData(dlpfc)@metadata[paste0("time_scSorter", id, "top", top)] <- (toc - tic)[3]

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
      

      # run CellAssign
      Sys.setenv(RETICULATE_PYTHON = "./miniconda2/envs/py37/bin/python")
      library(cellassign)
      Count <- assay(dlpfc, "counts")
      mat <- as.matrix(Count[sp_sig_genes, ])
      share_gene <- intersect(rownames(rho), rownames(mat))
      mat2 <- mat[share_gene, ]
      spot_idx <- which(colSums(mat2) != 0)
      mat3 <- mat2[, spot_idx]
      tic <- proc.time()
      fit <- cellassign(exprs_obj = t(mat3),
                          marker_gene_info = rho,
                          s = sizeFactors(dlpfc)[spot_idx], min_delta = 0.2,
                          learning_rate = 1e-2, shrinkage = TRUE, verbose = FALSE)
      toc <- proc.time()
      y_hat <- rep("Unknown", n)
      y_hat[spot_idx] <- fit$cell_type
      colData(dlpfc)[paste0("layer_cellassign", id, "top", top)] <- y_hat
      colData(dlpfc)@metadata[paste0("time_cellassign", id, "top", top)] <- (toc - tic)[3]

      # run SCINA
      library(SCINA)
      tic <- proc.time()
      results = SCINA(dat$dat, markers, max_iter = 100, convergence_n = 10, rm_overlap = 0)
      toc <- proc.time()
      colData(dlpfc)[paste0("layer_scina", id, "top", top)] <- results$cell_labels    
      colData(dlpfc)@metadata[paste0("time_scina", id, "top", top)] <- (toc - tic)[3]

      
      # SpatialAno
      tic <- proc.time()
      fit_s <- SpatialAnno(X = X, Adj_sp = Adj_sp, marker = markers, initial = "SCINA")
      toc <- proc.time()
      colData(dlpfc)@metadata[paste0("time_icmemlfc", id, "top", top)] <- (toc - tic)[3]
      colData(dlpfc)[paste0("icmemlfc", id, "top", top)] <- colnames(rho)[fit_s$type]

    }
  }
  tocAll <- Sys.time()
  print(tocAll - ticAll)

  outpath <- paste0("./Real_data_results/dataFiles/DLPFC/Brain12cross/")
  setwd(outpath)
  saveRDS(dlpfc, paste0(ID, "SCINA.rds"))
}
