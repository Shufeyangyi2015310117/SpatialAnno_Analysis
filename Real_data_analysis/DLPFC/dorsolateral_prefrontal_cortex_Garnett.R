
# run garnett
library(SingleCellExperiment) 
library(garnett)
setwd("./SSL/")
library("org.Hs.eg.db") 
 
dlpfcCross <- function(iter) {
  dir.spark <- "./Real_data_analysis/DLPFC/Datasets/brain12_spark/"
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

  dlpfc_Anno <- readRDS(paste0("./Real_data_results/dataFiles/DLPFC/Brain12cross/", ID, "SCINA.rds"))
  ## get the first 2000 genes in SVGs order
  dlpfc <- readRDS(paste0(dir.exp, ID, ".rds") )
  load(paste0(dir.spark,"brain_", ID,"_spark.Rdata") )
  set.seed(101)
  adjPval <- PvalDF[,2]
  names(adjPval) <- row.names(PvalDF)
  sort_adjPval <- sort(adjPval) 
  sp_sig_genes <- names(sort_adjPval)[1:num_cut]
  
   
  #------------------------------------------
  set.seed(20132014)
  #----------------------------------------------------------prepare a marker list
  #----------------------------------------------------------
  marker_dir <- "./Real_data_analysis/DLPFC/Datasets/markerfiles/"
  for (top in c(5, 10, 15)) {
    #for(id in name_ID12[samp==samp[iter]])  {
    for(id in name_ID12[samp != 2])  {
      marker_file <- paste0(id, "_DEGmarkerTop_", top, ".rds")
      markers <- readRDS(paste0(marker_dir, marker_file))
      if (ID %in% name_ID12[5:8]) {
        markers$Layer1 = NULL
        markers$Layer2 = NULL
      }  
       
    # make maker files for garnett

    markers2 <- markers
    if(!file.exists(paste0(marker_dir, id, "_DEGmarkerTop_", top, ".txt")))  {
      sink(paste0(marker_dir, id, "_DEGmarkerTop_", top, ".txt"))
      for (i in 1:length(markers2)){
        cat(paste0(">",names(markers2)[i],"\n"))
        cat("expressed: ")
        cat(paste(markers2[[i]], collapse = ", "))
        cat("\n")
        cat("\n")
      }
      sink()   
    }
# -----------
      Count <- assay(dlpfc, "counts")
      mat <- Count[sp_sig_genes, ]
      #symbol <- mapIds(org.Hs.eg.db, keys = rownames(mat), keytype = "ENSEMBL", column="SYMBOL")
      #rownames(mat) <- unname(symbol)


      fdata <- data.frame(gene_short_name = rownames(mat), num_cells_expressed = rowSums(mat != 0))
      #rownames(fdata) <- rownames(mat)
      pdata <- as.data.frame(matrix(1, ncol(mat), 2))
      rownames(pdata) = colnames(mat)


      # create a new CDS object
      pd <- new("AnnotatedDataFrame", data = pdata)
      fd <- new("AnnotatedDataFrame", data = fdata)
      cds <- newCellDataSet(as(mat, "dgCMatrix"),
                                  phenoData = pd,
                                  featureData = fd)


# -----------
    # generate size factors for normalization later
    cds <- estimateSizeFactors(cds)

    library(org.Hs.eg.db)
    marker_file_path <- paste0(marker_dir, id, "_DEGmarkerTop_", top, ".txt")
    marker_check <- check_markers(cds, marker_file_path,
                              db=org.Hs.eg.db,
                              cds_gene_id_type = "ENSEMBL",
                              marker_file_gene_id_type = "ENSEMBL")

    #plot_markers(marker_check); dev.off()
    set.seed(260)
    classifier <- train_cell_classifier(cds = cds,
                                         marker_file = marker_file_path,
                                         db=org.Hs.eg.db,
                                         cds_gene_id_type = "ENSEMBL",
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "ENSEMBL")



    feature_genes <- get_feature_genes(classifier,
                                      node = "root",
                                      db = org.Hs.eg.db)
    head(feature_genes)
    tic <- proc.time()
    cds <- classify_cells(cds, classifier,
                           db = org.Hs.eg.db,
                           cluster_extend = TRUE,
                           cds_gene_id_type = "ENSEMBL")
    toc <- proc.time()
    toc - tic
    
    colData(dlpfc_Anno)[paste0("garnett", id, "top", top)] <- pData(cds)[,5]
    }
  }
  tocAll <- Sys.time()
  print(tocAll - ticAll)

  outpath <- paste0("./Real_data_results/dataFiles/DLPFC/Brain12cross/")
  setwd(outpath)
  saveRDS(dlpfc_Anno, paste0(ID, "SCINA_Garnett.rds"))
}
