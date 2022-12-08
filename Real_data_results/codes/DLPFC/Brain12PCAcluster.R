##################### Supplementary Figure 2c.#####################
#######################################  
library(SingleCellExperiment) 
library(Seurat)
 

ARI <- matrix(0, nrow = 12, ncol = 3)

 
for (iter in (1:12)) {
  dir.exp <- "./Real_data_results/dataFiles/DLPFC/unsupervised/"
  dir.spark <- "./Real_data_analysis/DLPFC/Datasets/brain12_spark/"
  tic <- Sys.time()
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
  Adj_sp <- runAdj(pos, platform="Visium")
  #------------------------------------------
  #K <- trueK_set[iter]
  hq <- 15 # default as 15 factors
  set.seed(20132014)
  #----------------------------------------------------------prepare a marker list
  #----------------------------------------------------------
  marker_dir <- "./Real_data_analysis/DLPFC/Datasets/"
  marker_file <- paste0("151507", "_DEGmarkerTop_", 5, ".rds")
  markers <- readRDS(paste0(marker_dir, marker_file))

      if (ID %in% name_ID12[5:8]) {
        markers$Layer1 = NULL
        markers$Layer2 = NULL
      }  
      rho <- marker_list_to_mat(markers, TRUE)
      #rho <- rho[, which(colSums(rho) !=0 )]
      K <- ncol(rho)

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
      
   

## 1.  clustering on pca of marker genes
marker_logcount_pca <- prcomp(X_m)
marker_embedding = marker_logcount_pca$x[,1:15]
rownames(marker_embedding) = paste0("spot", 1:nrow(X_m))
colnames(marker_embedding) = paste0("gene", 1:15)
seu = CreateSeuratObject(
  counts = t(marker_embedding),
  project = "CreateSeuratObject",
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = NULL,
)
seu[["pca"]] <- CreateDimReducObject(embeddings = marker_embedding, key = "PC_", assay = DefaultAssay(seu))
seu <- FindNeighbors(seu,reduction = "pca", dims = 1:15)
seu = FindClusters(seu, resolution = 1)
ARI[iter, 1] = mclust::adjustedRandIndex(seu$seurat_clusters, y)


## 2. clustering on embeddings of non-marker (spatialAnno)
ID <- name_ID12[iter] 
fit_s <- readRDS(paste0("./Real_data_results/dataFiles/DLPFC/unsupervised/", ID, "_icmem.rds"))
hZ <- fit_s$Ez_u
rownames(hZ) = paste0("spot", 1:nrow(X_m))
colnames(hZ) = paste0("gene", 1:ncol(hZ)) 
seu = CreateSeuratObject(
  counts = t(hZ),
  project = "CreateSeuratObject",
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = NULL,
)
seu[["pca"]] <- CreateDimReducObject(embeddings = hZ, key = "PC_", assay = DefaultAssay(seu))
seu <- FindNeighbors(seu,reduction = "pca", dims = 1:15)
seu = FindClusters(seu, resolution = 1)
ARI[iter, 2] = mclust::adjustedRandIndex(seu$seurat_clusters, y)


## 3.  clustering on embeddings of marker and non-marker (spatialAnno)
PC <- cbind(marker_embedding, hZ)
rownames(PC) = paste0("spot", 1:nrow(X_m))
colnames(PC) = paste0("gene", 1:ncol(PC)) 

seu = CreateSeuratObject(
  counts = t(PC),
  project = "CreateSeuratObject",
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = NULL,
)
seu[["pca"]] <- CreateDimReducObject(embeddings = PC, key = "PC_", assay = DefaultAssay(seu))
seu <- FindNeighbors(seu,reduction = "pca", dims = 1:30)
seu = FindClusters(seu, resolution = 1)
ARI[iter, 3] = mclust::adjustedRandIndex(seu$seurat_clusters, y)

}

saveRDS(ARI, file="embeddingsARI.rds")



############################ Supplementary Figure 2c.#####################

library(ggplot2)
library(ggsci)
library(patchwork)
library(systemfonts)
theme_Publication <- function(base_size=14, base_family="Helvetica") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               #axis.text.x = element_blank(),
               axis.line = element_line(colour="black"),
               axis.ticks.x= element_blank(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.text = element_text(size=18),
               legend.position = "none",
               legend.direction = "horizontal",
               legend.key.size= unit(0.4, "cm"),
               legend.spacing  = unit(0, "cm"),
               legend.title = element_blank(),
               plot.margin=unit(c(10,5,5,0),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))
      
}
library(reshape)
library(ggplot2)
library(ggsci)
ARI <- readRDS("embeddingsARI.rds")
colnames(ARI) <- c("Marker", "Nonmarker", "Marker+\n nonmarker")

df <- melt(data=ARI)
df$X2 <- factor(df$X2, levels=colnames(ARI)) 

library(rstatix) 
library(ggpubr)

stat.test <- df %>% wilcox_test(value ~ X2) 
stat.test <- stat.test %>% add_xy_position(x="X2")


p <- ggplot(data=df, aes(x = X2, y = value)) +   
 geom_boxplot(fill=pal_simpsons()(3)) +
  stat_pvalue_manual(stat.test, label = "p") +
 coord_cartesian(ylim=c(0, .6)) +  theme_Publication(base_size=24) + 
   scale_fill_simpsons()    
   
pdf("dlpfc_embedARI.pdf") #保存为pdf
p
dev.off()


