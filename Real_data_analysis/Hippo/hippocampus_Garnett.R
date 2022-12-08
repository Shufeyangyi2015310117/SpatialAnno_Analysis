library(garnett)
 
setwd("./Real_data_analysis/Hippo/Datasets/")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
#slide.seq <- LoadData("ssHippo")
seulist <- readRDS("seulist_Hip2_minFeature15_preprocess.rds")
embro <- 1
slide.seq <- seulist[[embro]]
 
 
slide.seq <- SCTransform(slide.seq, assay = "RNA", return.only.var.genes = FALSE, verbose = FALSE) 
mat <- slide.seq@assays$RNA@counts#MOB_raw
rownames(mat) <- rownames(mat)

fdata <- data.frame(gene_short_name = rownames(mat), num_cells_expressed = rowSums(mat != 0))
pdata <- as.data.frame(matrix(1, ncol(mat), 2))
rownames(pdata) = colnames(mat)


# create a new CDS object
pd <- new("AnnotatedDataFrame", data = pdata)
fd <- new("AnnotatedDataFrame", data = fdata)
hp_cds <- newCellDataSet(as(mat, "dgCMatrix"),
                           phenoData = pd,
                           featureData = fd)


# generate size factors for normalization later
hp_cds <- estimateSizeFactors(hp_cds)


#data = ""
#write.table(data, file = "/Users/yiyang/Desktop/2020-9-19/2022/Xingjie/MOB2/MOB_markers.txt")

iter <- 1
dataSource <- 1
mk_n <-  c(5, 10)
 
markers <-  readRDS("./hippo_markers.rds")
markers2 <- markers
sink("hippo_markers.txt")
for (i in 1:length(markers2)){
  cat(paste0(">",names(markers2)[i],"\n"))
  if(names(markers2)[i] == "TH") cat(paste0("expressed: ", paste(markers2[[i]][1], markers2[[i]][2], markers2[[i]][3],
              sep = ", ")))
  if(names(markers2)[i] != "TH") cat(paste0("expressed: ", paste(markers2[[i]][1], markers2[[i]][2], markers2[[i]][3],
             markers2[[i]][4], markers2[[i]][5], sep = ", ")))
  cat("\n")
  cat("\n")
}
sink()
 

#file.show("hippo_markers.txt")

library(org.Mm.eg.db)
 
library(garnett)

marker_file_path <- "hippo_markers.txt"
marker_check <- check_markers(hp_cds, marker_file_path,
                              db=org.Mm.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")

plot_markers(marker_check)
dev.off()

 new_g <- garnett:::convert_gene_ids(row.names(fdata),
                            db=org.Mm.eg.db,
                            "SYMBOL",
                            "SYMBOL")



library(org.Mm.eg.db)
set.seed(260)
tic <- proc.time()
marker_file_path <- "hippo_markers.txt"
pbmc_classifier <- train_cell_classifier(cds = hp_cds,
                                         marker_file = marker_file_path,
                                         db=org.Mm.eg.db,
                                         cds_gene_id_type = "SYMBOL",
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "SYMBOL")



feature_genes <- get_feature_genes(pbmc_classifier,
                                   node = "root",
                                   db = org.Mm.eg.db)
head(feature_genes)


hp_cds <- classify_cells(hp_cds, pbmc_classifier,
                           db = org.Mm.eg.db,
                           cluster_extend = TRUE,
                           cds_gene_id_type = "SYMBOL")
toc <- proc.time()
toc - tic

head(pData(hp_cds))
table(pData(hp_cds)[,6])
table(pData(hp_cds)[,5])
hp_cds
saveRDS(hp_cds, file = paste0("hippo", embro, "_garnett.rds"))

 

Pre_cell_type = pData(hp_cds)[,5] 
saveRDS(Pre_cell_type, file = paste0("hippo", embro, "_anno_garnett.rds"))

