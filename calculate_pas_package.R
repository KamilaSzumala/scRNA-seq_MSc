library(Seurat)
library(Jasmine)
library(matrixStats)

load("data/data_pbmc_for_project.RData")
load("paths/reactome_all_path.RData")

paths <- reactome_path$gene_sets

# ---- preprocessing ----
cells.Seurat <- NormalizeData(cells.Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
cells.Seurat <- FindVariableFeatures(cells.Seurat, selection.method = "vst", nfeatures = round(dim(cells.Seurat)[1]*0.2))

data <- as.data.frame(GetAssayData(object = cells.Seurat, slot = "data"))
rm(cells, cells.Seurat, genes)

# ---- filter out pathways ---- 


filter_minmax <- function (pathway, filt_min, filt_max) {
  name <- c()
  for (i in 1:length(pathway)) {
    siz <- length(pathway[[i]])
    if (siz < filt_min | siz > filt_max) {
      name <- append(name, names(pathway[i]))
    }
  }
  pathway[name] = NULL
  return(pathway)
}


filter_cover <- function (X, pathway, filt_cov) {
  name <- c()
  X <- as.matrix(X)
  for (i in seq_along(pathway)) {
    df <- X[rownames(X) %in% pathway[[i]], ]
    if (is.matrix(df) && nrow(df) != 0) {
      cof <- round(nrow(df)/length(pathway[[i]]), digits = 2)
      if (cof < filt_cov) {
        name <- append(name, names(pathway[i]))
      }
    }
    else if (is.vector(df)) {
      cof <- round(length(df)/length(pathway[[i]]), digits = 2)
      if (cof < filt_cov) {
        name <- append(name, names(pathway[i]))
      }
    }
    else {
      name <- c(name, names(pathway[i]))
    }
  }
  pathway[name] = NULL
  return(pathway)
}


pathways <- filter_minmacells(paths, 15, 500)
pathways <- filter_cover(data, pathways, .65)

# ---- calculate PAS ----
extract_pathway <- function (cells, path) {
  path <- unlist(path)
  data_path <- cells[rownames(cells) %in% path, ]
  tmp <- rowVars(as.matrix(data_path))
  data_path <- data_path[!(is.na(tmp) == T | tmp == 0), ]
  return(data_path)
}


#* ---- Mean -----
pas_mean <- as.data.frame(do.call(rbind, lapply(pathways, 
                                                 function(path) {
                                                   colMeans(extract_pathway(data, path))
                                                 })))
pas_mean <- as.data.frame(pas_mean)
rownames(pas_mean) <- names(pathways)

#* ---- BINA ----
pas_bina <- do.call(rbind, lapply(pathways, function(path) {
  df <- extract_pathway(data, path)
  row <- colSums((df != 0)/nrow(df))
  row_logit <- log((row + 0.1)/(1 - row + 0.1))
  return(row_logit)
}))

rownames(pas_bina) <- names(pathways)
colnames(pas_bina) <- colnames(data)

#* ---- Zscore ----
X_normed <- scale_zscore(data)

pas_zscore <- do.call(rbind, lapply(pathways, function(path) {
  df_path <- extract_pathway(X_normed, path)
  row_score <- colSums(df_path)/sqrt(nrow(df_path))
  return(row_score)
}))

pas_zscore <- as.data.frame(pas_zscore)
rownames(pas_zscore) <- names(pathways)

#* ---- CERNO ----
X_ranked <- Rank_data(data)

pas_cerno <- do.call(rbind, lapply(pathways, function(path) {
  df_path <- extract_pathway(X_ranked, path)
  row_AUC <- Calc_AUC(X_ranked, df_path)
  return(row_AUC)
}))

pas_cerno <- as.data.frame(pas_cerno)
rownames(pas_cerno) <- names(pathways)


#* ---- AUCell ----
cells_rankings <- AUCell_buildRankings(as.matrix(data))
AUCs <- AUCell_calcAUC(pathways, cells_rankings, aucMaxRank = 0.25 * 
                         nrow(data), nCores = 1)
pas_aucell <- as.data.frame(AUCs@assays@data@listData$AUC)


#* ---- JASMINE ---- 
pas_jasmine_or <- as.data.frame(sapply(pathways, 
                                  function(path) {
                                    result <- as.vector(JASMINE(data, path, "oddsratio"))
                                    return(result$JAS_Scores)
                                    }))

pas_jasmine_or <- as.data.frame(t(pas_jasmine_or)) 

rownames(pas_jasmine_or) <- names(pathways)
colnames(pas_jasmine_or) <- colnames(data)

#*
pas_jasmine_lh <- as.data.frame(sapply(pathways, 
                                       function(path) {
                                         result <- as.vector(JASMINE(data, path, "likelihood"))
                                         return(result$JAS_Scores)
                                       }))

pas_jasmine_lh <- as.data.frame(t(pas_jasmine_lh)) 

rownames(pas_jasmine_lh) <- names(pathways)
colnames(pas_jasmine_lh) <- colnames(data)

#* ---- GSVA ----
library(GSVA)
pas_gsva <- gsva(as.matrix(data), pathways, kcdf = "Gaussian", mx.diff = F)
pas_gsva <- as.data.frame(pas_gsva)

#* ---- ssGSEA ----
pas_ssGSEA <- ssgseaParam(as.matrix(data), pathways)
pas_ssGSEA <- as.data.frame(gsva(pas_ssGSEA))
rownames(pas_ssGSEA) <- names(pathways)

#* ---- Sparse PCA ---
library(sparsepca)
library(matrixStats)

Path_df <- function(pathway, cells){
  x <- cells[rownames(cells) %in% pathway,]
  tmp <- rowVars(as.matrix(x))
  
  m <- x[!(is.na(tmp)==T | tmp == 0),]
  return(m)
}

calSpca <- function(pathways, cells){
  dfSpca <- data.frame()
  for (i in 1:length(pathways)){
    df <- extract_pathway(cells, pathways[[i]])
    pca.path <- spca(t(df), scale = T, center = T, max_iter = 800)
    pca.path.x <- as.data.frame(pca.path$scores)
    dfSpca <- rbind(dfSpca, t(pca.path.x$V1))
  }
  rownames(dfSpca) <- names(pathways)
  return (dfSpca)
}

pas_spca <- calSpca(pathways, data)
colnames(pas_spca) <- colnames(data)

#* ---- Plage ---- 
calPlage <- function(pathways, cells){
  dfPlage <- data.frame(); 
  for (i in 1:length(pathways)){
    df <- extract_pathway(cells, pathways[[i]])
    pca.path <- prcomp(t(df), scale = T, center = T)
    pca.path.x <- as.data.frame(pca.path$x)
    dfPlage <- rbind(dfPlage, pca.path.x$PC1)
  }
  rownames(dfPlage) <- names(pathways) 
  return (dfPlage)
}

pas_plage <- calPlage(pathways, data)
colnames(pas_plage) <- colnames(data)

#* ---- DropRatio ----
calDropRatio <- function(pathways, cells){
  dfDropRatio <- data.frame()
  for (i in 1:length(pathways)){
    df <- extract_pathway(cells, pathways[[i]])
    row <- colSums(((df!=0)/nrow(df)))
    dfDropRatio <- rbind(dfDropRatio, row)
  } 
  rownames(dfDropRatio) <- names(pathways)
  return (dfDropRatio)
}

pas_dropratio <- calDropRatio(pathways, data)
colnames(pas_dropratio) <- colnames(data)

#* ---- VISION ----
library(VISION)

n.umi <- colSums(cells)

scaled_counts <- t(t(cells) / n.umi) * mean(as.matrix(cells))

signatures <- list()
for(i in 1:length(pathways)){signatures[[i]] <- createGeneSignature(names(pathways[i]), sapply(pathways[[i]], function(y) y = 1))}
vis <- Vision(scaled_counts, signatures = signatures)
vis <- analyze(vis)


pas_vision <- as.data.frame(t(vis@SigScores))

#* ---- SiPSiC -----
library(SiPSiC)

pathScores <- lapply(pathways, function(x) getPathwayScores(as.matrix(raw_data), x))
pas_sipsic <- as.data.frame(t(sapply(pathScores, function(x) x$pathwayScores)))

#* ---- UCell ----
library(UCell)

scores <- ScoreSignatures_UCell(data, pathways)
pas_ucell <- as.data.frame(t(scores))

