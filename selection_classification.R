library(caret)
library(nnet)

tmp.env <- new.env()
load("data/data_pbmc_for_project.RData", tmp.env)

# load("data/pbmc_pas_go.Rdata")
# load("data/pbmc_pas_reactome.Rdata")
# load("data/pbmc_pas_kegg.Rdata")
# load("data/pbmc_pas_sig.Rdata")
load("data/pbmc_pas_hallmark.RData")

meta <- get("meta", tmp.env)
label <- as.factor(meta$CellType)
uq_class_names <- unique(label)

pas_jasmine_lh <- as.data.frame(t(pas_jasmine_lh))
pas_gsva <- as.data.frame(pas_gsva)
pas_transformation <- list(pas_aucell, pas_bina, pas_cerno, pas_dropratio, pas_gsva, 
                           pas_jasmine_lh, pas_jasmine_or, pas_mean, pas_plage, 
                           pas_sipsic, pas_spca, pas_ssGSEA, pas_ucell, pas_vision, pas_zscore)

names_vec <- c("AUCell", "BINA", "Cerno", "DropRatio", "GSVA", "JasmineLH", 
               "JasmineOR", "Mean", "Plage",
               "SiPSiC", "sPCA", "ssGSEA", "UCell", "Vision", "Zscore")
names(pas_transformation) <- names_vec


pas_comparision <- list()
for(pas in 5:length(pas_transformation)){
  class_comparision <- list()
  for(class in 1:length(uq_class_names)){
    lab <- ifelse(label == uq_class_names[class], 1, 0)
    result_each_path_mlr <- lapply(rownames(pas_transformation[[pas]]), function(x){
      data <- as.data.frame(t(pas_transformation[[pas]][x,]))
      data$lab <- as.factor(lab)
      
      sample <- sample(c(TRUE, FALSE), nrow(data), replace = T, prob = c(.8, .2))
      train <- data[sample,]
      test <- as.data.frame(data[!sample,1])
      colnames(test) <- colnames(train)[1]

      lr_mdl<- glm(lab~., train, family = "binomial")
      pred <- predict(lr_mdl, test, type="response")
      pred_lab <- factor(ifelse(pred > .5, 1, 0), levels = c(0, 1))
      mtx <- confusionMatrix(pred_lab, data[!sample,2], positive = "1")
      return(mtx)
    })
    names(result_each_path_mlr) <- rownames(pas_transformation[[pas]])
    class_comparision[[class]] <- result_each_path_mlr
  }
  names(class_comparision) <- uq_class_names
  pas_comparision[[pas]] <- class_comparision
}
names(pas_comparision) <- names_vec
results_mlr_selection <- pas_comparision



