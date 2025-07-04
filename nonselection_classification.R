library(caret)
library(nnet)

tmp.env <- new.env()
load("data/data_pbmc_for_project.RData", tmp.env)

load("data/pbmc_pas_go.Rdata")
# load("data/pbmc_pas_reactome.Rdata")
# load("data/pbmc_pas_kegg.Rdata")
# load("data/pbmc_pas_sig.Rdata")
# load("data/pbmc_pas_hallmark.RData")

meta <- get("meta", tmp.env)

pas_jasmine_lh <- as.data.frame(t(pas_jasmine_lh))
pas_transformation <- list(pas_aucell, pas_bina, pas_cerno, pas_dropratio, pas_gsva, 
                           pas_jasmine_lh, pas_jasmine_or, pas_mean, pas_plage, 
                           pas_sipsic, pas_spca, pas_ssGSEA, pas_ucell, pas_vision, pas_zscore)

names_vec <- c("AUCell", "BINA", "Cerno", "DropRatio", "GSVA", "JasmineLH", 
               "JasmineOR", "Mean", "Plage",
               "SiPSiC", "sPCA", "ssGSEA", "UCell", "Vision", "Zscore")

custom_summary <- function(data, lev = NULL, model = NULL) {
  cm <- confusionMatrix(data$pred, data$obs)
  balanced_accuracy <- cm[["byClass"]][,"Balanced Accuracy"]
  specificity <- cm[["byClass"]][,"Specificity"]
  sensitivity <- cm[["byClass"]][,"Sensitivity"]
  
  max <- mean(balanced_accuracy)
  names(max) <- "Balanced Accuracy"
  max
}

train_control <- trainControl(method = "cv", number = 10, summaryFunction = custom_summary, savePredictions = T)
tune_grid <- expand.grid(alpha = 0.5, lambda = 0.03)

result_noselection_mlr <- lapply(pas_transformation, function(x){
                              data <- as.data.frame(t(x))
                              data$lab <- as.factor(as.numeric(as.factor(meta$CellType)))
                              model <- train(lab ~ ., 
                                             data = data, 
                                             method = "glmnet",
                                             trace = F,
                                             metric = "Balanced Accuracy",
                                             trControl = train_control,
                                             tuneGrid = tune_grid)
                                return(model)
                          })
names(result_noselection_mlr) <- names_vec
  


result_noselection_rf <- lapply(pas_transformation, function(x){
                              data <- as.data.frame(t(x))
                              data$lab <- as.factor(as.numeric(as.factor(meta$CellType)))
                              model <- train(lab ~ ., 
                                             data = data, 
                                             method = "rf",
                                             trace = F,
                                             metric = "Balanced Accuracy",
                                             trControl = train_control
                                             )
                              return(model)
                            })

names(result_noselection_rf) <- names_vec

