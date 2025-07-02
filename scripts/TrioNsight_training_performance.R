
if (!requireNamespace("caret", quietly = TRUE)) {
  install.packages('caret')   
}
library(caret)
if (!requireNamespace("data.table", quietly = TRUE)) {
  BiocManager::install('data.table')
}
library(data.table)
if (!requireNamespace("utils", quietly = TRUE)) {
  BiocManager::install('utils')
}
library(utils)
if (!requireNamespace("randomForest", quietly = TRUE)) {
  install.packages('randomForest')   
}
library(randomForest)
if (!requireNamespace("kernlab", quietly = TRUE)) {
  BiocManager::install('kernlab')
}
library(kernlab)
if (!requireNamespace("gbm", quietly = TRUE)) {
  BiocManager::install('gbm')
}
library(gbm)
if (!requireNamespace("glmnet", quietly = TRUE)) {
  install.packages('glmnet')   
}
library(glmnet)
if (!requireNamespace("mboost", quietly = TRUE)) {
  BiocManager::install('mboost')
}
library(mboost)
if (!requireNamespace("naivebayes", quietly = TRUE)) {
  BiocManager::install('naivebayes')
}
library(naivebayes)
if (!requireNamespace("xgboost", quietly = TRUE)) {
  install.packages('xgboost')   
}
library(xgboost)
if (!requireNamespace("ranger", quietly = TRUE)) {
  BiocManager::install('ranger')
}
library(ranger)
if (!requireNamespace("earth", quietly = TRUE)) {
  BiocManager::install('earth')
}
library(earth)
if (!requireNamespace("C50", quietly = TRUE)) {
  install.packages('C50')   
}
library(C50)
if (!requireNamespace("mda", quietly = TRUE)) {
  BiocManager::install('mda')
}
library(mda)
if (!requireNamespace("adabag", quietly = TRUE)) {
  BiocManager::install('adabag')
}
library(adabag)
if (!requireNamespace("flextable", quietly = TRUE)) {
  BiocManager::install('flextable')
}
library(flextable)
if (!requireNamespace("magrittr", quietly = TRUE)) {
  BiocManager::install('magrittr')
}
library(magrittr)
if (!requireNamespace("dplyr", quietly = TRUE)) {
  BiocManager::install('dplyr')
}
library(dplyr)
if (!requireNamespace("tibble", quietly = TRUE)) {
  BiocManager::install('tibble')
}
library(tibble)

preprocess_training_data <- function(dt){
  
  if (all(c("REF", "ALT") %in% names(dt))) {
    setnames(dt, old = c("REF", "ALT"), new = c("wt", "mt"))
  }
  
  suppressWarnings(dt[, variant_by_inheritance := NULL])
  suppressWarnings(dt[, Accession := NULL])
  suppressWarnings(dt[, Source := NULL])
  suppressWarnings(dt[, CHROM_POS_REF_ALT := NULL])
  suppressWarnings(dt[, Variant_Seq := NULL])
  suppressWarnings(dt[, wt := NULL])
  suppressWarnings(dt[, mt := NULL])
  suppressWarnings(dt[, PolyPhen_Class := NULL])
  suppressWarnings(dt[, SIFT_Class := NULL])
  suppressWarnings(dt[, PolyPhen_Score := NULL])
  suppressWarnings(dt[, SIFT_Score := NULL])
  suppressWarnings(dt[, ESM_Score := NULL])
  ## suppressWarnings(dt[, PHACT_Score := NULL])
  suppressWarnings(dt[, AlphaMissense_Class := NULL])
  suppressWarnings(dt[, Variant_Identifier := NULL])
  suppressWarnings(dt[, Variant_Location := NULL])
  suppressWarnings(dt[, REF_1letter := NULL])
  suppressWarnings(dt[, ALT_1letter := NULL])
  suppressWarnings(dt[, REF_3letter := NULL])
  suppressWarnings(dt[, ALT_3letter := NULL])
  suppressWarnings(dt[, nativeSeqCheck := NULL])
  suppressWarnings(dt[, Model := NULL])
  suppressWarnings(dt[, Inherited_Location := NULL])
  
  # Function to replace NA values with mean if column is numeric
  replace_na_with_mean <- function(dt, col) {
    if (is.numeric(dt[[col]]) && anyNA(dt[[col]])) {
      # Calculate mean excluding NA values
      col_mean <- mean(dt[[col]], na.rm = TRUE)
      
      # Replace NA values with the mean
      dt[is.na(dt[[col]]), (col) := col_mean]
      
      # Print the action taken
      cat("Replaced NA values in column", col, "with the mean:", col_mean, "\n")
    } else if (!is.numeric(dt[[col]]) && anyNA(dt[[col]])) {
      cat("Skipped column", col, "because it is not numeric\n")
    }
  }
  # Iterate over columns and replace NA values with mean if column is numeric and has NA values
  for (col in names(dt)) {
    replace_na_with_mean(dt, col)
  }
  
  return(dt)
}

initialize_formulas <- function(dt){
  
  variant_clusters_names <- names(dt)[grepl('distance_to_cluster', names(dt))]
  variant_clusters_columns <- variant_clusters_names
  variant_clusters_names <- paste(variant_clusters_names, collapse = " - ")
  
  StructuralFeatures <- c("Interface", "Relaccess", "Impact",	"HBonds",	"SPhobic", "CPhilic", "BCharge", "SSGeom", "Voids",	"MLargest1",
                          "MLargest2", "MLargest3",	"MLargest4", "MLargest5", "MLargest6", "MLargest7", "MLargest8", "MLargest9", 
                          "MLargest10", "NLargest1", "NLargest2", "NLargest3", "NLargest4", "NLargest5", "NLargest6", "NLargest7", 
                          "NLargest8", "NLargest9", "NLargest10",	"Clash", "Glycine", "Proline", "CisPro", "foldxDdg", "foldx_plddt")
  
  StructuralFeatures <- intersect(names(dt), StructuralFeatures)
  StructuralFeatures <- c(StructuralFeatures, variant_clusters_columns)
  StructuralFeatures_columns <- StructuralFeatures
  StructuralFeatures <- paste(StructuralFeatures, collapse = " - ")
  
  SAAPdap <- c("Binding", "SProtFT_NAMES", "SProtFT_BOOL", "SProtFT_FEATURES", "SProtFTDist_BOOL", 
               "SProtFT_ACT_SITE", "SProtFT_BINDING", "SProtFT_CA_BIND", "SProtFT_DNA_BIND", "SProtFT_NP_BIND", 
               "SProtFT_METAL", "SProtFT_MOD_RES", "SProtFT_CARBOHYD", "SProtFT_MOTIF", "SProtFT_LIPID", 
               "Interface", "Relaccess", "Impact",	"HBonds",	"SPhobic", "CPhilic", "BCharge", "SSGeom", "Voids",	"MLargest1",
               "MLargest2", "MLargest3",	"MLargest4", "MLargest5", "MLargest6", "MLargest7", "MLargest8", "MLargest9", 
               "MLargest10", "NLargest1", "NLargest2", "NLargest3", "NLargest4", "NLargest5", "NLargest6", "NLargest7", 
               "NLargest8", "NLargest9", "NLargest10",	"Clash", "Glycine", "Proline", "CisPro")
  
  SAAPdap <- intersect(names(dt), SAAPdap)
  SAAPdap_columns <- SAAPdap
  SAAPdap <- paste(SAAPdap, collapse = " - ")
  
  SequenceFeatures <- c("Binding", "SProtFT_NAMES", "SProtFT_BOOL", "SProtFT_FEATURES", "SProtFTDist_BOOL", 
                        "SProtFT_ACT_SITE", "SProtFT_BINDING", "SProtFT_CA_BIND", "SProtFT_DNA_BIND", "SProtFT_NP_BIND", 
                        "SProtFT_METAL", "SProtFT_MOD_RES", "SProtFT_CARBOHYD", "SProtFT_MOTIF", "SProtFT_LIPID", 
                        'ScoreCons', 'BLOSUM80', 'Grantham', 'varSite')
  SequenceFeatures <- intersect(names(dt), SequenceFeatures)
  SequenceFeatures_columns <- SequenceFeatures
  SequenceFeatures <- paste(SequenceFeatures, collapse = " - ")
  
  CompPredictorsOnly <- c(SAAPdap_columns, variant_clusters_columns, 'ScoreCons', 'BLOSUM80', 'Grantham', 'varSite', 'foldxDdg', 'foldx_plddt')
  CompPredictorsOnly <- intersect(names(dt), CompPredictorsOnly)
  CompPredictorsOnly <- paste(CompPredictorsOnly, collapse = " - ")
  
  return(list(variant_clusters_names = variant_clusters_names, 
              StructuralFeatures = StructuralFeatures, 
              SAAPdap = SAAPdap, 
              SequenceFeatures = SequenceFeatures, 
              CompPredictorsOnly = CompPredictorsOnly))
  
}

train_all_methods_all_formulas <- function(dt, models, iterations, preTuned = FALSE, model_output){
  
  dt_save <- dt
  
  # Initialize .folds_
  num_folds <- length(unique(dt$MSALocation))
  fold_cols <- paste(".folds_", seq_len(num_folds), sep = "")
  
  idxList <- list()
  # LOOCV based on MSALocation
  unique_locs <- unique(dt$MSALocation)
  dt <- dt_save
  
  for (i in 1:length(unique_locs)) {
    idxList[[fold_cols[i]]] <- which(dt$MSALocation != unique_locs[i])
  }
  fitControl <- trainControl(index = idxList, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE, verboseIter = TRUE)
  
  trained_models <- list()
  
  formula_list <- initialize_formulas(dt = dt)
  ## dt[, rowIndex := NULL]
  
  for(model in models){
    
    for(i in 1:length(iterations)) {
      if(i == 1){
        formulaName <- 'Inclusive'
      }
      else if (i == 2){
        formulaName <- 'NofoldX'
      }
      else if (i == 3){
        formulaName <- 'NoVariantClusters'
      }
      else if (i == 4){
        formulaName <- 'NoSAAPdap'
      }
      else if (i == 5){
        formulaName <- 'NoStructuralFeatures'
      }
      else if (i == 6){
        formulaName <- 'NovarSite'
      }
      else if (i == 7){
        formulaName <- 'NoSequenceFeatures'
      }
      else if (i == 8){
        formulaName <- 'NoAlphaMissense'
      }
      else if (i == 9){
        formulaName <- 'NoVariPred'
      }
      else if (i == 12){
        formulaName <- 'NoPHACT'
      }
      else if (i == 11){
        formulaName <- 'NoPHACTboost'
      }
      else if (i == 12){
        formulaName <- 'Exclusionary'
      }
      else if (i == 13){
        formulaName <- 'CompPredictorsOnly'
      }
      
      cat(paste('Model: ', model, ', Iteration: ', formulaName, sep = ''), '\n')
      
      formula <- as.formula("Class ~ . - target_id - MSALocation - Inherited_Seq")
      
      if(formulaName != 'Inclusive'){
        if(formulaName == 'NoVariantClusters'){
          formula <- as.formula(paste(deparse(formula), "-", formula_list$variant_clusters_names))
        }
        else if(formulaName == 'NoSAAPdap'){
          formula <- as.formula(paste(deparse(formula), "-", formula_list$SAAPdap))
        }
        else if(formulaName == 'NoStructuralFeatures'){
          formula <- as.formula(paste(deparse(formula), "-", formula_list$StructuralFeatures))
        }
        else if(formulaName == 'NoSequenceFeatures'){
          formula <- as.formula(paste(deparse(formula), "-", formula_list$SequenceFeatures))
        }
        else if(formulaName == 'CompPredictorsOnly'){
          formula <- as.formula(paste(deparse(formula), "-", formula_list$CompPredictorsOnly))
        }
        else{
          formula <- as.formula(paste(deparse(formula), "-", iterations[i]))
        }
      }
      cat(paste('The formula is:', Reduce(paste, deparse(formula)), sep = ' '), '\n')
      
      modelsWithLessParams <- c('knn', 'gbm', 'xgbTree', 'xgbLinear', 'glmboost', 'ranger', 'earth', 'rpart', 'bagFDA')
      if(isFALSE(preTuned)){
        if(model %in% modelsWithLessParams){
          modelTrained <- train(formula, data = dt, 
                                method = model,
                                preProc = c("center", "scale"),
                                trControl = fitControl,
                                tuneLength = 15,  # Add tuning
                                metric = "ROC",  # Change this if not a classification problem
                                maximize = TRUE)  # TRUE if using accuracy or other metrics where higher is better
        }
        else{
          modelTrained <- train(formula, data = dt, 
                                method = model,
                                preProc = c("center", "scale"),
                                trControl = fitControl,
                                allowParallel = TRUE,
                                tuneLength = 15,  # Add tuning
                                metric = "ROC",  # Change this if not a classification problem
                                maximize = TRUE,  # TRUE if using accuracy or other metrics where higher is better
                                verbose = TRUE)
        }
      }
      else{
        preTuned <- preTuned_dt[[model]]
        preTuned <- expand.grid(preTuned)
        if(model %in% modelsWithLessParams){
          modelTrained <- train(formula, data = dt, 
                                method = model,
                                preProc = c("center", "scale"),
                                trControl = fitControl,
                                tuneGrid = preTuned,
                                metric = "ROC",  # Change this if not a classification problem
                                maximize = TRUE)  # TRUE if using accuracy or other metrics where higher is better
        }
        else{
          modelTrained <- train(formula, data = dt, 
                                method = model,
                                preProc = c("center", "scale"),
                                trControl = fitControl,
                                allowParallel = TRUE,
                                tuneGrid = preTuned,
                                metric = "ROC",  # Change this if not a classification problem
                                maximize = TRUE,  # TRUE if using accuracy or other metrics where higher is better
                                verbose = TRUE)
        }
      }
      trained_models[[model]][[formulaName]] <- modelTrained
    }
  }
  saveRDS(trained_models, file = paste('predictor/', model_output, '.rds', sep = ''))
  return(trained_models)
}

optimizing_hyperparameters_retrieve_predictions <- function(dt, model_input, model_output){
  
  predictions_models <- list()
  
  for(i in 1:length(model_input)){
    for(j in 1:length(model_input[[i]])){
      
      modelName <- names(model_input)[[i]]
      formulaName <- names(model_input[[i]])[[j]]
      
      tmpModel <- model_input[[modelName]][[formulaName]]
      tmpPred <- as.data.table(tmpModel$pred)
      
      bestTune_params <- tmpModel$bestTune
      
      cat(paste('Model:', modelName, 'Formula:', formulaName, '\n'))
      cat(paste('There are', length(bestTune_params), 'bestTune parameters\n'))
      cat(paste('Best tune parameters:', names(bestTune_params), '\n'))
      
      for(bestTune_param in names(bestTune_params)){
        tmpPred <- tmpPred[which(tmpPred[[bestTune_param]] == tmpModel$bestTune[[bestTune_param]]),]
        tmpPred[[bestTune_param]] <- NULL
      }
      tmpPred[['Resample']] <- NULL
      predictions_models[[modelName]][[formulaName]] <- tmpPred
    }
  }
  
  dt[, rowIndex := .I]
  for(i in 1:length(predictions_models)){
    for(j in 1:length(predictions_models[[i]])){
      predictions_models[[i]][[j]] <- as.data.table(merge(predictions_models[[i]][[j]], dt[,.(target_id, rowIndex)], by = 'rowIndex'))
      predictions_models[[i]][[j]][, rowIndex := NULL]
      predictions_models[[i]][[j]][, pred := NULL]
      
      observed <- unique(predictions_models[[i]][[j]][, .(target_id, obs)])
      colnames(observed)[2] <- 'obs'
      predictions_models[[i]][[j]][, obs := NULL]
      
      if(length(predictions_models[[i]][[j]][is.nan(Benign)]$Benign) > 0){
        cat(paste('There are', length(predictions_models[[i]][[j]][is.nan(Benign)]$Benign), 'NAs in Benign\n', sep = ' '))
        predictions_models[[i]][[j]][is.nan(Benign)]$Benign <- 0.50001
      }
      if(length(predictions_models[[i]][[j]][is.nan(Pathogenic)]$Pathogenic) > 0){
        cat(paste('There are', length(predictions_models[[i]][[j]][is.nan(Pathogenic)]$Pathogenic), 'NAs in Pathogenic\n', sep = ' '))
        predictions_models[[i]][[j]][is.nan(Pathogenic)]$Pathogenic <- 0.49999
      }
      
      predictions_models[[i]][[j]] <- aggregate(predictions_models[[i]][[j]], by = list(predictions_models[[i]][[j]]$target_id), mean)
      predictions_models[[i]][[j]] <- as.data.table(predictions_models[[i]][[j]])
      predictions_models[[i]][[j]][, target_id := NULL]
      colnames(predictions_models[[i]][[j]])[1] <- 'target_id'
      
      predictions_models[[i]][[j]][, pred := ifelse(Pathogenic > Benign, "Pathogenic", "Benign")]
      predictions_models[[i]][[j]] <- merge(predictions_models[[i]][[j]], observed, by = 'target_id')
    }
  }
  saveRDS(predictions_models, file = paste('predictor/predictions_', model_output, '.rds', sep = ''))
  return(predictions_models)
}

model_performance_calculation <- function(predictions_models, model_names, external_predictors = NULL, external_predictor_results = NULL, data_overlap = NULL, performance_output){
  
  formulas <- names(predictions_models[[1]])
  
  if(!is.null(external_predictors)){
    models_formulas <- c(external_predictors, c(outer(model_names, formulas, FUN = paste, sep = "~")))
  }
  else{
    models_formulas <- c(outer(model_names, formulas, FUN = paste, sep = "~"))
  }
  
  # Excludes Class column and TrioNSight column
  metrics <- data.table(Predictor = models_formulas)
  
  for(model in model_names){
    # Calculate metrics for each model + formula
    
    for(formula in formulas) {
      
      if(formula %in% names(predictions_models[[model]])) {
        
        data <- predictions_models[[model]][[formula]]$obs
        data <- as.data.table(as.character(data))
        colnames(data) <- 'Class'
        data$Class <- factor(data$Class, levels = c("Benign", "Pathogenic"))
        predicted <- predictions_models[[model]][[formula]]$pred
        predicted <- factor(predicted, levels = levels(data$Class))
        model_formula <- paste(model, '~', formula, sep = '')
        
        metrics[Predictor == model_formula, Precision := precision(predicted, data$Class, positive_class = "Pathogenic")]
        metrics[Predictor == model_formula, Recall := recall(predicted, data$Class, positive_class = "Pathogenic")]
        metrics[Predictor == model_formula, Accuracy := accuracy(predicted, data$Class)]
        metrics[Predictor == model_formula, False_Positive_Rate := false_positive_rate(predicted, data$Class, positive_class = "Pathogenic")]
        metrics[Predictor == model_formula, False_Negative_Rate := false_negative_rate(predicted, data$Class, positive_class = "Pathogenic")]
        metrics[Predictor == model_formula, F1_Score := f1_score(predicted, data$Class, positive_class = "Pathogenic")]
        metrics[Predictor == model_formula, MCC := mcc_score(predicted, data$Class, positive_class = "Pathogenic")]
      }
    }
  }
  
  if(!is.null(external_predictors)){
    external_predictor_results[Class == 'Ambigious', Class := 'Pathogenic']
    
    for(predictor in external_predictors) {
      
      data <- external_predictor_results
      data <- as.data.table(cbind(data[,.(target_id, Class)], data[[predictor]]))
      colnames(data)[3] <- predictor
      if(!is.null(data_overlap)){
        data <- data[!target_id %in% data_overlap[[predictor]]]
      }
      predicted <- data[[predictor]]
      observed <- data[['Class']]
      
      metrics[Predictor == predictor, Precision := precision(predicted, observed, positive_class = "Pathogenic")]
      metrics[Predictor == predictor, Recall := recall(predicted, observed, positive_class = "Pathogenic")]
      metrics[Predictor == predictor, Accuracy := accuracy(predicted, observed)]
      metrics[Predictor == predictor, False_Positive_Rate := false_positive_rate(predicted, observed, positive_class = "Pathogenic")]
      metrics[Predictor == predictor, False_Negative_Rate := false_negative_rate(predicted, observed, positive_class = "Pathogenic")]
      metrics[Predictor == predictor, F1_Score := f1_score(predicted, observed, positive_class = "Pathogenic")]
      metrics[Predictor == predictor, MCC := mcc_score(predicted, observed, positive_class = "Pathogenic")]
    }
    
    metrics[Predictor == 'PolyPhen2', Predictor := 'PolyPhen-2']
    metrics[Predictor == 'ESM1b', Predictor := 'ESM-1b']
  }
  
  saveRDS(metrics, file = paste('predictor/performance_', performance_output, '.rds', sep = ''))
  return(metrics)
}

###############################################################

## Change the code accordingly to train & test using the desired variant dataset, iterations, models.
iterations <- c('', 'foldxDdg - foldx_plddt', 'NoVariantClusters', 'NoSAAPdap', 'NoStructuralFeatures', 'varSite', 'NoSequenceFeatures', 'AlphaMissense_Score' , 'VariPred_Class', 'PHACT_Score', 'PHACTboost_Score', 'AlphaMissense_Score - VariPred_Class - PHACT_Score - PHACTboost_Score', 'CompPredictorsOnly')
# iterations <- c('')
models <- c("rf", "svmRadial", "svmLinear", "gbm", "glmnet", "glmboost", "knn", "nnet", "avNNet", "naive_bayes",
            "xgbTree", "xgbLinear", "ranger", "earth", "pcaNNet", "rpart", "C5.0", "bagFDA")
# models <- 'naive_bayes'
external_predictors <- c('CADD', 'ESM1b', 'AlphaMissense', 'VariPred', 'MutPred2', 'PolyPhen2', 'PHACT', 'PHACTboost')
preTuned_dt <- readRDS('resources/fineTuningParams_best.rds')

TrioNsight_features_dt <- fread('resources/TrioNsight_features_dataset.txt', sep = '\t')
TrioNsight_features_preprocessed_dt <- preprocess_training_data(TrioNsight_features_dt)
trainingSet_Overlaps <- readRDS('resources/trainingSet_Overlaps.rds')

all_models_all_methods_wDataOverlap <- train_all_methods_all_formulas(TrioNsight_features_preprocessed_dt, models = models, iterations = iterations, preTuned = preTuned_dt, model_output = 'all_models_all_methods_wDataOverlap')
predictions_all_models_all_methods_wDataOverlap <- optimizing_hyperparameters_retrieve_predictions(dt = TrioNsight_features_preprocessed_dt, model_input = all_models_all_methods_wDataOverlap, model_output = 'all_models_all_methods_wDataOverlap')

external_predictor_results <- fread('resources/external_predictions.txt', sep = '\t')
external_predictor_results_Bool <- external_predictor_results[, .(target_id, Class, AlphaMissense_Class, VariPred_Class, CADD_Class, ESM1b_Class, MutPred2_Class, PolyPhen2_Class, PHACT_Class, PHACTboost_Class)]
colnames(external_predictor_results_Bool) <- sapply(strsplit(colnames(external_predictor_results_Bool), '_', fixed = T), '[', 1)
colnames(external_predictor_results_Bool)[1] <- 'target_id'
external_predictor_results_Scores <- external_predictor_results[, .(target_id, Class, AlphaMissense_Score, CADD_Score, ESM1b_Score, MutPred2_Score, PolyPhen2_Score, PHACT_Score, PHACTboost_Score)]
colnames(external_predictor_results_Scores) <- sapply(strsplit(colnames(external_predictor_results_Scores), '_', fixed = T), '[', 1)
colnames(external_predictor_results_Scores)[1] <- 'target_id'

formulas <- names(predictions_all_models_all_methods_wDataOverlap[[1]])
models_all_methods_all_formulas_performances <- model_performance_calculation(predictions_models = predictions_all_models_all_methods_wDataOverlap, 
                                                                              model_names = models, 
                                                                              external_predictors = external_predictors, 
                                                                              external_predictor_results = external_predictor_results_Bool,
                                                                              data_overlap = NULL,
                                                                              performance_output = 'all_models_all_methods_wDataOverlap')

###############################################################

training_models_multipleruns <- function(dt, models, preTuned, variant_annotations, data_overlap, model_output){
  
  # LOOCV based on MSALocation
  unique_locs <- unique(dt$MSALocation)
  dt$target_id <- as.factor(dt$target_id)
  
  predictions_multiplerun_list <- list()
  
  formula <- as.formula("Class ~ . - target_id - MSALocation - Inherited_Seq")
  
  for(model in models){
    
    for(i in 1:20) {
      
      runNo <- paste('Training_', i, sep = '')
      cat(paste('At model: ', model, ', at run: ', i, sep = ''), '\n')
      
      preTuned_next <- preTuned_dt[[model]]
      preTuned_next <- expand.grid(preTuned_next)
      
      for (k in 1:length(unique_locs)) {
        # Current MSA location for testing
        current_loc <- unique_locs[k]
        
        # Training data: all except current MSA location
        train_data <- dt[dt$MSALocation != current_loc, ]
        
        # Testing data: only non-overlapping variants from current MSA location
        test_data <- dt[dt$MSALocation == current_loc & 
                          !dt$target_id %in% data_overlap$TrioNsight, ]
        
        # Skip if no non-overlapping variants to test
        if(nrow(test_data) == 0) next
        
        modelsWithLessParams <- c('knn', 'gbm', 'xgbTree', 'xgbLinear', 'glmboost', 'ranger', 'earth', 'rpart', 'bagFDA')
        
        if(model %in% modelsWithLessParams) {
          modelTrained <- train(formula,
                                data = train_data,
                                method = model,
                                preProc = c("center", "scale"),
                                trControl = trainControl(method = "none", classProbs = TRUE),
                                tuneGrid = preTuned_next,
                                metric = "ROC",  # Change this if not a classification problem
                                maximize = TRUE)  # TRUE if using accuracy or other metrics where higher is better
        }
        else {
          modelTrained <- train(formula,
                                data = train_data,
                                method = model,
                                preProc = c("center", "scale"),
                                trControl = trainControl(method = "none", classProbs = TRUE),
                                allowParallel = TRUE,
                                tuneGrid = preTuned_next,
                                metric = "ROC",  # Change this if not a classification problem
                                maximize = TRUE,  # TRUE if using accuracy or other metrics where higher is better
                                verbose = TRUE)
        }
        
        # Make predictions on test data
        predictions <- predict(modelTrained, newdata = test_data, type = "prob")
        predictions <- cbind(test_data$target_id, predictions)
        
        # Store results
        predictions_multiplerun_list[[model]][[runNo]][[current_loc]] <- predictions
        
      }
    }
  }
  
  predictions_multiplerun_parsed <- list()
  for(i in 1:length(predictions_multiplerun_list)){
    next_model <- predictions_multiplerun_list[[i]]
    model_name <- names(predictions_multiplerun_list)[i]
    for(j in 1:length(next_model)){
      next_model_iteration <- next_model[[j]]
      iteration_no <- names(next_model)[j]
      
      valid_indices <- which(!sapply(next_model_iteration, is.null))
      valid_predictions <- next_model_iteration[valid_indices]
      
      results <- do.call(rbind, valid_predictions)
      colnames(results)[1] <- 'target_id'
      
      final_results <- results %>%
        group_by(target_id) %>%
        summarise(
          Benign = mean(Benign),
          Pathogenic = mean(Pathogenic)
        ) %>%
        mutate(
          pred = ifelse(Benign > Pathogenic, "Benign", "Pathogenic")
        )
      
      final_results <- as.data.table(final_results)
      final_results <- merge(final_results, variant_annotations[,.(target_id, obs)], by = 'target_id')
      predictions_multiplerun_parsed[[model_name]][[iteration_no]] <- final_results
    }
  }
  
  saveRDS(predictions_multiplerun_parsed, file = paste('predictor/predictions_', model_output, '.rds', sep = ''))
  return(predictions_multiplerun_parsed)
}

###############################################################

TrioNsight_features_dt <- fread('resources/TrioNsight_features_dataset.txt', sep = '\t')
TrioNsight_features_preprocessed_dt <- preprocess_training_data(TrioNsight_features_dt)
trainingSet_Overlaps <- readRDS('resources/trainingSet_Overlaps.rds')

iterations <- ''
models <- c("rf", "svmRadial", "svmLinear", "gbm", "glmnet", "glmboost", "knn", "nnet", "avNNet", "naive_bayes",
            "xgbTree", "xgbLinear", "ranger", "earth", "pcaNNet", "rpart", "C5.0", "bagFDA")
external_predictors <- c('CADD', 'ESM1b', 'AlphaMissense', 'VariPred', 'MutPred2', 'PolyPhen2', 'PHACT', 'PHACTboost')

preTuned_dt <- readRDS('resources/fineTuningParams_best.rds')

external_predictor_results <- fread('resources/external_predictions.txt', sep = '\t')
variant_annotations <- external_predictor_results[, .(target_id, Class)]
colnames(variant_annotations)[2] <- 'obs'


predictions_all_models_multipleruns <- training_models_multipleruns(dt = TrioNsight_features_preprocessed_dt,
                                                                    models = models,
                                                                    preTuned = preTuned_dt,
                                                                    variant_annotations = variant_annotations,
                                                                    data_overlap = trainingSet_Overlaps, 
                                                                    model_output = 'all_models_multipleruns')
                                                                    
performance_models_all_methods_multipleruns <- model_performance_calculation(predictions_models = predictions_all_models_multipleruns, 
                                                                              model_names = models,
                                                                              external_predictors = NULL, 
                                                                              external_predictor_results = NULL,
                                                                              data_overlap = NULL,
                                                                              performance_output = 'all_models_multipleruns')

###############################################################

## Train a TrioNsight model by splitting training data into test and train sets
## This is a separate performance test from the previous one, which uses a leave-one-out cross-validation approach based on MSA locations.

TrioNsight_features_dt <- fread('resources/TrioNsight_features_dataset.txt', sep = '\t')
TrioNsight_features_dt$Class <- as.factor(TrioNsight_features_dt$Class)

# Split the variant dataset - ensure overlapping variants stay at the training set
split_indices <- list()

nonoverlapping_indices <- which(!TrioNsight_features_dt$target_id %in% trainingSet_Overlaps$TrioNsight)
test_dt <- TrioNsight_features_dt[nonoverlapping_indices,]
test_dt <- unique(test_dt[,.(target_id, Class)])
test_dt <- test_dt[, .SD[sample(.N, min(7, .N))], by = Class]

TrioNsight_features_dt[, Inherited_Seq := NULL]
TrioNsight_features_dt[, MSALocation := NULL]

# You can access the split datasets like this:
train_dt <- TrioNsight_features_dt[!target_id %in% test_dt$target_id, ]
test_dt <- TrioNsight_features_dt[target_id %in% test_dt$target_id, ]

train_dt_preprocessed <- preprocess_training_data(train_dt)
test_dt_preprocessed <- preprocess_training_data(test_dt)

train_target_ids <- train_dt_preprocessed$target_id
test_target_ids <- test_dt_preprocessed$target_id
train_dt_preprocessed[, target_id := NULL]
test_dt_preprocessed[, target_id := NULL]

indices <- split_indices$train  # Just one set of indices for training

fitControl <- trainControl(
  method = "none",  # Single train/test split
  index = indices,
  classProbs = TRUE,
  ## summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

param_grid <- expand.grid(
  laplace = 0,
  usekernel = FALSE,
  adjust = 1
)

# Modified train call without tuning parameters
random_split_model <- train(
  Class ~ . ,
  data = train_dt_preprocessed,
  trControl = fitControl,
  method = 'naive_bayes',
  preProc = c("center", "scale"),
  metric = 'ROC',
  tuneGrid = param_grid, 
  verbose = TRUE
)

predictions_random_split_model <- predict(random_split_model, newdata = test_dt_preprocessed, type = "prob")
predictions_random_split_model <- as.data.table(predictions_random_split_model)
predictions_random_split_model[, target_id := test_target_ids]

predictions_random_split_model <- aggregate(predictions_random_split_model, by = list(predictions_random_split_model$target_id), mean)
predictions_random_split_model <- as.data.table(predictions_random_split_model)
predictions_random_split_model[, target_id := NULL]
colnames(predictions_random_split_model)[1] <- 'target_id'

predictions_random_split_model[, Prediction := ifelse(Pathogenic > Benign, "Pathogenic", "Benign")]

observed <- unique(TrioNsight_features_dt[target_id %in% predictions_random_split_model$target_id][,.(target_id, Class)])
predictions_random_split_model <- merge(predictions_random_split_model, observed, by = 'target_id')
colnames(predictions_random_split_model)[4:5] <- c('pred', 'obs')
predictions_random_split_model <- list(list(predictions_random_split_model))
names(predictions_random_split_model) <- 'TrioNsight'
names(predictions_random_split_model[[1]]) <- 'random_split_model'

saveRDS(predictions_random_split_model, file = 'predictor/predictions_random_split_model.rds')

performance_random_split_model <- model_performance_calculation(predictions_models = predictions_random_split_model, 
                                                                                           model_names = 'TrioNsight', 
                                                                                           external_predictors = NULL, 
                                                                                           external_predictor_results = NULL,
                                                                                           data_overlap = NULL,
                                                                                           performance_output = 'random_split_model')

##################################################################

## Train TrioNsight on usual variant dataset and test on non-overlapping variant dataset

train_all_methods_all_formulas_woDataOverlap <- function(dt, models, iterations, preTuned, variant_annotations, data_overlap, model_output){
  
  formula_list <- initialize_formulas(dt = dt)
  # LOOCV based on MSALocation
  unique_locs <- unique(dt$MSALocation)
  
  # First, ensure target_id is handled properly - convert to character instead of factor
  dt$target_id <- factor(dt$target_id)
  
  # Initialize for storing predictions
  trained_models <- list()
  predictions_list <- list()
  
  ## dt[, rowIndex := NULL]
  
  for(model in models){
    
    for(i in 1:length(iterations)) {
      if(i == 1){
        formulaName <- 'Inclusive'
      }
      else if (i == 2){
        formulaName <- 'NofoldX'
      }
      else if (i == 3){
        formulaName <- 'NoVariantClusters'
      }
      else if (i == 4){
        formulaName <- 'NoSAAPdap'
      }
      else if (i == 5){
        formulaName <- 'NoStructuralFeatures'
      }
      else if (i == 6){
        formulaName <- 'NovarSite'
      }
      else if (i == 7){
        formulaName <- 'NoSequenceFeatures'
      }
      else if (i == 8){
        formulaName <- 'NoAlphaMissense'
      }
      else if (i == 9){
        formulaName <- 'NoVariPred'
      }
      else if (i == 10){
        formulaName <- 'NoPHACT'
      }
      else if (i == 11){
        formulaName <- 'NoPHACTboost'
      }
      else if (i == 12){
        formulaName <- 'Exclusionary'
      }
      else if (i == 13){
        formulaName <- 'CompPredictorsOnly'
      }
      
      cat(paste('Model: ', model, ', Iteration: ', formulaName, sep = ''), '\n')
      
      formula <- as.formula("Class ~ . - target_id - MSALocation - Inherited_Seq")
      
      if(formulaName != 'Inclusive'){
        if(formulaName == 'NoVariantClusters'){
          formula <- as.formula(paste(deparse(formula), "-", formula_list$variant_clusters_names))
        }
        else if(formulaName == 'NoSAAPdap'){
          formula <- as.formula(paste(deparse(formula), "-", formula_list$SAAPdap))
        }
        else if(formulaName == 'NoStructuralFeatures'){
          formula <- as.formula(paste(deparse(formula), "-", formula_list$StructuralFeatures))
        }
        else if(formulaName == 'NoSequenceFeatures'){
          formula <- as.formula(paste(deparse(formula), "-", formula_list$SequenceFeatures))
        }
        else if(formulaName == 'CompPredictorsOnly'){
          formula <- as.formula(paste(deparse(formula), "-", formula_list$CompPredictorsOnly))
        }
        else{
          formula <- as.formula(paste(deparse(formula), "-", iterations[i]))
        }
      }
      cat(paste('The formula is:', Reduce(paste, deparse(formula)), sep = ' '), '\n')
      
      for (k in 1:length(unique_locs)) {
        # Current MSA location for testing
        current_loc <- unique_locs[k]
        
        # Training data: all except current MSA location
        train_data <- dt[dt$MSALocation != current_loc, ]
        
        # Testing data: only non-overlapping variants from current MSA location
        test_data <- dt[dt$MSALocation == current_loc & 
                          !dt$target_id %in% data_overlap$TrioNsight, ]
        
        # Skip if no non-overlapping variants to test
        if(nrow(test_data) == 0) next
        
        modelsWithLessParams <- c('knn', 'gbm', 'xgbTree', 'xgbLinear', 'glmboost', 'ranger', 'earth', 'rpart', 'bagFDA')
        
        preTuned_next <- preTuned_dt[[model]]
        preTuned_next <- expand.grid(preTuned_next)
        
        if(model %in% modelsWithLessParams) {
          modelTrained <- train(formula,
                                data = train_data,
                                method = model,
                                preProc = c("center", "scale"),
                                trControl = trainControl(method = "none", classProbs = TRUE),
                                tuneGrid = preTuned_next,
                                metric = "ROC",  # Change this if not a classification problem
                                maximize = TRUE)  # TRUE if using accuracy or other metrics where higher is better
        }
        else {
          modelTrained <- train(formula,
                                data = train_data,
                                method = model,
                                preProc = c("center", "scale"),
                                trControl = trainControl(method = "none", classProbs = TRUE),
                                allowParallel = TRUE,
                                tuneGrid = preTuned_next,
                                metric = "ROC",  # Change this if not a classification problem
                                maximize = TRUE,  # TRUE if using accuracy or other metrics where higher is better
                                verbose = TRUE)
        }
        
        # Make predictions on test data
        predictions <- predict(modelTrained, newdata = test_data, type = "prob")
        predictions <- cbind(test_data$target_id, predictions)
        
        # Store results
        trained_models[[model]][[formulaName]][[current_loc]] <- modelTrained
        predictions_list[[model]][[formulaName]][[current_loc]] <- predictions
        
      }
      
    }
    
  }
  
  predictions_parsed <- list()
  models_parsed <- list()
  for(i in 1:length(predictions_list)){
    next_model <- predictions_list[[i]]
    model_name <- names(predictions_list)[i]
    for(j in 1:length(next_model)){
      next_model_formula <- next_model[[j]]
      formula_name <- names(next_model)[j]
      
      valid_indices <- which(!sapply(next_model_formula, is.null))
      valid_predictions <- next_model_formula[valid_indices]
      
      results <- do.call(rbind, valid_predictions)
      colnames(results)[1] <- 'target_id'
      
      final_results <- results %>%
        group_by(target_id) %>%
        summarise(
          Benign = mean(Benign),
          Pathogenic = mean(Pathogenic)
        ) %>%
        mutate(
          pred = ifelse(Benign > Pathogenic, "Benign", "Pathogenic")
        )
      
      final_results <- as.data.table(final_results)
      final_results <- merge(final_results, variant_annotations[,.(target_id, obs)], by = 'target_id')
      predictions_parsed[[model_name]][[formula_name]] <- final_results
      models_parsed[[model_name]][[formula_name]] <- trained_models[[model_name]][[formula_name]][valid_indices]
    }
  }
  
  saveRDS(predictions_parsed, file = paste('predictor/predictions_', model_output, '.rds', sep = ''))
  saveRDS(models_parsed, file = paste('predictor/models_', model_output, '.rds', sep = ''))
  return(predictions_parsed)
}

##########################################

TrioNsight_features_dt <- fread('resources/TrioNsight_features_dataset.txt', sep = '\t')
TrioNsight_features_preprocessed_dt <- preprocess_training_data(TrioNsight_features_dt)
trainingSet_Overlaps <- readRDS('resources/trainingSet_Overlaps.rds')
preTuned_dt <- readRDS('resources/fineTuningParams_best.rds')

external_predictor_results <- fread('resources/external_predictions.txt', sep = '\t')

external_predictor_results_Bool <- external_predictor_results[, .(target_id, Class, AlphaMissense_Class, VariPred_Class, CADD_Class, ESM1b_Class, MutPred2_Class, PolyPhen2_Class, PHACT_Class, PHACTboost_Class)]
colnames(external_predictor_results_Bool) <- sapply(strsplit(colnames(external_predictor_results_Bool), '_', fixed = T), '[', 1)
colnames(external_predictor_results_Bool)[1] <- 'target_id'

variant_annotations <- external_predictor_results[, .(target_id, Class)]
colnames(variant_annotations)[2] <- 'obs'

predictions_all_methods_all_formulas_woDataOverlap <- train_all_methods_all_formulas_woDataOverlap(dt = TrioNsight_features_preprocessed_dt, 
                                                      models = models, 
                                                      iterations = iterations,
                                                      preTuned = preTuned_dt,
                                                      variant_annotations = variant_annotations,
                                                      data_overlap = trainingSet_Overlaps,
                                                      model_output = 'TrioNsight')

performance_models_all_methods_all_formulas_woDataOverlap <- model_performance_calculation(predictions_models = predictions_all_methods_all_formulas_woDataOverlap, 
                                                                              model_names = models, 
                                                                              external_predictors = external_predictors, 
                                                                              external_predictor_results = external_predictor_results_Bool,
                                                                              data_overlap = trainingSet_Overlaps,
                                                                              performance_output = 'all_models_all_methods')
