
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

train_all_methods_all_formulas <- function(dt, models, iterations, model_output){
  
  # Initialize .folds_
  num_folds <- length(unique(dt$MSALocation))
  fold_cols <- paste(".folds_", seq_len(num_folds), sep = "")
  
  idxList <- list()
  # LOOCV based on MSALocation
  unique_locs <- unique(dt$MSALocation)
  
  for (i in 1:length(unique_locs)) {
    idxList[[fold_cols[i]]] <- which(dt$MSALocation != unique_locs[i])
  }
  fitControl <- trainControl(index = idxList, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE, verboseIter = TRUE)
  
  trained_models <- list()
  
  variant_clusters_names <- names(dt)[grepl('distance_to_cluster', names(dt))]
  variant_clusters_names <- paste(variant_clusters_names, collapse = " - ")
  
  SAAPdap <- c("Binding", "SProtFT_NAMES", "SProtFT_BOOL", "SProtFT_FEATURES", "SProtFTDist_BOOL", 
               "SProtFT_ACT_SITE", "SProtFT_BINDING", "SProtFT_CA_BIND", "SProtFT_DNA_BIND", "SProtFT_NP_BIND", 
               "SProtFT_METAL", "SProtFT_MOD_RES", "SProtFT_CARBOHYD", "SProtFT_MOTIF", "SProtFT_LIPID", 
               "Interface", "Relaccess", "Impact",	"HBonds",	"SPhobic", "CPhilic", "BCharge", "SSGeom", "Voids",	"MLargest1",
               "MLargest2", "MLargest3",	"MLargest4", "MLargest5", "MLargest6", "MLargest7", "MLargest8", "MLargest9", 
               "MLargest10", "NLargest1", "NLargest2", "NLargest3", "NLargest4", "NLargest5", "NLargest6", "NLargest7", 
               "NLargest8", "NLargest9", "NLargest10",	"Clash", "Glycine", "Proline", "CisPro")
  
  SAAPdap <- intersect(names(dt), SAAPdap)
  SAAPdap <- paste(SAAPdap, collapse = " - ")
  
  ## dt[, rowIndex := NULL]
  
  for(model in models){
    
    for(i in 1:length(iterations)) {
      
      if(i == 1){
        formulaName <- 'Inclusive'
      }
      else if(i == 2){
        formulaName <- 'NofoldX'
      }
      else if (i == 3){
        formulaName <- 'NovarSite'
      }
      else if (i == 4){
        formulaName <- 'NoAlphaMissense'
      }
      else if (i == 5){
        formulaName <- 'NoVariPred'
      }
      else if (i == 6){
        formulaName <- 'NoPHACTboost'
      }
      else if (i == 7){
        formulaName <- 'NoSAAPdap'
      }
      else if (i == 8){
        formulaName <- 'NoVariantClusters'
      }
      else if (i == 9){
        formulaName <- 'Exclusionary'
      }
      
      cat(paste('Model: ', model, ', Iteration: ', formulaName, sep = ''), '\n')
      
      formula <- as.formula("Class ~ . - target_id - MSALocation - Inherited_Seq")
      
      if(formulaName != 'Inclusive'){
        if(formulaName == 'NoVariantClusters'){
          formula <- as.formula(paste(deparse(formula), "-", variant_clusters_names))
        }
        else if(formulaName == 'NoSAAPdap'){
          formula <- as.formula(paste(deparse(formula), "-", SAAPdap))
        }
        else{
          formula <- as.formula(paste(deparse(formula), "-", iterations[i]))
        }
      }
      cat(paste('The formula is:', Reduce(paste, deparse(formula)), sep = ' '), '\n')
      
      modelsWithLessParams <- c('knn', 'gbm', 'xgbTree', 'xgbLinear', 'glmboost', 'ranger', 'earth', 'rpart', 'bagFDA')
      if(model %in% modelsWithLessParams){
        modelTrained <- train(formula, data = dt, 
                              method = model,
                              preProc = c("center", "scale"),
                              trControl = fitControl,
                              tuneLength = 5,  # Add tuning
                              metric = "ROC",  # Change this if not a classification problem
                              maximize = TRUE)  # TRUE if using accuracy or other metrics where higher is better
      }
      else{
        modelTrained <- train(formula, data = dt, 
                              method = model,
                              preProc = c("center", "scale"),
                              trControl = fitControl,
                              allowParallel = TRUE,
                              tuneLength = 5,  # Add tuning
                              metric = "ROC",  # Change this if not a classification problem
                              maximize = TRUE,  # TRUE if using accuracy or other metrics where higher is better
                              verbose = TRUE)
      }
      trained_models[[model]][[formulaName]] <- modelTrained
    }
  }
  saveRDS(trained_models, file = paste('predictor/', model_output, '.rds', sep = ''))
  return(trained_models)
}

optimizing_hyperparameters <- function(model_input, model_output){
  
  optimized_models <- list()
  
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
      optimized_models[[modelName]][[formulaName]] <- tmpPred
    }
  }
  
  dt[, rowIndex := .I]
  for(i in 1:length(optimized_models)){
    for(j in 1:length(optimized_models[[i]])){
      optimized_models[[i]][[j]] <- as.data.table(merge(optimized_models[[i]][[j]], dt[,.(target_id, rowIndex)], by = 'rowIndex'))
      optimized_models[[i]][[j]][, rowIndex := NULL]
      optimized_models[[i]][[j]][, pred := NULL]
      
      observed <- unique(optimized_models[[i]][[j]][, .(target_id, obs)])
      colnames(observed)[2] <- 'obs'
      optimized_models[[i]][[j]][, obs := NULL]
      
      if(length(optimized_models[[i]][[j]][is.nan(Benign)]$Benign) > 0){
        cat(paste('There are', length(optimized_models[[i]][[j]][is.nan(Benign)]$Benign), 'NAs in Benign\n', sep = ' '))
        optimized_models[[i]][[j]][is.nan(Benign)]$Benign <- 0.50001
      }
      if(length(optimized_models[[i]][[j]][is.nan(Pathogenic)]$Pathogenic) > 0){
        cat(paste('There are', length(optimized_models[[i]][[j]][is.nan(Pathogenic)]$Pathogenic), 'NAs in Pathogenic\n', sep = ' '))
        optimized_models[[i]][[j]][is.nan(Pathogenic)]$Pathogenic <- 0.49999
      }
      
      optimized_models[[i]][[j]] <- aggregate(optimized_models[[i]][[j]], by = list(optimized_models[[i]][[j]]$target_id), mean)
      optimized_models[[i]][[j]] <- as.data.table(optimized_models[[i]][[j]])
      optimized_models[[i]][[j]][, target_id := NULL]
      colnames(optimized_models[[i]][[j]])[1] <- 'target_id'
      
      optimized_models[[i]][[j]][, pred := ifelse(Pathogenic > Benign, "Pathogenic", "Benign")]
      optimized_models[[i]][[j]] <- merge(optimized_models[[i]][[j]], observed, by = 'target_id')
    }
  }
  saveRDS(optimized_models, file = paste('predictor/', model_output, '_optimized.rds', sep = ''))
  return(optimized_models)
}

model_performance_calculation <- function(optimized_models, model_names, external_predictors = NULL, external_predictor_results = NULL, performance_output){
  
  formulas <- names(optimized_models[[1]])
  
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
      
      if(formula %in% names(optimized_models[[model]])) {
        data <- optimized_models[[model]][[formula]]$obs
        data <- as.data.table(as.character(optimized_models[[model]][[formula]]$obs))
        colnames(data) <- 'Class'
        predicted <- optimized_models[[model]][[formula]]$pred
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
    data <- external_predictor_results
    
    for(predictor in external_predictors) {
      predicted <- data[[predictor]]
      
      metrics[Predictor == predictor, Precision := precision(predicted, data$Class, positive_class = "Pathogenic")]
      metrics[Predictor == predictor, Recall := recall(predicted, data$Class, positive_class = "Pathogenic")]
      metrics[Predictor == predictor, Accuracy := accuracy(predicted, data$Class)]
      metrics[Predictor == predictor, False_Positive_Rate := false_positive_rate(predicted, data$Class, positive_class = "Pathogenic")]
      metrics[Predictor == predictor, False_Negative_Rate := false_negative_rate(predicted, data$Class, positive_class = "Pathogenic")]
      metrics[Predictor == predictor, F1_Score := f1_score(predicted, data$Class, positive_class = "Pathogenic")]
      metrics[Predictor == predictor, MCC := mcc_score(predicted, data$Class, positive_class = "Pathogenic")]
    }
    
    metrics[Predictor == 'PolyPhen2', Predictor := 'PolyPhen-2']
    metrics[Predictor == 'ESM1b', Predictor := 'ESM-1b']
  }
  
  saveRDS(metrics, file = paste('predictor/', performance_output, '_performances.rds', sep = ''))
  return(metrics)
}

###############################################################

iterations <- c('', 'foldxDdg - foldx_plddt', 'varSite', 'AlphaMissense_Score' , 'VariPred_Class', 'PHACT_Score - PHACTboost_Score', 'SAAPdap', 'distance_to_clusters', 'AlphaMissense_Score - VariPred_Class - PHACT_Score - PHACTboost_Score')
models <- c("rf", "svmRadial", "svmLinear", "gbm", "glmnet", "glmboost", "knn", "nnet", "avNNet", "naive_bayes",
            "xgbTree", "xgbLinear", "ranger", "earth", "pcaNNet", "rpart", "C5.0", "bagFDA")
external_predictors <- c('CADD', 'ESM1b', 'AlphaMissense', 'VariPred', 'MutPred2', 'PolyPhen2', 'PHACTboost')


dt <- fread('resources/TrioNsight_features_dataset.txt', sep = '\t')
dt_preprocessed <- preprocess_training_data(dt)
models_all_methods_all_formulas <- train_all_methods_all_formulas(dt_preprocessed, models = models, iterations = iterations, model_output = 'models_all_methods_all_formulas')
models_all_methods_all_formulas_optimized <- optimizing_hyperparameters(model_input = models_all_methods_all_formulas, model_output = 'models_all_methods_all_formulas')

external_predictor_results <- fread('resources/external_predictions.txt', sep = '\t')
external_predictor_results_Bool <- external_predictor_results[, .(Class, AlphaMissense_Class, VariPred_Class, CADD_Class, ESM1b_Class, MutPred2_Class, PolyPhen2_Class, PHACTboost_Class)]
colnames(external_predictor_results_Bool) <- sapply(strsplit(colnames(external_predictor_results_Bool), '_', fixed = T), '[', 1)
external_predictor_results_Scores <- external_predictor_results[, .(Class, AlphaMissense_Score, CADD_Score, ESM1b_Score, MutPred2_Score, PolyPhen2_Score, PHACTboost_Score)]
colnames(external_predictor_results_Scores) <- sapply(strsplit(colnames(external_predictor_results_Scores), '_', fixed = T), '[', 1)

formulas <- names(models_all_methods_all_formulas_optimized[[1]])
models_all_methods_all_formulas_performances <- model_performance_calculation(optimized_models = models_all_methods_all_formulas_optimized, 
                                                                              model_names = models, 
                                                                              external_predictors = external_predictors, 
                                                                              external_predictor_results = external_predictor_results_Bool,
                                                                              performance_output = 'models_all_methods_all_formulas')

###############################################################

training_models_multipleruns <- function(dt, models, model_output){
  
  # Initialize .folds_
  num_folds <- length(unique(dt$MSALocation))
  fold_cols <- paste(".folds_", seq_len(num_folds), sep = "")
  
  idxList <- list()
  # LOOCV based on MSALocation
  unique_locs <- unique(dt$MSALocation)
  
  for (i in 1:length(unique_locs)) {
    idxList[[fold_cols[i]]] <- which(dt$MSALocation != unique_locs[i])
  }
  fitControl <- trainControl(index = idxList, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE, verboseIter = TRUE)
  
  trained_models_multiplerun <- list()
  
  variant_clusters_names <- names(dt)[grepl('distance_to_cluster', names(dt))]
  variant_clusters_names <- paste(variant_clusters_names, collapse = " - ")
  
  SAAPdap <- c("Binding", "SProtFT_NAMES", "SProtFT_BOOL", "SProtFT_FEATURES", "SProtFTDist_BOOL", 
               "SProtFT_ACT_SITE", "SProtFT_BINDING", "SProtFT_CA_BIND", "SProtFT_DNA_BIND", "SProtFT_NP_BIND", 
               "SProtFT_METAL", "SProtFT_MOD_RES", "SProtFT_CARBOHYD", "SProtFT_MOTIF", "SProtFT_LIPID", 
               "Interface", "Relaccess", "Impact",	"HBonds",	"SPhobic", "CPhilic", "BCharge", "SSGeom", "Voids",	"MLargest1",
               "MLargest2", "MLargest3",	"MLargest4", "MLargest5", "MLargest6", "MLargest7", "MLargest8", "MLargest9", 
               "MLargest10", "NLargest1", "NLargest2", "NLargest3", "NLargest4", "NLargest5", "NLargest6", "NLargest7", 
               "NLargest8", "NLargest9", "NLargest10",	"Clash", "Glycine", "Proline", "CisPro")
  
  SAAPdap <- intersect(names(dt), SAAPdap)
  SAAPdap <- paste(SAAPdap, collapse = " - ")
  
  formula <- as.formula("Class ~ . - target_id - MSALocation - Inherited_Seq")
  
  for(model in models){
    
    for(i in 1:20) {
      
      runNo <- paste('Training_', i, sep = '')
      cat(paste('At model: ', model, ', at run: ', i, sep = ''), '\n')
      
      modelsWithLessParams <- c('knn', 'gbm', 'xgbTree', 'xgbLinear', 'glmboost', 'ranger', 'earth', 'rpart', 'bagFDA')
      if(model %in% modelsWithLessParams){
        modelTrained <- train(formula, data = dt, 
                              method = model,
                              preProc = c("center", "scale"),
                              trControl = fitControl,
                              tuneLength = 5,  # Add tuning
                              metric = "ROC",  # Change this if not a classification problem
                              maximize = TRUE)  # TRUE if using accuracy or other metrics where higher is better
      }
      else{
        modelTrained <- train(formula, data = dt, 
                              method = model,
                              preProc = c("center", "scale"),
                              trControl = fitControl,
                              allowParallel = TRUE,
                              tuneLength = 5,  # Add tuning
                              metric = "ROC",  # Change this if not a classification problem
                              maximize = TRUE,  # TRUE if using accuracy or other metrics where higher is better
                              verbose = TRUE)
      }
      trained_models_multiplerun[[model]][[runNo]] <- modelTrained
    }
    saveRDS(trained_models_multiplerun, file = paste('predictor/', model_output, '.rds', sep = ''))
  }
  return(trained_models_multiplerun)
}

###############################################################

iterations <- c('', 'foldxDdg - foldx_plddt', 'varSite', 'AlphaMissense_Score' , 'VariPred_Class', 'PHACT_Score - PHACTboost_Score', 'SAAPdap', 'distance_to_clusters', 'AlphaMissense_Score - VariPred_Class - PHACT_Score - PHACTboost_Score')
models <- c("rf", "svmRadial", "svmLinear", "gbm", "glmnet", "glmboost", "knn", "nnet", "avNNet", "naive_bayes",
            "xgbTree", "xgbLinear", "ranger", "earth", "pcaNNet", "rpart", "C5.0", "bagFDA")
external_predictors <- c('CADD', 'ESM1b', 'AlphaMissense', 'VariPred', 'MutPred2', 'PolyPhen2', 'PHACTboost')

models_all_methods_multipleruns <- training_models_multipleruns(dt_preprocessed, models = models, model_output = 'models_all_methods_multipleruns')
models_all_methods_multipleruns_optimized <- optimizing_hyperparameters(model_input = models_all_methods_multipleruns, model_output = 'models_all_methods_multipleruns')
models_all_methods_multipleruns_performances <- model_performance_calculation(optimized_models = models_all_methods_multipleruns_optimized, 
                                                                              model_names = models,
                                                                              external_predictors = NULL, 
                                                                              external_predictor_results = NULL,
                                                                              performance_output = 'models_all_methods_multipleruns')

###############################################################

split_variants_by_msa <- function(df, test_size = 0.2, seed = 42) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # Calculate target counts for test set
  total_pathogenic <- sum(df$Class == "Pathogenic")
  total_benign <- sum(df$Class == "Benign")
  target_test_pathogenic <- floor(total_pathogenic * test_size)
  target_test_benign <- floor(total_benign * test_size)
  
  # Get MSA location info
  msa_info <- df %>%
    group_by(MSALocation) %>%
    summarise(
      pathogenic_count = sum(Class == "Pathogenic"),
      benign_count = sum(Class == "Benign"),
      total_count = n(),
      .groups = "drop"
    ) %>%
    # Add random numbers for shuffling
    mutate(random = runif(n()))
  
  # Initialize vectors for tracking
  test_indices <- c()
  train_indices <- c()
  current_test_pathogenic <- 0
  current_test_benign <- 0
  
  # Sort MSA locations randomly
  msa_info <- arrange(msa_info, random)
  
  # Iterate through MSA locations and assign to train/test
  for(i in 1:nrow(msa_info)) {
    msa_loc <- msa_info$MSALocation[i]
    path_count <- msa_info$pathogenic_count[i]
    benign_count <- msa_info$benign_count[i]
    
    # Check if we can add to test while maintaining balance
    can_add_to_test <- (
      current_test_pathogenic < target_test_pathogenic &&
        current_test_benign < target_test_benign &&
        (current_test_pathogenic + path_count <= target_test_pathogenic) &&
        (current_test_benign + benign_count <= target_test_benign)
    )
    
    # Get indices for current MSA location
    current_indices <- which(df$MSALocation == msa_loc)
    
    if(can_add_to_test) {
      test_indices <- c(test_indices, current_indices)
      current_test_pathogenic <- current_test_pathogenic + path_count
      current_test_benign <- current_test_benign + benign_count
    } else {
      train_indices <- c(train_indices, current_indices)
    }
  }
  
  return(list(train = train_indices, test = test_indices))
}

# Function to print split statistics
print_split_stats <- function(df, split_indices) {
  train_df <- df[split_indices$train, ]
  test_df <- df[split_indices$test, ]
  
  cat("\nSplit Statistics:\n")
  cat("Training set size:", nrow(train_df), "variants\n")
  cat("Test set size:", nrow(test_df), "variants\n")
  
  cat("\nTraining set:\n")
  print(table(train_df$Class))
  
  cat("\nTest set:\n")
  print(table(test_df$Class))
  
  cat("\nUnique MSA locations:\n")
  cat("Training:", length(unique(train_df$MSALocation)), "\n")
  cat("Test:", length(unique(test_df$MSALocation)), "\n")
}

dt <- fread('resources/TrioNsight_features_dataset.txt', sep = '\t')
dt$Class <- as.factor(dt$Class)
# Perform the split
split_indices <- split_variants_by_msa(dt, test_size = 0.2)

# Print statistics
print_split_stats(dt, split_indices)

dt[, Inherited_Seq := NULL]
dt[, MSALocation := NULL]

# You can access the split datasets like this:
train_dt <- dt[split_indices$train, ]
test_dt <- dt[split_indices$test, ]

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

saveRDS(random_split_model, file = 'predictor/random_split_model_optimized.rds')
split_model_performance <- predict(random_split_model, newdata = test_dt_preprocessed, type = "prob")
split_model_performance <- as.data.table(split_model_performance)
split_model_performance[, target_id := test_target_ids]

split_model_performance <- aggregate(split_model_performance, by = list(split_model_performance$target_id), mean)
split_model_performance <- as.data.table(split_model_performance)
split_model_performance[, target_id := NULL]
colnames(split_model_performance)[1] <- 'target_id'

split_model_performance[, Prediction := ifelse(Pathogenic > Benign, "Pathogenic", "Benign")]

observed <- unique(dt[target_id %in% split_model_performance$target_id][,.(target_id, Class)])
random_split_model_performance <- merge(split_model_performance, observed, by = 'target_id')
saveRDS(random_split_model_performance, file = 'predictor/random_split_model_performance.rds')

###############################################################

dt <- fread('resources/TrioNsight_features_dataset.txt', sep = '\t')
dt_preprocessed <- preprocess_training_data(dt)

model_features <- setdiff(names(dt_preprocessed), c('target_id', 'MSALocation', 'Inherited_Seq'))
dt_clean <- dt_preprocessed[, ..model_features]
modelTrained <- train(Class ~ ., data = dt_clean,
                      method = 'naive_bayes',
                      preProc = c("center", "scale"),
                      trControl = fitControl,
                      allowParallel = TRUE,
                      tuneLength = 5,  # Add tuning
                      metric = "ROC",  # Change this if not a classification problem
                      maximize = TRUE,  # TRUE if using accuracy or other metrics where higher is better
                      verbose = TRUE)

TrioNsight_impactMap_features <- readRDS('resources//TrioNsight_impactMap_features.rds')
TrioNsight_impactMap_features_cleaned <- preprocess_training_data(TrioNsight_impactMap_features)

TrioNsight_impactMap_predicted <- predict(modelTrained, newdata = TrioNsight_impactMap_features, type = "prob")
TrioNsight_impactMap_predicted <- as.data.table(TrioNsight_impactMap_predicted)
TrioNsight_impactMap_predicted[, target_id := TrioNsight_impactMap_features$target_id]

# Calculate mean probabilities for each variant
TrioNsight_impactMap_predicted <- TrioNsight_impactMap_predicted[, .(
  mean_benign = mean(Benign),
  mean_pathogenic = mean(Pathogenic),
  sd_pathogenic = sd(Pathogenic)  # Standard deviation to check consistency
), by = target_id]

# Sort by pathogenic probability
setorder(TrioNsight_impactMap_predicted, -mean_pathogenic)

# Calculate some potential cutoffs
cutoff_stats <- TrioNsight_impactMap_predicted[, .(
  q90 = quantile(mean_pathogenic, 0.90),
  q95 = quantile(mean_pathogenic, 0.95),
  q99 = quantile(mean_pathogenic, 0.99)
)]

# Look at variants with high disagreement between predictions
inconsistent_preds <- TrioNsight_impactMap_predicted[sd_pathogenic > 0.2]

# Create confidence levels
TrioNsight_impactMap_predicted[, Prediction := case_when(
  sd_pathogenic < 0.05 & mean_pathogenic > 0.75 ~ "Pathogenic",
  sd_pathogenic >= 0.05 & mean_pathogenic > 0.75 ~ "Likely pathogenic",
  sd_pathogenic < 0.05 & mean_pathogenic < 0.25 ~ "Benign",
  sd_pathogenic >= 0.05 & mean_pathogenic < 0.25 ~ "Likely benign",
  sd_pathogenic > 0.2 ~ "Uncertain",
  TRUE ~ "Uncertain"
)]

TrioNsight_impactMap_predicted[, Accession := sapply(strsplit(target_id, '_', fixed = T), '[', 1)]
TrioNsight_impactMap_predicted[, Residue := sapply(strsplit(target_id, '_', fixed = T), '[', 2)]
TrioNsight_impactMap_predicted[, Residue := extractNumbersFromString(Residue)]
TrioNsight_impactMap_predicted[, Mutant := sapply(strsplit(target_id, '_', fixed = T), '[', 2)]
TrioNsight_impactMap_predicted[, Mutant := substrRight(TrioNsight_impactMap_predicted$Mutant, 1)]

aa_codes_ordered <- as.data.table(names(aa_codes))
aa_codes_ordered[, Order := 1:20]
aa_codes_ordered[, Order := as.numeric(Order)]
colnames(aa_codes_ordered)[1] <- 'Mutant'

TrioNsight_impactMap_predicted <- merge(TrioNsight_impactMap_predicted, aa_codes_ordered, by = 'Mutant')
TrioNsight_impactMap_predicted <- TrioNsight_impactMap_predicted[with(TrioNsight_impactMap_predicted, order(Accession, Residue, Order, decreasing = F))]

TrioNsight_impactMap_predicted[, Order := NULL]
TrioNsight_impactMap_predicted[, Residue := NULL]
TrioNsight_impactMap_predicted[, Accession := NULL]
TrioNsight_impactMap_predicted[, Mutant := NULL]

fwrite(TrioNsight_impactMap_predicted, file = 'TrioNsight_all_predictions.txt', sep = '\t')