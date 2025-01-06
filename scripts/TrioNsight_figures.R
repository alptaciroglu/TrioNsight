# Load required libraries
if (!requireNamespace("data.table", quietly = TRUE)) {
  BiocManager::install('data.table')
}
library(data.table)
if (!requireNamespace("flextable", quietly = TRUE)) {
  BiocManager::install('flextable')
}
library(flextable)
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  BiocManager::install('ggplot2')
}
library(ggplot2)
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  BiocManager::install('openxlsx')
}
library(openxlsx)

##############################################
################# Figure 1 ###################
##############################################

## Part A: Please refer to https://github.com/alptaciroglu/drawProteins for the code to generate the image
## Part B: Figure created using the PyMOL software

##############################################
################# Figure 2 ###################
##############################################

## Figure created using the draw.io web application 

##############################################
################# Figure 3 ###################
##############################################

## Figure created using the draw.io web application 

##############################################
################# Figure 4 ###################
##############################################

## Figure created using the PyMOL software

##############################################
################# Figure 5 ###################
##############################################

## Figure created using the draw.io web application 

##############################################
################# Table 1 ####################
##############################################

# Create the data.table
models <- data.table(Methods = c("rf", "svmRadial", "svmLinear", "gbm", "glmnet", "glmboost", "knn", "nnet", "avNNet", "naive_bayes", "xgbTree", "xgbLinear", "ranger", "earth", "pcaNNet", "rpart", "C5.0", "bagFDA"))

# Add numbering to the Methods
models[, Numbered_Methods := paste0(seq_len(.N), ". ", Methods)]

# Create a matrix with 4 columns and fill it with the numbered methods
matrix_data <- matrix(c(models$Numbered_Methods, rep("", 2)), ncol = 4, byrow = TRUE)

# Create a flextable
flextable_object <- flextable(as.data.frame(matrix_data)) %>% 
  autofit() %>% 
  theme_booktabs() %>% 
  bold(part = "header")

# Remove column names
flextable_object <- delete_part(flextable_object, part = "header")

# Add a single centralized "Methods" header
flextable_object <- add_header_row(flextable_object, values = "Table 1: Evaluated Machine Learning Methods for TrioNsight", colwidths = 4)

# Add a line just after the header
flextable_object <- hline(flextable_object, i = 1, border = fp_border(color = "black", width = 2), part = "header")

# Set equal width for all columns
flextable_object <- width(flextable_object, width = 2)

# Set font family and size for the entire table
flextable_object <- font(flextable_object, fontname = "Cambria", part = "all")
flextable_object <- fontsize(flextable_object, size = 11, part = "all")
flextable_object <- fontsize(flextable_object, size = 13, part = "header")

# Center align the "Methods" header
flextable_object <- align(flextable_object, align = "left", part = "header")

# Make the "Methods" header bold
flextable_object <- bold(flextable_object, bold = TRUE, part = "header")

save_as_image(flextable_object, path = "Table 1.png", width = 12, height = 8, units = "in", res = 1000)
flextable::save_as_html(flextable_object, path = "Table 1.html")

##############################################
################# Figure 6 ###################
##############################################

## Figure created using the draw.io web application 

##############################################
################# Figure 7 ###################
##############################################

## Figure created using the PyMOL software

##############################################
################# Figure 8 ###################
##############################################

## Figure created using the PyMOL software

##############################################
################# Figure 9 ###################
##############################################

## Figure created using the PyMOL software

##############################################
################# Table 2.1 ##################
##############################################

best_n_analysis <- readRDS("C:/Users/AlperenTaciroglu/Desktop/Work/Evolution/Results/2024.09 - September/RDataFiles/best_n_analysis.rds")

best_clustering <- rbind(best_n_analysis$best_n_results$TrioN_DHPH_1290_1594_A_1nty$p_values, 
                         best_n_analysis$best_n_results$TrioN_DHPH_1290_1594_A_1nty$x_squared)
best_clustering <- as.data.frame(best_clustering)
rownames(best_clustering) <- c("p_values", "x_squared")
best_clustering <- best_clustering[,1:7]
best_clustering <- t(best_clustering)
best_clustering <- as.data.frame(best_clustering)
rownames(best_clustering) <- 2:8

best_clustering$p_values <- round(best_clustering$p_values, 4)
best_clustering$x_squared <- round(best_clustering$x_squared, 4)
colnames(best_clustering) <- c('p-values', 'x-squared')

flextable_object <- flextable(best_clustering %>% tibble::rownames_to_column("Number of Clusters")) %>% 
  autofit() %>% 
  theme_booktabs() %>% 
  bold(part = "header") %>%
  align(align = "center", part = "all")

flextable_object <- add_header_row(flextable_object, values = "Table 2.1: Best clustering by chi-square test", colwidths = 3)
flextable_object <- hline_top(flextable_object, border = fp_border_default(width = 0), part = "header")

# Left align the "Methods" header
flextable_object <- align(flextable_object, i = 1, j = NULL, align = "left", part = "header")
# Make the "Methods" header bold
flextable_object <- bold(flextable_object, bold = TRUE, part = "header")

flextable_object <- font(flextable_object, fontname = "Cambria", part = "all")
flextable_object <- fontsize(flextable_object, size = 11, part = "all")
flextable_object <- fontsize(flextable_object, size = 12, part = "header")
flextable_object <- fontsize(flextable_object, size = 10, part = "footer")
flextable_object <- line_spacing(flextable_object, i = NULL, j = NULL, space = 0.7, part = "body")

toHighlight <- which.min(as.numeric(unlist(best_clustering[,1])))
flextable_object %<>% flextable::highlight(toHighlight, 2, color = 'wheat')

save_as_image(flextable_object, path = "Table 2.1.png", width = 12, height = 6, units = "in", res = 1000)
flextable::save_as_html(flextable_object, path = "Table 2.1.html")

##############################################
################# Table 2.2 ##################
##############################################

TrioN_DHPH_1290_1594_A_1nty_strucclus_n5 <- fread('StartUpFiles/Current/StructureAnalysis/TrioN_like_allAlg_TrioN_DH_strucclus_res_Unambigious/TrioN_DHPH_1290_1594_A_1nty_n5.txt', sep = '\t', skip = 62, header = F)

Clusters <- TrioN_DHPH_1290_1594_A_1nty_strucclus_n5$V1[grepl('Cluster', TrioN_DHPH_1290_1594_A_1nty_strucclus_n5$V1, fixed = T)]
ClusterNames <- sapply(strsplit(Clusters, ' (', fixed = T), '[', 1)
ClusterCounts <- sapply(strsplit(Clusters, ' (', fixed = T), '[', 2)
ClusterCounts <- as.numeric(gsub(' members)', '', ClusterCounts, fixed = T))

xx <- sapply(strsplit(sapply(strsplit(TrioN_DHPH_1290_1594_A_1nty_strucclus_n5$V1, ' # A ', fixed = T), '[', 2), ' ', fixed = T), '[', 1)
xx <- as.numeric(xx[!is.na(xx)])
yy <- sapply(strsplit(sapply(strsplit(TrioN_DHPH_1290_1594_A_1nty_strucclus_n5$V1, ' # A ', fixed = T), '[', 2), ' ', fixed = T), '[', 4)
yy <- yy[!is.na(yy)]
zz <- c(rep(ClusterNames[1], ClusterCounts[1]), rep(ClusterNames[2], ClusterCounts[2]), 
        rep(ClusterNames[3], ClusterCounts[3]), rep(ClusterNames[4], ClusterCounts[4]), rep(ClusterNames[5], ClusterCounts[5]))

TrioN_DHPH_1290_1594_A_1nty_strucclus_n5 <- as.data.table(cbind(xx, yy, zz))
colnames(TrioN_DHPH_1290_1594_A_1nty_strucclus_n5) <- c('PDB_Loc', 'Class', 'Cluster')
TrioN_DHPH_1290_1594_A_1nty_strucclus_n5 <- TrioN_DHPH_1290_1594_A_1nty_strucclus_n5[!is.na(TrioN_DHPH_1290_1594_A_1nty_strucclus_n5$PDB_Loc)]
TrioN_DHPH_1290_1594_A_1nty_strucclus_n5[, PDB_Loc := NULL]

clusteringTable <- table(TrioN_DHPH_1290_1594_A_1nty_strucclus_n5$Cluster, TrioN_DHPH_1290_1594_A_1nty_strucclus_n5$Class)
clusteringTable <- table2dataframe(clusteringTable)
colnames(clusteringTable)[1] <- 'Clusters'
clusteringTable$ClusterSize <- clusteringTable$Pathogenic + clusteringTable$Benign

flextable_object <- flextable(clusteringTable) %>% 
  autofit() %>% 
  theme_booktabs() %>% 
  bold(part = "header") %>%
  align(align = "center", part = "all")

flextable_object <- add_header_row(flextable_object, values = "Table 2.2: Variant classes for clustering n=5", colwidths = 4)
flextable_object <- hline_top(flextable_object, border = fp_border(color = "black", width = 0), part = "header")

# Left align the "Methods" header
flextable_object <- align(flextable_object, i = 1, j = NULL, align = "left", part = "header")
# Make the "Methods" header bold
flextable_object <- bold(flextable_object, bold = TRUE, part = "header")

flextable_object <- font(flextable_object, fontname = "Cambria", part = "all")
flextable_object <- fontsize(flextable_object, size = 12, part = "header")
flextable_object <- fontsize(flextable_object, size = 11, part = "body")
flextable_object <- fontsize(flextable_object, size = 10, part = "footer")

save_as_image(flextable_object, path = "Table 2.2.png", width = 12, height = 6, units = "in", res = 1000)
flextable::save_as_html(flextable_object, path = "Table 2.2.html")

##############################################
################# Table 3 ####################
##############################################

model_performances <- readRDS('predictor/models_all_methods_all_formulas_performances.rds')

metrics <- model_performances
selected_methods <- c('CADD', 'ESM-1b', 'AlphaMissense', 'VariPred', 'MutPred2', 'PolyPhen-2', 'PHACTboost', 'naive_bayes~Inclusive')
metrics <- metrics[Predictor %in% selected_methods,]
metrics$Predictor[8] <- 'TrioNsight'

# Create a visually appealing table with flextable
flextable_object <- flextable(metrics) %>% 
  autofit() %>% 
  theme_booktabs() %>% 
  bold(part = "header") %>%
  align(align = "center", part = "all") 

flextable_object <- add_header_row(flextable_object, values = "Table 3: Performance evaluation of TrioNsight with its closest competitors", colwidths = 8)
flextable_object <- hline_top(flextable_object, border = fp_border_default(width = 0), part = "header")

# Left align the "Methods" header
flextable_object <- align(flextable_object, i = 1, j = NULL, align = "left", part = "header")
# Make the "Methods" header bold
flextable_object <- bold(flextable_object, bold = TRUE, part = "header")
flextable_object <- font(flextable_object, fontname = "Cambria", part = "all")
flextable_object <- fontsize(flextable_object, size = 11, part = "all")
flextable_object <- fontsize(flextable_object, size = 12, part = "header")
flextable_object <- fontsize(flextable_object, size = 10, part = "footer")

metrics <- as.data.frame(metrics)

for(i in 2:ncol(metrics)){
  
  if(colnames(metrics)[i] == 'Precision'){
    toHighligth <- which.max(as.numeric(unlist(metrics[,i])))
  }
  else if(colnames(metrics)[i] == 'Recall'){
    toHighligth <- which.max(as.numeric(unlist(metrics[,i])))
  }
  else if(colnames(metrics)[i] == 'Accuracy'){
    toHighligth <- which.max(as.numeric(unlist(metrics[,i])))
  }
  else if(colnames(metrics)[i] == 'False_Positive_Rate'){
    toHighligth <- which.min(as.numeric(unlist(metrics[,i])))
  }
  else if(colnames(metrics)[i] == 'False_Negative_Rate'){
    toHighligth <- which.min(as.numeric(unlist(metrics[,i])))
  }
  else if(colnames(metrics)[i] == 'F1_Score'){
    toHighligth <- which.max(as.numeric(unlist(metrics[,i])))
  }
  else if(colnames(metrics)[i] == 'MCC'){
    toHighligth <- which.max(as.numeric(unlist(metrics[,i])))
  }
  flextable_object %<>% flextable::highlight(toHighligth, i, color = 'wheat')
}

save_as_image(flextable_object, path = "Table 3.png", width = 16, height = 6, units = "in", res = 1000)
save_as_html(flextable_object, path = "Table 3.html")

##############################################
################# Table 4 ####################
##############################################

TrioNsight_distinguishedPredictions <- readRDS('resources/TrioNsight_distinguishedPredictions.rds')

flextable_object <- flextable(TrioNsight_distinguishedPredictions) %>% 
  autofit() %>% 
  theme_booktabs() %>% 
  bold(part = "header") %>%
  align(align = "center", part = "all")

flextable_object <- set_caption(flextable_object, as_paragraph(
  as_chunk("Table 4: TrioNsight correctly classifies variants that are misclassified by other predictors", props = fp_text_default(font.family = "Cambria", font.size = 12, bold = TRUE))
), word_stylename = "Table Caption", align_with_table = FALSE)

flextable_object <- style(flextable_object, j = "Variation ID", part = "body", pr_t = fp_text_default(underlined = TRUE))

flextable_object <- align(flextable_object, i = 1, j = NULL, align = "center", part = "header")
flextable_object <- bold(flextable_object, bold = TRUE, part = "header")
flextable_object <- vline(flextable_object, i = NULL, j = "Classification", border = NULL, part = "all")

flextable_object <- compose(
  x = flextable_object, j = "Variation ID",
  value = as_paragraph(
    hyperlink_text(x = `Variation ID`, url = url
    )
  )
)
flextable_object <- delete_columns(flextable_object, j = "url")

flextable_object <- font(flextable_object, fontname = "Cambria", part = "all")
flextable_object <- fontsize(flextable_object, size = 12, part = "header")
flextable_object <- fontsize(flextable_object, size = 11, part = "body")
flextable_object <- fontsize(flextable_object, size = 10, part = "footer")

toHighlight <- which(TrioNsight_distinguishedPredictions$Classification != TrioNsight_distinguishedPredictions$AlphaMissense)
flextable_object %<>% flextable::highlight(toHighlight, i = toHighlight, j = "AlphaMissense", color = 'wheat', part = "body")
toHighlight <- which(TrioNsight_distinguishedPredictions$Classification != TrioNsight_distinguishedPredictions$PHACTboost)
flextable_object %<>% flextable::highlight(toHighlight, i = toHighlight, j = "PHACTboost", color = 'wheat', part = "body")
toHighlight <- which(TrioNsight_distinguishedPredictions$Classification != TrioNsight_distinguishedPredictions$VariPred)
flextable_object %<>% flextable::highlight(toHighlight, i = toHighlight, j = "VariPred", color = 'wheat', part = "body")
toHighlight <- which(TrioNsight_distinguishedPredictions$Classification != TrioNsight_distinguishedPredictions$MutPred2)
flextable_object %<>% flextable::highlight(toHighlight, i = toHighlight, j = "MutPred2", color = 'wheat', part = "body")
toHighlight <- which(TrioNsight_distinguishedPredictions$Classification != TrioNsight_distinguishedPredictions$CADD)
flextable_object %<>% flextable::highlight(toHighlight, i = toHighlight, j = "CADD", color = 'wheat', part = "body")
toHighlight <- which(TrioNsight_distinguishedPredictions$Classification != TrioNsight_distinguishedPredictions$`ESM-1b`)
flextable_object %<>% flextable::highlight(toHighlight, i = toHighlight, j = "ESM-1b", color = 'wheat', part = "body")
toHighlight <- which(TrioNsight_distinguishedPredictions$Classification != TrioNsight_distinguishedPredictions$`PolyPhen-2`)
flextable_object %<>% flextable::highlight(toHighlight, i = toHighlight, j = "PolyPhen-2", color = 'wheat', part = "body")
toHighlight <- which(TrioNsight_distinguishedPredictions$Classification != TrioNsight_distinguishedPredictions$TrioNsight)
flextable_object %<>% flextable::highlight(toHighlight, i = toHighlight, j = "TrioNsight", color = 'wheat', part = "body")

save_as_image(flextable_object, path = "Table 4.png", width = 16, height = 8, units = "in", res = 1000)
flextable::save_as_html(flextable_object, path = "Table 4.html")

##############################################
################# Table 5 ####################
##############################################

external_predictor_results <- fread('resources/external_predictions.txt', sep = '\t')
external_predictor_results_Bool <- external_predictor_results[, .(target_id, Class, AlphaMissense_Class, VariPred_Class, CADD_Class, ESM1b_Class, MutPred2_Class, PolyPhen2_Class, PHACTboost_Class)]
colnames(external_predictor_results_Bool) <- sapply(strsplit(colnames(external_predictor_results_Bool), '_', fixed = T), '[', 1)

trainingSet_Overlaps <- readRDS('resources/trainingSet_Overlaps.rds')

external_predictor_results_Bool[Class == 'Ambigious', Class := 'Pathogenic']
predictors <- c('CADD', 'ESM1b', 'AlphaMissense', 'VariPred', 'MutPred2', 'PolyPhen2', 'PHACTboost')
full_data <- external_predictor_results_Bool
metrics <- data.table(Predictor = predictors)

for(predictor in predictors) {
  data <- full_data[!target %in% trainingSet_Overlaps[[predictor]]]
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
colnames(metrics) <- c('Predictor', 'Precision', 'Recall', 'Accuracy', 'False Positive Rate', 'False Negative Rate', 'F1 Score', 'MCC')

# Create a visually appealing table with flextable
flextable_object <- flextable(metrics) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  bold(part = "header")

flextable_object <- add_header_row(flextable_object, values = "Table 5: Benchmarking Methods Performance (Non-Overlapping Training Dataset)", colwidths = 8)
flextable_object <- hline_top(flextable_object, border = fp_border_default(width = 0), part = "header")
flextable_object <- align(flextable_object, i = 1, j = NULL, align = "left", part = "header")
flextable_object <- bold(flextable_object, bold = TRUE, part = "header")

flextable_object <- font(flextable_object, fontname = "Cambria", part = "all")
flextable_object <- fontsize(flextable_object, size = 11, part = "all")
flextable_object <- fontsize(flextable_object, size = 12, part = "header")
flextable_object <- fontsize(flextable_object, size = 10, part = "footer")

metrics <- as.data.frame(metrics)

for(i in 2:ncol(metrics)){
  
  if(colnames(metrics)[i] == 'Precision'){
    toHighligth <- which.max(as.numeric(unlist(metrics[,i])))
  }
  else if(colnames(metrics)[i] == 'Recall'){
    toHighligth <- which.max(as.numeric(unlist(metrics[,i])))
  }
  else if(colnames(metrics)[i] == 'Accuracy'){
    toHighligth <- which.max(as.numeric(unlist(metrics[,i])))
  }
  else if(colnames(metrics)[i] == 'False_Positive_Rate'){
    toHighligth <- which.min(as.numeric(unlist(metrics[,i])))
  }
  else if(colnames(metrics)[i] == 'False_Negative_Rate'){
    toHighligth <- which.min(as.numeric(unlist(metrics[,i])))
  }
  else if(colnames(metrics)[i] == 'F1_Score'){
    toHighligth <- which.max(as.numeric(unlist(metrics[,i])))
  }
  else if(colnames(metrics)[i] == 'MCC'){
    toHighligth <- which.max(as.numeric(unlist(metrics[,i])))
  }
  flextable_object %<>% flextable::highlight(toHighligth, i, color = 'wheat')
}

save_as_image(flextable_object, path = "Table 5.png", width = 16, height = 6, units = "in", res = 1000)
flextable::save_as_html(flextable_object, path = "Table 5.html")

##############################################
############ Supplementary Figure 1 ##########
##############################################

variantDistribution <- readRDS('resources/variant_distributions.rds')

variantDistribution <- variantDistribution[protName != 'TRIO_1',]
variantDistribution_Normalized <- variantDistribution[,c('Pathogenic_Rest_Norm', 'Benign_Rest_Norm', 'Pathogenic_DH_PH_Norm', 'Benign_DH_PH_Norm')]
variantDistribution_Normalized <- as.data.table(t(as.data.frame(colMeans(variantDistribution_Normalized))))
variantDistribution_Normalized[, protName := 'All']

variantDistribution_TRIO <- variantDistribution[protName == 'TRIO',]
variantDistribution_TRIO <- variantDistribution_TRIO[,c('protName', 'Pathogenic_Rest_Norm', 'Benign_Rest_Norm', 'Pathogenic_DH_PH_Norm', 'Benign_DH_PH_Norm')]

variantDistribution <- rbind(variantDistribution_TRIO, variantDistribution_Normalized)

variantDistribution.mlt <- melt(variantDistribution, id.vars = 'protName', 
                                measure.vars = c('Benign_Rest_Norm', 'Benign_DH_PH_Norm', 'Pathogenic_Rest_Norm', 'Pathogenic_DH_PH_Norm'))

colnames(variantDistribution.mlt) <- c('protName', 'Variant_Type', 'Number_of_Variants_Normalized')

variantDistribution.mlt$onDomain <- TRUE
variantDistribution.mlt[grepl('Rest', Variant_Type),]$onDomain <- FALSE

variantDistribution.mlt$Type <- 'Benign'
variantDistribution.mlt[grepl('Pathogenic', Variant_Type),]$Type <- 'Pathogenic'

variantDistribution.mlt[Variant_Type == "Pathogenic_Rest_Norm" & protName != 'TRIO',]$Variant_Type <- "All Proteins with DH/PH\nOther domains"
variantDistribution.mlt[Variant_Type == "Pathogenic_DH_PH_Norm" & protName != 'TRIO',]$Variant_Type <- "All Proteins with DH/PH\nDH or PH domains"
variantDistribution.mlt[Variant_Type == "Benign_Rest_Norm" & protName != 'TRIO',]$Variant_Type <- "All Proteins with DH/PH\nOther domains"
variantDistribution.mlt[Variant_Type == "Benign_DH_PH_Norm" & protName != 'TRIO',]$Variant_Type <- "All Proteins with DH/PH\nDH or PH domains"

variantDistribution.mlt[Variant_Type == "Pathogenic_Rest_Norm" & protName == 'TRIO',]$Variant_Type <- "TRIO Only\nOther domains"
variantDistribution.mlt[Variant_Type == "Pathogenic_DH_PH_Norm" & protName == 'TRIO',]$Variant_Type <- "TRIO Only\nDH or PH domains"
variantDistribution.mlt[Variant_Type == "Benign_Rest_Norm" & protName == 'TRIO',]$Variant_Type <- "TRIO Only\nOther domains"
variantDistribution.mlt[Variant_Type == "Benign_DH_PH_Norm" & protName == 'TRIO',]$Variant_Type <- "TRIO Only\nDH or PH domains"

p_All <- ggplot(data=variantDistribution.mlt, aes(x = Variant_Type, y = Number_of_Variants_Normalized, fill = Type)) +
  geom_bar(position = "dodge", stat="identity") + xlab('Region') + ylab('Number of Variants per Residue') + 
  scale_fill_brewer(palette="Accent") +
  scale_x_discrete(labels=c("TRIO Only\nOther domains" = "Other regions", "TRIO Only\nDH or PH domains" = "DH or PH domains",
                            "All Proteins with DH/PH\nOther domains" = "Other regions", "All Proteins with DH/PH\nDH or PH domains" = "DH or PH domains")) +
  geom_label(aes(x = 1.5, y = 0.015,
                 label = "Human Proteins\nwith DH/PH domains"),
             stat = "unique", fill = 'cornsilk3',
             size = 6) +
  geom_label(aes(x = 3.5, y = 0.015,
                 label = "TRIO"),
             stat = "unique", fill = 'cornsilk3',
             size = 7) +
  labs(fill="Variant Class") + 
  ggtitle('Variant localisation in human proteins containing a DH domain') + 
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size = 14))

tiff('Figure S1.tiff', units="in", width=10, height=7, res=1000)
p_All
dev.off()

png('Figure S1.png', units="in", width=10, height=7, res=1000)
p_All
dev.off()

##############################################
############ Supplementary Figure 2 ##########
##############################################

All_VarDB_VariantCounts <- readRDS('resources/All_VarDB_VariantCounts.rds')
TrioNlike_RegionLengths <- readRDS('resources/TrioNlike_RegionLengths.rds')

All_VarDB_VariantCounts_Plot <- unique(All_VarDB_VariantCounts[,.(target_id, ClinicalSignificance, VarMapped2Domain)])
All_VarDB_VariantCounts_Plot[, Accession := sapply(strsplit(target_id, '_', fixed = T), '[', 1)]
All_VarDB_VariantCounts_Plot[, target_id := NULL]
All_VarDB_VariantCounts_Plot <- merge(All_VarDB_VariantCounts_Plot, TrioNlike_RegionLengths, by = 'Accession')
All_VarDB_VariantCounts_Plot[, DH_Domains_Length := NULL]
All_VarDB_VariantCounts_Plot[, PH_Domains_Length := NULL]
All_VarDB_VariantCounts_Plot[, Non_DHPH_Length := NULL]

result <- All_VarDB_VariantCounts_Plot[, .(
  Benign_DH = sum(ClinicalSignificance == "Benign" & VarMapped2Domain == "DH"),
  Benign_PH = sum(ClinicalSignificance == "Benign" & VarMapped2Domain == "PH"),
  Benign_Not_DH_and_PH = sum(ClinicalSignificance == "Benign" & VarMapped2Domain == "Not_DH_and_PH"),
  Pathogenic_DH = sum(ClinicalSignificance == "Pathogenic" & VarMapped2Domain == "DH"),
  Pathogenic_PH = sum(ClinicalSignificance == "Pathogenic" & VarMapped2Domain == "PH"),
  Pathogenic_Not_DH_and_PH = sum(ClinicalSignificance == "Pathogenic" & VarMapped2Domain == "Not_DH_and_PH")
), by = Accession]
result <- merge(result, TrioNlike_RegionLengths[,.(Accession, DH_Domains_Length, PH_Domains_Length, Non_DHPH_Length)], by = 'Accession')

normalized_result <- result[, .(
  Accession = Accession,
  Normalized_Benign_DH = Benign_DH / DH_Domains_Length,
  Normalized_Benign_PH = Benign_PH / PH_Domains_Length,
  Normalized_Benign_Not_DH_and_PH = Benign_Not_DH_and_PH / Non_DHPH_Length,
  Normalized_Pathogenic_DH = Pathogenic_DH / DH_Domains_Length,
  Normalized_Pathogenic_PH = Pathogenic_PH / PH_Domains_Length,
  Normalized_Pathogenic_Not_DH_and_PH = Pathogenic_Not_DH_and_PH / Non_DHPH_Length
)]

merged_counts <- as.data.table(colSums(normalized_result[,2:7]))/12
merged_counts[, Name := colnames(normalized_result)[2:7]]
colnames(merged_counts)[1] <- 'Number_of_Variants_Normalized'
merged_counts[, Name := gsub('Normalized_', '', Name, fixed = T)]
merged_counts[, Location := gsub('Benign_', '', Name, fixed = T)]
merged_counts[, Location := gsub('Pathogenic_', '', Location, fixed = T)]
merged_counts[, Class := sapply(strsplit(Name, '_', fixed = T), '[', 1)]
merged_counts[, Name := NULL]
merged_counts[Location == 'Not_DH_and_PH']$Location <- 'Other Regions'

p_All <- ggplot(data=merged_counts, aes(x = Location, y = Number_of_Variants_Normalized, fill = Class)) +
  geom_bar(position = "dodge", stat="identity") + xlab('Region') + ylab('Number of Variants per Residue') + 
  scale_fill_brewer(palette="Accent") +
  labs(fill="Variant Class") + 
  ggtitle('Variant localisation in human proteins containing a TrioN-like DH domain') + 
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size = 14))

tiff('Figure S2.tiff', units="in", width=10, height=7, res=1000)
p_All
dev.off()

png('Figure S2.png', units="in", width=10, height=7, res=1000)
p_All
dev.off()

##############################################
############ Supplementary Figure 3 ##########
##############################################

## Figure created using the PyMOL software

##############################################
############ Supplementary Figure 4 ##########
##############################################

## Please refer to https://github.com/alptaciroglu/drawProteins for the code to generate the image

##############################################
############ Supplementary Movie 1 ###########
##############################################

## Figure created using the PyMOL software

##############################################
############ Supplementary Movie 2 ###########
##############################################

## Figure created using the PyMOL software

##############################################
############ Supplementary Table 1 ###########
##############################################

SSAP_Human_DH <- readRDS('resources/SSAP_Human_DH_straln.rds')

flextable_object <- flextable(SSAP_Human_DH) %>% 
  autofit() %>% 
  theme_booktabs() %>% 
  bold(part = "header") %>%
  align(align = "center", part = "all") 

flextable_object <- add_footer_lines(flextable_object, "Target: TrioN DH Domain (TRIO_1284_1470)")
flextable_object <- add_footer_lines(flextable_object, "Target Length: 187 residues")

flextable_object <- add_header_row(flextable_object, values = "Supplementary Table 1: Structural Alignment of Human DH Domains to TrioN DH Domain", colwidths = 8)
flextable_object <- hline_top(flextable_object, border = fp_border_default(width = 0), part = "header")

flextable_object <- align(flextable_object, i = 1, j = NULL, align = "left", part = "header")
# Make the "Methods" header bold
flextable_object <- bold(flextable_object, bold = TRUE, part = "header")
flextable_object <- font(flextable_object, fontname = "Cambria", part = "all")
flextable_object <- fontsize(flextable_object, size = 11, part = "all")
flextable_object <- fontsize(flextable_object, size = 12, part = "header")
flextable_object <- fontsize(flextable_object, size = 10, part = "footer")
flextable_object <- line_spacing(flextable_object, i = NULL, j = NULL, space = 0.5, part = "body")

save_as_image(flextable_object, path = "Table S1.png", units = "in", res = 1000)
save_as_html(flextable_object, path = "Table S1.html")

##############################################
############ Supplementary Figure 5 ##########
##############################################

RMSD_Values_SSAP_HumanDH <- read.xlsx('resources/SSAP_Human_DH_straln.xlsx')
RMSD_Values_SSAP_HumanDH <- as.data.table(RMSD_Values_SSAP_HumanDH)

RMSD_Values_SSAP_HumanDH$protName[1:2] <- c('TrioN', 'KalrnN')
RMSD_Values_SSAP_HumanDH$protName[8] <- 'TrioC'
RMSD_Values_SSAP_HumanDH$protName[9] <- 'KalrnC'

sp <- ggplot(RMSD_Values_SSAP_HumanDH, aes(x=RMSd, y=Equivalenced_Perc, color=RMSd)) + 
  geom_point() + xlab('RMSD') + ylab('Equivalenced Percentages') + geom_text_repel(aes(label = protName), size = 3.5) +
  scale_y_continuous(limits = c(70, 100)) +
  scale_x_reverse(limits = c(5, 0)) + 
  ggtitle('Structural Similarity of Human DH Domains (Reference - TrioN)') +
  theme(text = element_text(size=15),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = 'none')

tiff('Figure S5.tiff', units="in", width = 12, height = 6, res=1000)
sp
dev.off()

png('Figure S5.png', units="in", width = 12, height = 6, res=1000)
sp
dev.off()

##############################################
############ Supplementary Table 2 ###########
##############################################

model_performances <- readRDS('predictor/models_all_methods_all_formulas_performances.rds')
metrics <- model_performances

# Create a visually appealing table with flextable
flextable_object <- flextable(metrics) %>% 
  autofit() %>% 
  theme_booktabs() %>% 
  bold(part = "header") %>%
  align(align = "center", part = "all") 

flextable_object <- add_header_row(flextable_object, values = "Supplementary Table 2: Performance evaluation of 18 ML Methods tested for method and feature optimisation", colwidths = 8)
flextable_object <- hline_top(flextable_object, border = fp_border_default(width = 0), part = "header")

# Left align the "Methods" header
flextable_object <- align(flextable_object, i = 1, j = NULL, align = "left", part = "header")
# Make the "Methods" header bold
flextable_object <- bold(flextable_object, bold = TRUE, part = "header")
flextable_object <- font(flextable_object, fontname = "Cambria", part = "all")
flextable_object <- fontsize(flextable_object, size = 11, part = "all")
flextable_object <- fontsize(flextable_object, size = 12, part = "header")
flextable_object <- fontsize(flextable_object, size = 10, part = "footer")

metrics <- as.data.frame(metrics)

for(i in 2:ncol(metrics)){
  
  if(colnames(metrics)[i] == 'Precision'){
    toHighligth <- which.max(as.numeric(unlist(metrics[,i])))
  }
  else if(colnames(metrics)[i] == 'Recall'){
    toHighligth <- which.max(as.numeric(unlist(metrics[,i])))
  }
  else if(colnames(metrics)[i] == 'Accuracy'){
    toHighligth <- which.max(as.numeric(unlist(metrics[,i])))
  }
  else if(colnames(metrics)[i] == 'False_Positive_Rate'){
    toHighligth <- which.min(as.numeric(unlist(metrics[,i])))
  }
  else if(colnames(metrics)[i] == 'False_Negative_Rate'){
    toHighligth <- which.min(as.numeric(unlist(metrics[,i])))
  }
  else if(colnames(metrics)[i] == 'F1_Score'){
    toHighligth <- which.max(as.numeric(unlist(metrics[,i])))
  }
  else if(colnames(metrics)[i] == 'MCC'){
    toHighligth <- which.max(as.numeric(unlist(metrics[,i])))
  }
  flextable_object %<>% flextable::highlight(toHighligth, i, color = 'wheat')
}

save_as_image(flextable_object, path = "Table S2.png", width = 16, height = 96, units = "in", res = 1000)
save_as_html(flextable_object, path = "Table S2.html")

##############################################
############ Supplementary Figure 6 ##########
##############################################

metrics <- readRDS('predictor/models_all_methods_multipleruns_performances.rds')

data_multitrain_metrics <- metrics
data_multitrain_metrics <- as.data.table(data_multitrain_metrics)
data_multitrain_metrics[, c("Method", "Formula") := tstrsplit(Predictor, "~", fixed = TRUE)]

data_multitrain_metrics[is.na(Formula), Formula := "None"]

data_multitrain_metrics$Formula <- factor(data_multitrain_metrics$Formula)

data_multitrain_metrics$Method <- factor(data_multitrain_metrics$Method, levels = c(
  sort(unique(data_multitrain_metrics$Method))
))

ML_multiTrain_Methods_MCC <- ggplot(data_multitrain_metrics, aes(x = Method, y = MCC, fill = Method)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "MCC Score Distribution Across 20 Training Runs for 18 ML Algorithms in Variant Effect Prediction", 
    x = "Method", 
    y = "MCC"
  ) +
  theme(
    title = element_text(size = 15),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),  # Increase x-axis text size
    axis.title.x = element_text(size = 18),                        # Increase x-axis title size
    axis.text.y = element_text(size = 16),                         # Increase y-axis text size
    axis.title.y = element_text(size = 18),                         # Increase y-axis title size
    legend.position = "none"                                       # Remove legend
  ) +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.1))  # Increase the frequency of y-axis ticks

tiff("Figure_S6.tiff", width = 12, height = 8, units = "in", res = 1000)
print(ML_multiTrain_Methods_MCC)
dev.off()

png("Figure_S6.png", width = 12, height = 8, units = "in", res = 1000)
print(ML_multiTrain_Methods_MCC)
dev.off()

##############################################
############ Supplementary Figure 7 ##########
##############################################

model_performances <- readRDS('predictor/models_all_methods_all_formulas_performances.rds')
metrics <- model_performances

metrics[, c("Method", "Formula") := tstrsplit(Predictor, "~", fixed = TRUE)]

metrics[is.na(Formula), Formula := "None"]
metrics$Formula[1:7] <- metrics$Predictor[1:7]

metrics$Formula <- factor(metrics$Formula, levels = c(
  "Inclusive", "NoAlphaMissense", "NofoldX", "NovarSite", "NoVariPred", "NoPHACT", "NoPHACTboost", "NoSAAPdap", "NoVariantClusters", "Exclusionary", 
  unique(metrics$Formula[!(metrics$Formula %in% c("Inclusive", "NoAlphaMissense", "NofoldX", "NovarSite", "NoVariPred", "NoPHACT", "NoPHACTboost", "NoSAAPdap", "NoVariantClusters", "Exclusionary"))])
))

new_labels <- c("Inclusive" = "Inclusive", "NoAlphaMissense" = "w/o AlphaMissense", "NofoldX" = "w/o foldX", "NovarSite" = "w/o varSite", "NoVariPred" = "w/o VariPred", "NoPHACT" = "w/o PHACT", "NoPHACTboost" = 'w/o PHACTboost', "NoSAAPdap" = "w/o SAAPdap", "NoVariantClusters" = "w/o Variant Clusters", "Exclusionary" = "w/o External\nPredictors")

metrics <- metrics[-c(1:7),]

# Now, create the updated plot with the reordered factor levels and new x-axis labels
ML_Methods_MCC_Formula <- ggplot(metrics, aes(x = Formula, y = MCC, fill = Formula)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Performance of ML Methods by Formula through MCC", 
    x = "Formula", 
    y = "MCC"
  ) +
  theme(
    title = element_text(size = 20),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),  # Increase x-axis text size
    axis.title.x = element_text(size = 18),                        # Increase x-axis title size
    axis.text.y = element_text(size = 16),                         # Increase y-axis text size
    axis.title.y = element_text(size = 18),                         # Increase y-axis title size
    legend.position = "none"                                       # Remove legend
  ) +
  scale_x_discrete(labels = new_labels) +  # Modify x-axis labels
  scale_y_continuous(breaks = seq(-1, 1, by = 0.1)) +  # Increase the frequency of y-axis ticks
  stat_summary(
    fun.min = min, 
    fun.max = max, 
    geom = "text", 
    aes(label = sprintf("Min: %.2f\nMax: %.2f", after_stat(ymin), after_stat(ymax))),  # Corrected notation
    position = position_dodge(width = 0.75),
    vjust = -0.5  # Single value for vertical adjustment to avoid mismatch error
  )

tiff("Figure S7.tiff", width = 12, height = 8, units = "in", res = 1000)
ML_Methods_MCC_Formula
dev.off()

png("Figure S7.png", width = 12, height = 8, units = "in", res = 1000)
ML_Methods_MCC_Formula
dev.off()

##############################################
############ Supplementary Figure 8 ##########
##############################################

model_performances <- readRDS('predictor/models_all_methods_all_formulas_performances.rds')
metrics <- model_performances
metrics <- metrics[-c(1:7),]
metrics[, Method := sapply(strsplit(Predictor, '~', fixed = T), '[', 1)]

metrics$Method <- factor(metrics$Method, levels = c(
  "avNNet", "bagFDA", "C5.0", "earth", "gbm", "glmboost", "glmnet", "knn", "naive_bayes", "nnet", "pcaNNet", "ranger", "rf", "rpart", "svmLinear", "svmRadial", "xgbLinear", "xgbTree", 
  unique(metrics$Method[!(metrics$Method %in% c("avNNet", "bagFDA", "C5.0", "earth", "gbm", "glmboost", "glmnet", "knn", "naive_bayes", "nnet", "pcaNNet", "ranger", "rf", "rpart", "svmLinear", "svmRadial", "xgbLinear", "xgbTree"))])
))

ML_Methods_MCC <- ggplot(metrics, aes(x = Method, y = MCC, fill = Method)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Performance of ML Methods by MCC", 
    x = "Method", 
    y = "MCC"
  ) +
  theme(
    title = element_text(size = 20),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),  # Increase x-axis text size
    axis.title.x = element_text(size = 18),                        # Increase x-axis title size
    axis.text.y = element_text(size = 16),                         # Increase y-axis text size
    axis.title.y = element_text(size = 18),                         # Increase y-axis title size
    legend.position = "none"                                       # Remove legend
  ) +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.1))  # Increase the frequency of y-axis ticks

tiff("Figure S8.tiff", width = 12, height = 8, units = "in", res = 1000)
ML_Methods_MCC
dev.off()

png("Figure S8.png", width = 12, height = 8, units = "in", res = 1000)
ML_Methods_MCC
dev.off()
