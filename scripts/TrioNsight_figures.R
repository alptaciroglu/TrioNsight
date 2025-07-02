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
if (!requireNamespace("magrittr", quietly = TRUE)) {
  BiocManager::install('magrittr')
}
library(magrittr)
if (!requireNamespace("officer", quietly = TRUE)) {
  BiocManager::install('officer')
}
library(officer)
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  BiocManager::install('ggrepel')
}
library(ggrepel)

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

## Figure created using the PyMOL software

##############################################
################# Figure 4 ###################
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

# save_as_image(flextable_object, path = "Table 1.png", width = 12, height = 8, units = "in", res = 1000)
flextable::save_as_html(flextable_object, path = "Table 1.html")
write.xlsx(as.data.frame(matrix_data), file = 'Table 1.xlsx', asTable = TRUE, overwrite = TRUE)

##############################################
################# Figure 5 ###################
##############################################

## Figure created using the draw.io web application 

##############################################
################# Figure 6 ###################
##############################################

## Figure created using the PyMOL software

##############################################
################# Figure 7 ###################
##############################################

## Figure created using the PyMOL software

##############################################
################# Figure 8 ###################
##############################################

## Figure created using the PyMOL software

##############################################
################# Table 2.1 ##################
##############################################

best_n_analysis <- readRDS("stable_run/best_n_analysis.rds")

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
colnames(best_clustering) <- c('p-values', 'ꭓ-squared')

flextable_object <- flextable(best_clustering %>% tibble::rownames_to_column("Number of Clusters")) %>% 
  autofit() %>% 
  theme_booktabs() %>% 
  bold(part = "header") %>%
  align(align = "center", part = "all")

flextable_object <- add_header_row(flextable_object, values = "Table 2.1: Best clustering by ꭓ-squared test", colwidths = 3)
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

# toHighlight <- which.min(as.numeric(unlist(best_clustering[,1])))
# ## flextable_object %<>% flextable::highlight(toHighlight, 2, color = 'wheat')
# flextable_object %<>% flextable::bold(toHighlight, 2)

# save_as_image(flextable_object, path = "Table 2.1.png", width = 12, height = 6, units = "in", res = 1000)
flextable::save_as_html(flextable_object, path = "Table 2.1.html")

best_clustering <- as.data.table(cbind("Number of Clusters" = rownames(best_clustering), best_clustering))
best_clustering$`Number of Clusters` <- as.numeric(best_clustering$`Number of Clusters`)
write.xlsx(best_clustering, file = 'Table 2_1.xlsx', asTable = TRUE, overwrite = TRUE)

##############################################
################# Table 2.2 ##################
##############################################

TrioN_DHPH_1290_1594_A_1nty_strucclus_n5 <- fread('stable_run/strucclus_res/TrioN_DHPH_1290_1594_A_1nty_n5.txt', sep = '\t', skip = 62, header = F)

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

# save_as_image(flextable_object, path = "Table 2.2.png", width = 12, height = 6, units = "in", res = 1000)
flextable::save_as_html(flextable_object, path = "Table 2.2.html")
write.xlsx(clusteringTable, file = 'Table 2_2.xlsx', asTable = TRUE, overwrite = TRUE)

##############################################
################# Table 3 ####################
##############################################

model_performances <- readRDS('predictor/performance_all_models_all_methods.rds')

metrics <- model_performances
selected_methods <- c('CADD', 'ESM-1b', 'AlphaMissense', 'VariPred', 'MutPred2', 'PolyPhen-2', 'PHACT', 'PHACTboost', 'naive_bayes~Inclusive')
metrics <- metrics[Predictor %in% selected_methods,]
metrics$Predictor[9] <- 'TrioNsight'

metrics <- as.data.frame(metrics)
metrics[,-1] <- sapply(metrics[,-1], function(x){round(x,3)})
metrics <- as.data.table(metrics)

# Create a visually appealing table with flextable
flextable_object <- flextable(metrics) %>% 
  autofit() %>% 
  theme_booktabs() %>% 
  bold(part = "header") %>%
  align(align = "center", part = "all") 

flextable_object <- add_header_row(flextable_object, values = "Table 3: Benchmarking TrioNsight performance against variant pathogenicity prediction algorithms (Non-Overlapping Training Dataset)", colwidths = 8)
flextable_object <- hline_top(flextable_object, border = fp_border_default(width = 0), part = "header")

# Left align the "Methods" header
flextable_object <- align(flextable_object, i = 1, j = NULL, align = "left", part = "header")
# Make the "Methods" header bold
flextable_object <- bold(flextable_object, bold = TRUE, part = "header")
flextable_object <- font(flextable_object, fontname = "Cambria", part = "all")
flextable_object <- fontsize(flextable_object, size = 11, part = "all")
flextable_object <- fontsize(flextable_object, size = 12, part = "header")
flextable_object <- fontsize(flextable_object, size = 10, part = "footer")

for(i in 2:ncol(metrics)){
  
  metric_column <- as.numeric(unlist(metrics[,..i]))
  
  if(colnames(metrics)[i] == 'Precision'){
    toBold <- which(metric_column == min(metric_column, na.rm = TRUE))
  }
  else if(colnames(metrics)[i] == 'Recall'){
    toBold <- which(metric_column == max(metric_column, na.rm = TRUE))
  }
  else if(colnames(metrics)[i] == 'Accuracy'){
    toBold <- which(metric_column == max(metric_column, na.rm = TRUE))
  }
  else if(colnames(metrics)[i] == 'False_Positive_Rate'){
    toBold <- which(metric_column == min(metric_column, na.rm = TRUE))
  }
  else if(colnames(metrics)[i] == 'False_Negative_Rate'){
    toBold <- which(metric_column == min(metric_column, na.rm = TRUE))
  }
  else if(colnames(metrics)[i] == 'F1_Score'){
    toBold <- which(metric_column == max(metric_column, na.rm = TRUE))
  }
  else if(colnames(metrics)[i] == 'MCC'){
    toBold <- which(metric_column == max(metric_column, na.rm = TRUE))
  }
  flextable_object %<>% flextable::bold(toBold, i, bold = T, part = "body")
}

# save_as_image(flextable_object, path = "Table 3_1.png", width = 16, height = 6, units = "in", res = 1000)
save_as_html(flextable_object, path = "Table 3.html")
write.xlsx(metrics, file = 'Table 3.xlsx', asTable = TRUE, overwrite = TRUE)

##############################################
################# Figure 9 ###################
##############################################

model_performances <- readRDS('predictor/performance_all_models_all_methods.rds')

metrics <- model_performances

metrics[, c("Method", "Formula") := tstrsplit(Predictor, "~", fixed = TRUE)]

metrics[is.na(Formula), Formula := "None"]
metrics$Formula[1:7] <- metrics$Predictor[1:7]

metrics$Formula <- factor(metrics$Formula, levels = c(
  "Inclusive", "NofoldX", "NoVariantClusters", "NoSAAPdap", "NoStructuralFeatures", "NovarSite", "NoSequenceFeatures", "NoAlphaMissense", "NoVariPred", "NoPHACT",  "NoPHACTboost",  "Exclusionary", "CompPredictorsOnly",
  unique(metrics$Formula[!(metrics$Formula %in% c("Inclusive", "NofoldX", "NoVariantClusters", "NoSAAPdap", "NoStructuralFeatures", "NovarSite", "NoSequenceFeatures", "NoAlphaMissense", "NoVariPred", "NoPHACT",  "NoPHACTboost", "Exclusionary", "CompPredictorsOnly"))])
))

new_labels <- c("Inclusive" = "Inclusive", "NofoldX" = "w/o foldX", "NoVariantClusters" = "w/o Variant Clusters", "NoSAAPdap" = "w/o SAAPdap", "NoStructuralFeatures" = "w/o Structural\nFeatures", "NovarSite" = "w/o varSite", "NoSequenceFeatures" = "w/o Sequence\nFeatures", "NoAlphaMissense" = "w/o AlphaMissense", "NoVariPred" = "w/o VariPred", "NoPHACT" = 'w/o PHACT', "NoPHACTboost" = 'w/o PHACTboost', "Exclusionary" = "w/o External\nPredictors", "CompPredictorsOnly" = "Only Component\nPredictors")

metrics <- metrics[-c(1:8),]
metrics <- metrics[grepl('naive', Predictor, fixed = T)]

# Now, create the updated plot with the reordered factor levels and new x-axis labels
TrioNsight_FeatureAblation <- ggplot(metrics, aes(x = Formula, y = MCC, fill = Formula)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "TrioNsight - Feature Ablation Study", 
    x = "Feature Sets", 
    y = "MCC"
  ) +
  theme(
    title = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  
    ## axis.title.x = element_text(size = 18),                        
    axis.text.y = element_text(size = 14),                         
    ## axis.title.y = element_text(size = 18),                         
    legend.position = "none"                                       
  ) +
  scale_x_discrete(labels = new_labels) +  
  scale_y_continuous(breaks = seq(0.5, 1, by = 0.05), limits = c(0.5, 1.0)) +  
  stat_summary(
    fun = max,
    geom = "text", 
    aes(label = sprintf("%.3f", after_stat(y))),  # Display max value with 3 decimal places
    vjust = -0.5,  # Position above the whisker
    size = 4       # Adjust text size as needed
  )

# tiff("Figure S9.tiff", width = 14, height = 8, units = "in", res = 1000)
# ML_Methods_MCC_Formula
# dev.off()

png("Figure 9.png", width = 10, height = 7, units = "in", res = 1000)
TrioNsight_FeatureAblation
dev.off()

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
flextable_object <- vline(flextable_object, i = NULL, j = "Annotation", border = NULL, part = "all")

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

toBold <- which(TrioNsight_distinguishedPredictions$Annotation == TrioNsight_distinguishedPredictions$AlphaMissense)
flextable_object %<>% flextable::bold(i = toBold, j = "AlphaMissense", bold = TRUE, part = "body")
toBold <- which(TrioNsight_distinguishedPredictions$Annotation == TrioNsight_distinguishedPredictions$PHACTboost)
flextable_object %<>% flextable::bold(i = toBold, j = "PHACTboost", bold = TRUE, part = "body")
toBold <- which(TrioNsight_distinguishedPredictions$Annotation == TrioNsight_distinguishedPredictions$VariPred)
flextable_object %<>% flextable::bold(i = toBold, j = "VariPred", bold = TRUE, part = "body")
toBold <- which(TrioNsight_distinguishedPredictions$Annotation == TrioNsight_distinguishedPredictions$MutPred2)
flextable_object %<>% flextable::bold(i = toBold, j = "MutPred2", bold = TRUE, part = "body")
toBold <- which(TrioNsight_distinguishedPredictions$Annotation == TrioNsight_distinguishedPredictions$CADD)
flextable_object %<>% flextable::bold(i = toBold, j = "CADD", bold = TRUE, part = "body")
toBold <- which(TrioNsight_distinguishedPredictions$Annotation == TrioNsight_distinguishedPredictions$`ESM-1b`)
flextable_object %<>% flextable::bold(i = toBold, j = "ESM-1b", bold = TRUE, part = "body")
toBold <- which(TrioNsight_distinguishedPredictions$Annotation == TrioNsight_distinguishedPredictions$`PolyPhen-2`)
flextable_object %<>% flextable::bold(i = toBold, j = "PolyPhen-2", bold = TRUE, part = "body")
toBold <- which(TrioNsight_distinguishedPredictions$Annotation == TrioNsight_distinguishedPredictions$TrioNsight)
flextable_object %<>% flextable::bold(i = toBold, j = "TrioNsight", bold = TRUE, part = "body")

# save_as_image(flextable_object, path = "Table 4.png", width = 16, height = 8, units = "in", res = 1000)
flextable::save_as_html(flextable_object, path = "Table 4.html")
write.xlsx(TrioNsight_distinguishedPredictions, file = 'Table 4.xlsx', asTable = TRUE, overwrite = TRUE)

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

variantDistribution.mlt[Variant_Type == "Pathogenic_Rest_Norm" & protName == 'TRIO',]$Variant_Type <- "TrioN Only\nOther domains"
variantDistribution.mlt[Variant_Type == "Pathogenic_DH_PH_Norm" & protName == 'TRIO',]$Variant_Type <- "TrioN Only\nDH or PH domains"
variantDistribution.mlt[Variant_Type == "Benign_Rest_Norm" & protName == 'TRIO',]$Variant_Type <- "TrioN Only\nOther domains"
variantDistribution.mlt[Variant_Type == "Benign_DH_PH_Norm" & protName == 'TRIO',]$Variant_Type <- "TrioN Only\nDH or PH domains"

p_All <- ggplot(data=variantDistribution.mlt, aes(x = Variant_Type, y = Number_of_Variants_Normalized, fill = Type)) +
  geom_bar(position = "dodge", stat="identity") + xlab('Region') + ylab('Number of Variants per Residue') + 
  scale_fill_brewer(palette="Accent") +
  scale_x_discrete(labels=c("TrioN Only\nOther domains" = "Other regions", "TrioN Only\nDH or PH domains" = "DH or PH domains",
                            "All Proteins with DH/PH\nOther domains" = "Other regions", "All Proteins with DH/PH\nDH or PH domains" = "DH or PH domains")) +
  geom_label(aes(x = 1.5, y = 0.015,
                 label = "Human Proteins\nwith DH/PH domains"),
             stat = "unique", fill = NA,
             size = 6) +
  geom_segment(aes(x = 0.55, xend = 0.85, y = 0.015, yend = 0.015),
               colour = "black", linetype ="solid", size = 1) +
  geom_segment(aes(x = 2.15, xend = 2.45, y = 0.015, yend = 0.015),
               colour = "black", linetype ="solid", size = 1) +
  geom_label(aes(x = 3.5, y = 0.015,
                 label = "TrioN"),
             stat = "unique", fill = NA,
             size = 7) +
  geom_segment(aes(x = 2.50, xend = 3.25, y = 0.015, yend = 0.015),
               colour = "black", linetype ="dashed", size = 1) +
  geom_segment(aes(x = 3.75, xend = 4.45, y = 0.015, yend = 0.015),
               colour = "black", linetype ="dashed", size = 1) +
  labs(fill="Variant Class") + 
  ggtitle('Variant localisation in human proteins containing a DH domain') + 
  theme_minimal() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

# tiff('Figure S3.tiff', units="in", width=10, height=7, res=1000)
# p_All
# dev.off()

png('Figure S1.png', units="in", width=10, height=7, res=1000)
p_All
dev.off()

##############################################
############ Supplementary Figure 2 ##########
##############################################

## Figure created using the draw.io web application 

##############################################
############ Supplementary Figure 3 ##########
##############################################

## Figure created using the draw.io web application 

##############################################
############ Supplementary Table 1 ###########
##############################################

## Generate TrioNsight variant impact map ##

TrioNsight_features_dt <- fread('resources/TrioNsight_features_dataset.txt', sep = '\t')
TrioNsight_features_preprocessed_dt <- preprocess_training_data(TrioNsight_features_dt)

unique_locs <- unique(TrioNsight_features_preprocessed_dt$MSALocation)
TrioNsight_features_preprocessed_dt$target_id <- factor(TrioNsight_features_preprocessed_dt$target_id)

TrioNsight_impactMap_features <- readRDS('resources/TrioNsight_impactMap_features.rds')
TrioNsight_impactMap_features_cleaned <- preprocess_training_data(TrioNsight_impactMap_features)

preTuned_dt <- readRDS('resources/fineTuningParams_best.rds')

impactmap_predictionslist <- list()

for (k in 1:length(unique_locs)) {
  
  track_progress(k, length(unique_locs))
  # Current MSA location for testing
  current_loc <- unique_locs[k]
  
  # Training data: all except current MSA location
  train_data <- TrioNsight_features_preprocessed_dt[TrioNsight_features_preprocessed_dt$MSALocation != current_loc, ]
  training_variants <- unique(train_data$target_id)
  
  # Testing data: only non-overlapping variants from current MSA location
  test_data <- TrioNsight_impactMap_features_cleaned[!target_id %in% training_variants]
  testing_variants <- test_data$target_id
  
  # Skip if no non-overlapping variants to test
  if(nrow(test_data) == 0) next
  
  train_data[, target_id := NULL]
  test_data[, target_id := NULL]
  train_data[, MSALocation := NULL]
  test_data[, MSALocation := NULL]
  train_data[, Inherited_Seq := NULL]
  test_data[, Inherited_Seq := NULL]
  
  modelTrained <- train(Class ~ . ,
                        data = train_data,
                        method = 'naive_bayes',
                        preProc = c("center", "scale"),
                        trControl = trainControl(method = "none", classProbs = TRUE),
                        allowParallel = TRUE,
                        tuneGrid = preTuned_dt[['naive_bayes']],
                        metric = "ROC",  # Change this if not a classification problem
                        maximize = TRUE,  # TRUE if using accuracy or other metrics where higher is better
                        verbose = TRUE)
  
  
  
  impactmap_predictionslist[[k]] <- predict(modelTrained, newdata = test_data, type = "prob")
  impactmap_predictionslist[[k]] <- as.data.table(impactmap_predictionslist[[k]])
  impactmap_predictionslist[[k]][, target_id := testing_variants]
  
  
}

TrioNsight_impactMap_predicted <- do.call(rbind, impactmap_predictionslist)

# Calculate mean probabilities for each variant
TrioNsight_impactMap_predicted <- TrioNsight_impactMap_predicted[, .(
  Benign = mean(Benign),
  Pathogenic = mean(Pathogenic),
  sd = sd(Pathogenic)  # Standard deviation to check consistency
), by = target_id]

# Sort by pathogenic probability
setorder(TrioNsight_impactMap_predicted, -Pathogenic)

# Calculate some potential cutoffs
cutoff_stats <- TrioNsight_impactMap_predicted[, .(
  q90 = quantile(Pathogenic, 0.90),
  q95 = quantile(Pathogenic, 0.95),
  q99 = quantile(Pathogenic, 0.99)
)]

# Look at variants with high disagreement between predictions
inconsistent_preds <- TrioNsight_impactMap_predicted[sd > 0.2]

# Create confidence levels
TrioNsight_impactMap_predicted[, Prediction := case_when(
  sd < 0.05 & Pathogenic > 0.75 ~ "Pathogenic",
  sd >= 0.05 & Pathogenic > 0.75 ~ "Likely pathogenic",
  sd < 0.05 & Pathogenic < 0.25 ~ "Benign",
  sd >= 0.05 & Pathogenic < 0.25 ~ "Likely benign",
  sd > 0.2 ~ "Uncertain",
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

## fwrite(TrioNsight_impactMap_predicted, file = 'TrioNsight_all_predictions.txt', sep = '\t')
write.xlsx(TrioNsight_impactMap_predicted, file = 'Table S1.xlsx', asTable = TRUE, overwrite = TRUE)

##############################################
############ Supplementary Figure 4 ##########
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
  theme_minimal() +  
  theme(text = element_text(size=16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

# tiff('Figure S4.tiff', units="in", width=10, height=7, res=1000)
# p_All
# dev.off()

png('Figure S4.png', units="in", width=10, height=7, res=1000)
p_All
dev.off()

##############################################
############ Supplementary Figure 5 ##########
##############################################

## Please refer to https://github.com/alptaciroglu/drawProteins for the code to generate the image

##############################################
############ Supplementary Movie 1 ###########
##############################################

## Figure created using the PyMOL software

##############################################
############ Supplementary Figure 6 ##########
##############################################

## Figure created using the PyMOL software

##############################################
############ Supplementary Figure 7 ##########
##############################################

RMSD_Values_SSAP_HumanDH <- readRDS('resources/SSAP_Human_DH_straln.rds')
RMSD_Values_SSAP_HumanDH <- as.data.table(RMSD_Values_SSAP_HumanDH)

RMSD_Values_SSAP_HumanDH[, protName := sapply(strsplit(Protein, ' ', fixed = T), '[', 1)]
RMSD_Values_SSAP_HumanDH$protName[1:2] <- c('TrioN', 'KalrnN')
RMSD_Values_SSAP_HumanDH$protName[8] <- 'TrioC'
RMSD_Values_SSAP_HumanDH$protName[9] <- 'KalrnC'

sp <- ggplot(RMSD_Values_SSAP_HumanDH, aes(x=RMSD, y=`Equivalenced (%)`, color=RMSD)) + 
  geom_point() + xlab('RMSD') + ylab('Equivalenced Percentages') + geom_text_repel(aes(label = protName), size = 3.5) +
  scale_y_continuous(limits = c(70, 100)) +
  scale_x_reverse(limits = c(5, 0)) + 
  ggtitle('Structural Similarity of Human DH Domains (Reference - TrioN)') +
  theme_minimal() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = 'none')

# tiff('Figure S7.tiff', units="in", width = 12, height = 6, res=1000)
# sp
# dev.off()

png('Figure S7.png', units="in", width = 10, height = 7, res=1000)
sp
dev.off()

##############################################
############ Supplementary Table 2 ###########
##############################################

TrioN_like_DH_Domains <- readRDS('resources/TrioN_like_DH_Domains.rds')
write.xlsx(TrioN_like_DH_Domains, file = 'Table S2.xlsx', asTable = TRUE, overwrite = TRUE)

##############################################
############ Supplementary Table 3 ###########
##############################################

SSAP_Human_DH <- readRDS('resources/SSAP_Human_DH_straln.rds')

flextable_object <- flextable(SSAP_Human_DH) %>% 
  autofit() %>% 
  theme_booktabs() %>% 
  bold(part = "header") %>%
  align(align = "center", part = "all") 

flextable_object <- add_footer_lines(flextable_object, "Target: TrioN DH Domain (TRIO_1284_1470)")
flextable_object <- add_footer_lines(flextable_object, "Target Length: 187 residues")

flextable_object <- add_header_row(flextable_object, values = "Supplementary Table 3: Structural Alignment of Human DH Domains to TrioN DH Domain", colwidths = 8)
flextable_object <- hline_top(flextable_object, border = fp_border_default(width = 0), part = "header")

flextable_object <- align(flextable_object, i = 1, j = NULL, align = "left", part = "header")
# Make the "Methods" header bold
flextable_object <- bold(flextable_object, bold = TRUE, part = "header")
flextable_object <- font(flextable_object, fontname = "Cambria", part = "all")
flextable_object <- fontsize(flextable_object, size = 11, part = "all")
flextable_object <- fontsize(flextable_object, size = 12, part = "header")
flextable_object <- fontsize(flextable_object, size = 10, part = "footer")
flextable_object <- line_spacing(flextable_object, i = NULL, j = NULL, space = 0.5, part = "body")

# save_as_image(flextable_object, path = "Table S3.png", units = "in", res = 1000)
save_as_html(flextable_object, path = "Table S3.html")
write.xlsx(SSAP_Human_DH, file = 'Table S3.xlsx', asTable = TRUE, overwrite = TRUE)

##############################################
############ Supplementary Movie 2 ###########
##############################################

## Figure created using the PyMOL software

##############################################
############ Supplementary Figure 8 ##########
##############################################

metrics <- readRDS('predictor/performance_all_models_multipleruns.rds')

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
    title = "Distribution of MCC Scores from 20 Training Iterations for 18 ML Methods", 
    x = "Method", 
    y = "MCC"
  ) +
  theme(
    title = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # Increase x-axis text size
    ## axis.title.x = element_text(size = 18),                        # Increase x-axis title size
    axis.text.y = element_text(size = 14),                         # Increase y-axis text size
    ## axis.title.y = element_text(size = 18),                         # Increase y-axis title size
    legend.position = "none"                                       # Remove legend
  ) +
  scale_y_continuous(breaks = seq(0.5, 1, by = 0.05), limits = c(0.5, 1.0))

# tiff("Figure S8.tiff", width = 12, height = 8, units = "in", res = 1000)
# print(ML_multiTrain_Methods_MCC)
# dev.off()

png("Figure S8.png", width = 10, height = 7, units = "in", res = 1000)
print(ML_multiTrain_Methods_MCC)
dev.off()

##############################################
########### Supplementary Figure 9 ##########
##############################################

model_performances <- readRDS('predictor/performance_all_models_all_methods.rds')
metrics <- model_performances
metrics <- metrics[-c(1:8),]
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
    title = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # Increase x-axis text size
    ## axis.title.x = element_text(size = 18),                        # Increase x-axis title size
    axis.text.y = element_text(size = 14),                         # Increase y-axis text size
    ## axis.title.y = element_text(size = 18),                         # Increase y-axis title size
    legend.position = "none"                                       # Remove legend
  ) +
  scale_y_continuous(breaks = seq(0.5, 1, by = 0.05), limits = c(0.5, 1.0))
  
# tiff("Figure S9.tiff", width = 12, height = 8, units = "in", res = 1000)
# ML_Methods_MCC
# dev.off()

png("Figure S9.png", width = 10, height = 7, units = "in", res = 1000)
ML_Methods_MCC
dev.off()

##############################################
############ Supplementary Table 4 ###########
##############################################

model_performances <- readRDS('predictor/performance_all_models_all_methods.rds')
metrics <- model_performances

metrics <- as.data.frame(metrics)
metrics[,-1] <- sapply(metrics[,-1], function(x){round(x,3)})
metrics <- as.data.table(metrics)

benchmarking_methods <- c('CADD', 'ESM-1b', 'AlphaMissense', 'VariPred', 'MutPred2', 'PolyPhen-2', 'PHACT', 'PHACTboost')
metrics <- metrics[!Predictor %in% benchmarking_methods,]

# Create a visually appealing table with flextable
flextable_object <- flextable(metrics) %>% 
  autofit() %>% 
  theme_booktabs() %>% 
  bold(part = "header") %>%
  align(align = "center", part = "all") 

flextable_object <- add_header_row(flextable_object, values = "Supplementary Table 4: Performance of 18 ML Methods evaluated using systematic feature ablation", colwidths = 8)
flextable_object <- hline_top(flextable_object, border = fp_border_default(width = 0), part = "header")

# Left align the "Methods" header
flextable_object <- align(flextable_object, i = 1, j = NULL, align = "left", part = "header")
# Make the "Methods" header bold
flextable_object <- bold(flextable_object, bold = TRUE, part = "header")
flextable_object <- font(flextable_object, fontname = "Cambria", part = "all")
flextable_object <- fontsize(flextable_object, size = 11, part = "all")
flextable_object <- fontsize(flextable_object, size = 12, part = "header")
flextable_object <- fontsize(flextable_object, size = 10, part = "footer")

for(i in 2:ncol(metrics)){
  
  metric_column <- as.numeric(unlist(metrics[,..i]))
  
  if(colnames(metrics)[i] == 'Precision'){
    toBold <- which(metric_column == min(metric_column, na.rm = TRUE))
  }
  else if(colnames(metrics)[i] == 'Recall'){
    toBold <- which(metric_column == max(metric_column, na.rm = TRUE))
  }
  else if(colnames(metrics)[i] == 'Accuracy'){
    toBold <- which(metric_column == max(metric_column, na.rm = TRUE))
  }
  else if(colnames(metrics)[i] == 'False_Positive_Rate'){
    toBold <- which(metric_column == min(metric_column, na.rm = TRUE))
  }
  else if(colnames(metrics)[i] == 'False_Negative_Rate'){
    toBold <- which(metric_column == min(metric_column, na.rm = TRUE))
  }
  else if(colnames(metrics)[i] == 'F1_Score'){
    toBold <- which(metric_column == max(metric_column, na.rm = TRUE))
  }
  else if(colnames(metrics)[i] == 'MCC'){
    toBold <- which(metric_column == max(metric_column, na.rm = TRUE))
  }
  flextable_object %<>% flextable::bold(toBold, i, bold = T, part = "body")
}

# save_as_image(flextable_object, path = "Table S3.png", width = 16, height = 96, units = "in", res = 600)
save_as_html(flextable_object, path = "Table S4.html")
write.xlsx(metrics, file = 'Table S4.xlsx', asTable = TRUE, overwrite = TRUE)

##############################################
############ Supplementary Table 5 ###########
##############################################

model_performances <- readRDS('predictor/performance_all_models_all_methods_wDataOverlap.rds')

metrics <- model_performances
selected_methods <- c('CADD', 'ESM-1b', 'AlphaMissense', 'VariPred', 'MutPred2', 'PolyPhen-2', 'PHACT', 'PHACTboost', 'naive_bayes~Inclusive')
metrics <- metrics[Predictor %in% selected_methods,]
metrics$Predictor[9] <- 'TrioNsight'

metrics <- as.data.frame(metrics)
metrics[,-1] <- sapply(metrics[,-1], function(x){round(x,3)})
metrics <- as.data.table(metrics)

# Create a visually appealing table with flextable
flextable_object <- flextable(metrics) %>% 
  autofit() %>% 
  theme_booktabs() %>% 
  bold(part = "header") %>%
  align(align = "center", part = "all") 

flextable_object <- add_header_row(flextable_object, values = "Supplementary Table 5: Benchmarking TrioNsight performance against variant pathogenicity prediction algorithms (Overlapping Training Dataset)", colwidths = 8)
flextable_object <- hline_top(flextable_object, border = fp_border_default(width = 0), part = "header")

# Left align the "Methods" header
flextable_object <- align(flextable_object, i = 1, j = NULL, align = "left", part = "header")
# Make the "Methods" header bold
flextable_object <- bold(flextable_object, bold = TRUE, part = "header")
flextable_object <- font(flextable_object, fontname = "Cambria", part = "all")
flextable_object <- fontsize(flextable_object, size = 11, part = "all")
flextable_object <- fontsize(flextable_object, size = 12, part = "header")
flextable_object <- fontsize(flextable_object, size = 10, part = "footer")

for(i in 2:ncol(metrics)){
  
  metric_column <- as.numeric(unlist(metrics[,..i]))
  
  if(colnames(metrics)[i] == 'Precision'){
    toBold <- which(metric_column == min(metric_column, na.rm = TRUE))
  }
  else if(colnames(metrics)[i] == 'Recall'){
    toBold <- which(metric_column == max(metric_column, na.rm = TRUE))
  }
  else if(colnames(metrics)[i] == 'Accuracy'){
    toBold <- which(metric_column == max(metric_column, na.rm = TRUE))
  }
  else if(colnames(metrics)[i] == 'False_Positive_Rate'){
    toBold <- which(metric_column == min(metric_column, na.rm = TRUE))
  }
  else if(colnames(metrics)[i] == 'False_Negative_Rate'){
    toBold <- which(metric_column == min(metric_column, na.rm = TRUE))
  }
  else if(colnames(metrics)[i] == 'F1_Score'){
    toBold <- which(metric_column == max(metric_column, na.rm = TRUE))
  }
  else if(colnames(metrics)[i] == 'MCC'){
    toBold <- which(metric_column == max(metric_column, na.rm = TRUE))
  }
  flextable_object %<>% flextable::bold(toBold, i, bold = T, part = "body")
}

# save_as_image(flextable_object, path = "Table S4.png", width = 16, height = 6, units = "in", res = 1000)
flextable::save_as_html(flextable_object, path = "Table S5.html")
write.xlsx(metrics, file = 'Table S5.xlsx', asTable = TRUE, overwrite = TRUE)