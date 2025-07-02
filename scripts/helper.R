####################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages('BiocManager')   
}
library(BiocManager)
if (!requireNamespace("data.table", quietly = TRUE)) {
  BiocManager::install('data.table')
}
library(data.table)
if (!requireNamespace("httr", quietly = TRUE)) {
  BiocManager::install('httr')
}
library(httr)
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  BiocManager::install('jsonlite')
}
library(jsonlite)
if (!requireNamespace("rjson", quietly = TRUE)) {
  BiocManager::install('rjson')
}
library(rjson)
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install('Biostrings')
}
library(Biostrings)
if (!requireNamespace("stringr", quietly = TRUE)) {
  install.packages('stringr')   
}
library(stringr)
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)
if (!requireNamespace("stats", quietly = TRUE)) {
  install.packages("stats")
}
library(stats)
if (!requireNamespace("bio3d", quietly = TRUE)) {
  install.packages("bio3d")
}
library(bio3d)
if (!requireNamespace("grantham", quietly = TRUE)) {
  BiocManager::install('grantham')
}
library(grantham)
if (!requireNamespace("adabag", quietly = TRUE)) {
  BiocManager::install('adabag')
}
library(adabag)
if (!requireNamespace("caret", quietly = TRUE)) {
  BiocManager::install('caret')
}
library(caret)
if (!requireNamespace("flextable", quietly = TRUE)) {
  BiocManager::install('flextable')
}
library(flextable)
if (!requireNamespace("magrittr", quietly = TRUE)) {
  BiocManager::install('magrittr')
}
library(magrittr)
if (!requireNamespace("officer", quietly = TRUE)) {
  BiocManager::install('officer')
}
library(officer)

source('scripts/variant_inheritance.R')
source('scripts/get_sequence_based_features.R')
source('scripts/connect_protVar_API.R')
source('scripts/VariPred.R')
source('scripts/SAAPdap.R')
source('scripts/distance_to_variant_clusters.R')

####################################################

aa_codes <- c(
  "A" = "Ala", "R" = "Arg", "N" = "Asn", "D" = "Asp", "C" = "Cys",
  "E" = "Glu", "Q" = "Gln", "G" = "Gly", "H" = "His", "I" = "Ile",
  "L" = "Leu", "K" = "Lys", "M" = "Met", "F" = "Phe", "P" = "Pro",
  "S" = "Ser", "T" = "Thr", "W" = "Trp", "Y" = "Tyr", "V" = "Val"
)
reverse_aa_codes <- setNames(names(aa_codes), aa_codes)
AAs <- readRDS('resources/AAs_cheatSheet.rds')
Prots2PDBsAnnot <- readRDS('resources/Prots2PDBsAnnot.rds')
PHACT_TrioN_like <- readRDS('stable_run/PHACT_TrioN_like.rds')
DHDomains2PDBsAnnot <- readRDS('resources/DHDomains2PDBsAnnot.rds')
reference_clustering <- readRDS('stable_run/reference_clustering.rds')
model_domainboundaries <- readRDS('resources/model_domainboundaries.rds')
trainingSet_Overlaps <- readRDS('resources/trainingSet_Overlaps.rds')

####################################################

Sys.setlocale("LC_TIME", "en_US.utf8")
todaysDate <- format(Sys.Date(), format = "%d%m%y")

####################################################

track_progress <- function(i, total_rows, bar_length = 50) {
  progress <- i / total_rows
  completed_length <- round(progress * bar_length)
  remaining_length <- bar_length - completed_length
  
  cat("Working [", paste(rep("=", completed_length), collapse = ""), ">",
      paste(rep("-", remaining_length), collapse = ""), "] ",
      round(progress * 100), "%\r", sep = "")
  
  if (i == total_rows) {
    cat('\n -> Done\n\n')
  }
  
  # Flush output to console to show progress immediately
  flush.console()
  
}

####################################################

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

####################################################

table2dataframe <- function(table){
  table_df <- as.data.frame.matrix(table)
  table_df <- as.data.frame(cbind(rownames(table_df), table_df))
  colnames(table_df)[1] <- 'Genes'
  rownames(table_df) <- NULL
  return(table_df)
}

####################################################

extractNumbersFromString <- function(my_strings){
  
  # Use gsub to remove all non-numeric characters
  my_numbers <- gsub("[^[:digit:].-]", "", my_strings)
  
  # Convert the resulting vector of strings to numeric
  my_numbers <- as.numeric(my_numbers)
  
  # View the resulting vector of numbers
  return(my_numbers)
  
}

####################################################

extractLinesWithKeyword <- function(file_path, keyword) {
  # Read lines from the file
  lines <- readLines(file_path)
  
  # Use grep to find lines starting with the keyword
  keyword_lines <- grep(paste0("^", keyword), lines, value = TRUE)
  
  return(t(as.data.table(strsplit(keyword_lines, '  ', fixed = T)[[1]])))
}

####################################################

# Function to calculate precision
precision <- function(predicted, actual, positive_class) {
  true_positives <- sum(predicted == positive_class & predicted == actual)
  false_positives <- sum(predicted == positive_class & predicted != actual)
  
  if(true_positives == 0) {
    return(0)
  } else if (true_positives == 0 & false_positives == 0) {
    return(0)
  } else {
    return(round((true_positives / (true_positives + false_positives)), 4))
  }
}

# Function to calculate recall
recall <- function(predicted, actual, positive_class) {
  true_positives <- sum(predicted == positive_class & predicted == actual)
  false_negatives <- sum(predicted != positive_class & predicted != actual)
  
  if(true_positives == 0) {
    return(0)
  } else if (true_positives == 0 & false_negatives == 0) {
    return(0)
  } else {
    return(round((true_positives / (true_positives + false_negatives)), 4))
  }
}

# Function to calculate accuracy
accuracy <- function(predicted, actual) {
  correct_predictions <- sum(predicted == actual)
  total_predictions <- length(predicted)
  
  return(round((correct_predictions / total_predictions), 4))
}

# Function to calculate false positive rate
false_positive_rate <- function(predicted, actual, positive_class) {
  false_positives <- sum(predicted == positive_class & predicted != actual)
  true_negatives <- sum(predicted != positive_class & predicted == actual)
  
  if(false_positives == 0) {
    return(0)
  } else if (false_positives == 0 & true_negatives == 0) {
    return(0)
  } else {
    return(round((false_positives / (false_positives + true_negatives)), 4))
  }
}

# Function to calculate false negative rate
false_negative_rate <- function(predicted, actual, positive_class) {
  false_negatives <- sum(predicted != positive_class & predicted != actual)
  true_positives <- sum(predicted == positive_class & predicted == actual)
  
  if(false_negatives == 0) {
    return(0)
  } else if (false_negatives == 0 & true_positives == 0) {
    return(0)
  } else {
    return(round((false_negatives / (false_negatives + true_positives)), 4))
  }
}

# Function to calculate F1 score
f1_score <- function(predicted, actual, positive_class) {
  prec <- precision(predicted, actual, positive_class)
  rec <- recall(predicted, actual, positive_class)
  
  if((prec + rec) == 0) {
    return(0)
  } else {
    return(round((2 * (prec * rec) / (prec + rec)), 4))
  }
}

mcc_score <- function(predicted, actual, positive_class) {
  true_positives <- sum(predicted == positive_class & predicted == actual)
  true_negatives <- sum(predicted != positive_class & predicted == actual)
  false_positives <- sum(predicted == positive_class & predicted != actual)
  false_negatives <- sum(predicted != positive_class & predicted != actual)
  
  numerator <- (true_positives * true_negatives) - (false_positives * false_negatives)
  denominator <- sqrt(true_positives + false_positives) * sqrt(true_positives + false_negatives) * sqrt(true_negatives + false_positives) * sqrt(true_negatives + false_negatives)
  
  if(denominator == 0) {
    return(0)
  } else {
    return(round(numerator / denominator, 4))
  }
}

####################################################

warning_sign <- "\u26A0"
calculation_sign <- "\u23F3"
completion_sign <- "\u2714"
error_sign <- "\U274C"
info_sign <- "\U24D8"

underscore_start <- function(message = "Start", length = 40) {
  cat(strrep("\u2500", length), "\n")
  cat(message, "\n")
  cat(strrep("\u2500", length), "\n")
}

underscore_finish <- function(message = "Finish", length = 40) {
  cat(strrep("\u2501", length), "\n")
  cat(message, "\n")
  cat(strrep("\u2501", length), "\n")
}

print_underlined <- function(text, underline_char = "\u2500", padding = 2) {
  # Print the text with padding
  cat(paste0(rep(" ", padding), collapse = ""), text, "\n")
  
  # Calculate the length of the text plus padding
  total_length <- nchar(text) + (2 * padding)
  
  # Print the underline
  cat(strrep(underline_char, total_length), "\n")
}

####################################################

# Function to generate all possible single amino acid substitutions from data.frame
generate_substitutions_from_df <- function(df, position_field, residue_field) {
  aa_list <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  
  residues <- df[[position_field]]
  original_aa <- df[[residue_field]]
  
  substitutions <- data.frame(
    Residue = rep(residues, each = length(aa_list)),
    original = rep(original_aa, each = length(aa_list)),
    substitution = rep(aa_list, times = length(residues))
  )
  
  substitutions <- substitutions[substitutions$original != substitutions$substitution, ]
  substitutions <- as.data.table(substitutions)
  
  return(substitutions)
}
