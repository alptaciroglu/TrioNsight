#' Retrieve protein variant mappings from the ProtVar API
#'
#' @param variants4protVar in "CHROM POS REF ALT" or "UniProtID NativeAAResidueNoMutantAA" format
#' @param include_function logical indicating whether to include functional annotations
#' @param include_population logical indicating whether to include population annotations
#' @param include_structure logical indicating whether to include structural annotations
#' @return results from the ProtVar API for the given variants in JSON format
#' @export
#' 
get_protVar_features <- function(variants4protVar, include_function = TRUE, include_population = TRUE, include_structure = TRUE) {
  
  # API endpoint
  url <- "https://www.ebi.ac.uk/ProtVar/api/mappings"
  
  # Query parameters
  params <- list(
    `function` = tolower(as.character(include_function)),
    population = tolower(as.character(include_population)),
    structure = tolower(as.character(include_structure))
  )
  
  # Request body
  body <- jsonlite::toJSON(variants4protVar, auto_unbox = TRUE)
  
  # Make the POST request
  response <- httr::POST(
    url = url,
    query = params,
    body = body,
    httr::add_headers(
      "accept" = "application/json",
      "Content-Type" = "application/json"
    )
  )
  
  # Check if the request was successful
  if (httr::status_code(response) == 200) {
    # Parse the JSON response
    result <- httr::content(response, "parsed")
    return(result)
  } else {
    stop(paste("Error:", httr::status_code(response), "\nMessage:", httr::content(response, "text")))
  }
}

#' Retrieve protein variant mappings from the ProtVar API and parse the results
#'
#' @param protVar_results results from the ProtVar API in JSON format
#' @param variantFormat character string indicating the format of the variants. protVar accepts either genomic coordinates in "CHROM POS REF ALT" format or UniProt identfiers in "UniProtID NativeAAResidueNoMutantAA" format
#' @return parsed results from the ProtVar API in data.frame format matching the variant data
#' @export
#' 
parse_protVar_features <- function(protVar_results, variantFormat){
  
  # Define the paths for each field
  paths <- list(
    Query = c('inputStr'),
    varSite = c("mappings", 1, "genes", 1, "isoforms", 1, "conservScore", "score"),
    AlphaMissense_Score = c("mappings", 1, "genes", 1, "isoforms", 1, "amScore", "amPathogenicity"),
    AlphaMissense_Class = c("mappings", 1, "genes", 1, "isoforms", 1, "amScore", "amClass"),
    ESM_Score = c("mappings", 1, "genes", 1, "isoforms", 1, "esmScore", "score"),
    foldxDdg = c("mappings", 1, "genes", 1, "isoforms", 1, "referenceFunction", "foldxs", 1, "foldxDdg"),
    foldx_plddt = c("mappings", 1, "genes", 1, "isoforms", 1, "referenceFunction", "foldxs", 1, "plddt"),
    ext_pred1_name = c("mappings", 1, "genes", 1, "isoforms", 1, "populationObservations", "proteinColocatedVariant", 1, "predictions", 1, "predAlgorithmNameType"),
    ext_pred1_score = c("mappings", 1, "genes", 1, "isoforms", 1, "populationObservations", "proteinColocatedVariant", 1, "predictions", 1, "score"),
    ext_pred1_type = c("mappings", 1, "genes", 1, "isoforms", 1, "populationObservations", "proteinColocatedVariant", 1, "predictions", 1, "predictionValType"),
    ext_pred2_name = c("mappings", 1, "genes", 1, "isoforms", 1, "populationObservations", "proteinColocatedVariant", 1, "predictions", 2, "predAlgorithmNameType"),
    ext_pred2_score = c("mappings", 1, "genes", 1, "isoforms", 1, "populationObservations", "proteinColocatedVariant", 1, "predictions", 2, "score"),
    ext_pred2_type = c("mappings", 1, "genes", 1, "isoforms", 1, "populationObservations", "proteinColocatedVariant", 1, "predictions", 2, "predictionValType")
  )
  
  # Define the paths for each field
  paths_uniprot <- list(
    protVar_identifiers = c('inputStr'),
    varSite = c("derivedGenomicInputs", 1, "mappings", 1, "genes", 1, "isoforms", 1, "conservScore", "score"),
    AlphaMissense_Score = c("derivedGenomicInputs", 1, "mappings", 1, "genes", 1, "isoforms", 1, "amScore", "amPathogenicity"),
    AlphaMissense_Class = c("derivedGenomicInputs", 1, "mappings", 1, "genes", 1, "isoforms", 1, "amScore", "amClass"),
    ESM_Score = c("derivedGenomicInputs", 1, "mappings", 1, "genes", 1, "isoforms", 1, "esmScore", "score"),
    foldxDdg = c("derivedGenomicInputs", 1, "mappings", 1, "genes", 1, "isoforms", 1, "referenceFunction", "foldxs", 1, "foldxDdg"),
    foldx_plddt = c("derivedGenomicInputs", 1, "mappings", 1, "genes", 1, "isoforms", 1, "referenceFunction", "foldxs", 1, "plddt"),
    ext_pred1_name = c("derivedGenomicInputs", 1, "mappings", 1, "genes", 1, "isoforms", 1, "populationObservations", "proteinColocatedVariant", 1, "predictions", 1, "predAlgorithmNameType"),
    ext_pred1_score = c("derivedGenomicInputs", 1, "mappings", 1, "genes", 1, "isoforms", 1, "populationObservations", "proteinColocatedVariant", 1, "predictions", 1, "score"),
    ext_pred1_type = c("derivedGenomicInputs", 1, "mappings", 1, "genes", 1, "isoforms", 1, "populationObservations", "proteinColocatedVariant", 1, "predictions", 1, "predictionValType"),
    ext_pred2_name = c("derivedGenomicInputs", 1, "mappings", 1, "genes", 1, "isoforms", 1, "populationObservations", "proteinColocatedVariant", 1, "predictions", 2, "predAlgorithmNameType"),
    ext_pred2_score = c("derivedGenomicInputs", 1, "mappings", 1, "genes", 1, "isoforms", 1, "populationObservations", "proteinColocatedVariant", 1, "predictions", 2, "score"),
    ext_pred2_type = c("derivedGenomicInputs", 1, "mappings", 1, "genes", 1, "isoforms", 1, "populationObservations", "proteinColocatedVariant", 1, "predictions", 2, "predictionValType")
  )
  
  # Define a helper function to safely extract nested values
  safe_extract <- function(x, path) {
    
    for (i in seq_along(path)) {
      p <- path[i]
      save_path <- p
      
      tryCatch({
        p <- as.numeric(p)
      },
      error = function(e) {
        stop(paste(error_sign, ' protVar returned an unexpected error', sep = ''))
      },
      warning = function(w){
        p <- save_path
      })
      
      tryCatch({
        if (i == length(path)) {
          return(if (!is.null(x[[p]])) x[[p]] else NA)
        }
        if (is.list(x) && length(x) == 0) {
          return(NA)
        }

        if (is.list(x) && !is.null(x[[p]])) {
          x <- x[[p]]
        } else if (is.numeric(p) && length(x) >= p) {
          x <- x[[p]]
        } else {
          return(NA)
        }
      },
      error = function(e) {
        return(NA)
      })
    }
    return(NA)
  }

  # Extract all fields at once
  if(variantFormat == 'UniProt') {
    paths <- paths_uniprot
  }
  protVar_results_parsed <- lapply(paths, function(path) {
    sapply(protVar_results$inputs, function(x) safe_extract(x, path))
  })
  
  # Convert the list to a data table
  protVar_results_parsed <- as.data.table(protVar_results_parsed)
  protVar_results_parsed[, PolyPhen_Score := 0]
  protVar_results_parsed[, PolyPhen_Class := '']
  protVar_results_parsed[, SIFT_Score := 0]
  protVar_results_parsed[, SIFT_Class := '']
  
  protVar_results_parsed[ext_pred1_name == 'PolyPhen']$PolyPhen_Score <- 
    protVar_results_parsed[ext_pred1_name == 'PolyPhen']$ext_pred1_score
  protVar_results_parsed[ext_pred1_name == 'PolyPhen']$PolyPhen_Class <-
    protVar_results_parsed[ext_pred1_name == 'PolyPhen']$ext_pred1_type
  protVar_results_parsed[ext_pred2_name == 'PolyPhen']$PolyPhen_Score <-
    protVar_results_parsed[ext_pred2_name == 'PolyPhen']$ext_pred2_score
  protVar_results_parsed[ext_pred2_name == 'PolyPhen']$PolyPhen_Class <-
    protVar_results_parsed[ext_pred2_name == 'PolyPhen']$ext_pred2_type
  protVar_results_parsed[ext_pred1_name == 'SIFT']$SIFT_Score <-
    protVar_results_parsed[ext_pred1_name == 'SIFT']$ext_pred1_score
  protVar_results_parsed[ext_pred1_name == 'SIFT']$SIFT_Class <-
    protVar_results_parsed[ext_pred1_name == 'SIFT']$ext_pred1_type
  protVar_results_parsed[ext_pred2_name == 'SIFT']$SIFT_Score <-
    protVar_results_parsed[ext_pred2_name == 'SIFT']$ext_pred2_score
  protVar_results_parsed[ext_pred2_name == 'SIFT']$SIFT_Class <-
    protVar_results_parsed[ext_pred2_name == 'SIFT']$ext_pred2_type
  protVar_results_parsed[PolyPhen_Score == 0,]$PolyPhen_Score <- NA
  protVar_results_parsed[PolyPhen_Class == '',]$PolyPhen_Class <- NA
  protVar_results_parsed[SIFT_Score == 0,]$SIFT_Score <- NA
  protVar_results_parsed[SIFT_Class == '',]$SIFT_Class <- NA

  protVar_results_parsed[, ext_pred1_name := NULL]
  protVar_results_parsed[, ext_pred1_score := NULL]
  protVar_results_parsed[, ext_pred1_type := NULL]
  protVar_results_parsed[, ext_pred2_name := NULL]
  protVar_results_parsed[, ext_pred2_score := NULL]
  protVar_results_parsed[, ext_pred2_type := NULL]
  
  protVar_results_parsed <- as.data.frame(protVar_results_parsed)
  protVar_results_parsed$AlphaMissense_Class <- str_to_title(protVar_results_parsed$AlphaMissense_Class)
  return(protVar_results_parsed)
}

#' Wrapper function to query the ProtVar API and parse the results
#'
#' @param variant_data data.frame containing the variant data
#' @param variants character vector of variants to query protVar for
#' @param varFormat character string specifying the format of the variants. Must be either 'UniProt' or 'Genomic'
#' @param outputFolder character string specifying the output folder to save the results. Default is NULL. If NULL, the results will saved to tmp_protVar_results folder in the working directory
#' @param chunk_size integer specifying the number of variants to query ProtVar at a time. Default is 100
#' @param max_attempts integer specifying the maximum number of attempts to query ProtVar for a chunk. Default is 10
#' @return variant_data data.frame containing the variant data with the ProtVar features appended
#' @export
#' 
protVar_features_wrapper <- function(variant_data, variants, varFormat = 'UniProt', outputFolder = NULL, chunk_size = 100, max_attempts = 10) {
  
  if(varFormat != 'UniProt' & varFormat != 'Genomic'){
    stop(paste(error_sign, ' varFormat must be either "UniProt" i.e. O43307 Trp192Gly or "Genomic" i.e. 13 113081271 G A', sep = ''))
  }
  
  underscore_start("ProtVar connection is being established")
  
  cat(paste(' ', info_sign, ' At the end of this process, the following features will be obtained ', info_sign, '\n', sep = ''))
  cat(paste(' ', info_sign, ' AlphaMissense Predictions - foldX stability calculations ', info_sign, '\n', sep = ''))
  cat(paste(' ', info_sign, ' ESM scores for deleteriousness - VarSite conservation scores ', info_sign, '\n', sep = ''))
  
  cat('\n')
  cat(paste(' ', warning_sign, ' ProtVar data will be saved to a temporary folder ', warning_sign, '\n', sep = ''))
  cat(paste(' ', warning_sign, ' Please do not modify the folder while the script is running ', warning_sign, '\n', sep = ''))
  cat(paste(' ', warning_sign, ' After the feature extraction is complete, this folder will be deleted ', warning_sign, '\n', sep = ''))
  cat('\n')
  
  # Function to process a single chunk
  process_chunk <- function(chunk) {
    attempts <- 0
    while (attempts < max_attempts) {
      attempts <- attempts + 1
      tryCatch({
        result <- get_protVar_features(variants4protVar = chunk)
        cat(paste(' Chunk processed successfully on attempt ', attempts, ' ', completion_sign, '\n', sep = ''))
        return(result)
      }, error = function(e) {
        if (attempts < max_attempts) {
          cat(paste(' ', calculation_sign, ' Attempt ', attempts, ' failed. Retrying ', calculation_sign, '\n', sep = ''))
        } else {
          cat(paste(' ', warning_sign, ' All ', max_attempts ,' attempts failed for this chunk. Moving on to the next chunk ', warning_sign, '\n', sep = ''))
        }
      })
    }
    return(NULL)  # Return NULL if all attempts fail
  }
  
  # Split variants into chunks
  chunks <- split(variants, ceiling(seq_along(variants) / chunk_size))
  
  if(is.null(outputFolder)){
    outputFolder <- 'tmp_protVar_results'
  }
  dir.create(outputFolder, showWarnings = FALSE)
  setwd(outputFolder)
  
  # Process each chunk
  failedVarCount <- 0
  for (i in seq_along(chunks)){
    cat(paste(' ', calculation_sign, ' Processing chunk ', i, ' of ', length(chunks), ' ', calculation_sign, '\n', sep = ''))
    chunk_result <- process_chunk(chunks[[i]])
    if (is.null(chunk_result)) {
      failedVarCount <- failedVarCount + length(chunks[[i]])
      next
    }
    saveRDS(chunk_result, file = paste('protVar_chunk_', i ,'.rds', sep = ''))
  }
  
  if(failedVarCount > 0){
    cat(paste(' ', warning_sign, ' ProtVar data could not be obtained for ', failedVarCount, '/', length(variants), ' variants ', warning_sign, '\n', sep = ''))
  }
  else{
    cat(paste(' ProtVar data successfully obtained for all ', length(variants), ' ', completion_sign, '\n', sep = ''))
  }
  
  read_files <- function(folder){
    tryCatch({
      chunk_files <- list.files(pattern = 'protVar_chunk_')
      if (length(chunk_files) == 0) {
        stop(paste(' ', error_sign, ' No protVar results could be found in the temporary folder!', sep = ''))
      }
      protVar_features_extracted <- lapply(chunk_files, function(f) parse_protVar_features(protVar_results = readRDS(f), variantFormat = varFormat))
      protVar_features_extracted <- as.data.table(do.call(rbind, protVar_features_extracted))
      return(protVar_features_extracted)
    },error = function(e){
      cat(paste(' ', warning_sign, ' An error occurred while reading the protVar data from temporary folder: ', conditionMessage(e), '  ', warning_sign, '\n', sep = ''))
      return(FALSE)
    }, finally = {
      # Always return to the parent directory, even if an error occurred
      setwd("..")
    })
  }
  
  protVar_features_extracted <- read_files(outputFolder)
  if(isFALSE(protVar_features_extracted)) {
    cat(paste(' ', warning_sign, ' Processing failed. Temporary folder ', outputFolder, ' is deleted. ', warning_sign, '\n', sep = ''))
    unlink(outputFolder, recursive = TRUE)
    cat(paste(' ', warning_sign, ' Try running the script at a later time. ', warning_sign, '\n', sep = ''))
  } else {
    # Delete the folder and its contents
    cat(paste(' Processing complete. Cleaning up. ', completion_sign, '\n', sep = ''))
    variant_data <- merge(variant_data, protVar_features_extracted, by = 'protVar_identifiers', all.x = TRUE)
    variant_data[, protVar_identifiers := NULL]
    unlink(outputFolder, recursive = TRUE)
    underscore_finish("ProtVar features obtained. Closing the connection")
    return(variant_data)
  }
}
