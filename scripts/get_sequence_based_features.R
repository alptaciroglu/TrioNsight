
#' Function to convert individual amino acid codes
#'
#' @param aa single amino acid in character format
#' @param aa_codes Vector to convert 1-letter amino acid codes to 3-letter amino acid codes
#' @param reverse_aa_codes Vector to convert 3-letter amino acid codes to 1-letter amino acid codes
#' @return list(one, three): One-letter and three-letter converted amino acid codes in list format.
#' @export
#' 
check_aa <- function(aa, aa_codes, reverse_aa_codes) {
  if (nchar(aa) == 1 && aa %in% names(aa_codes)) {
    return(list(one = aa, three = aa_codes[aa]))
  } else if (nchar(aa) == 3 && aa %in% aa_codes) {
    return(list(one = reverse_aa_codes[aa], three = aa))
  } else {
    return(list(one = NA, three = NA))  # Return NA for invalid codes
  }
}

#' Function to detect and convert amino acid codes
#'
#' @param variant_data Variant data in data.frame format
#' @param aa_codes Vector to convert 1-letter amino acid codes to 3-letter amino acid codes
#' @param reverse_aa_codes Vector to convert 3-letter amino acid codes to 1-letter amino acid codes
#' @param column_name Column name as string, matching the column in variant_data object that holds amino acid codes.
#' @return variant_data: with new columns for one-letter and three-letter amino acid codes
#' @export
convert_aa_codes <- function(variant_data, aa_codes, reverse_aa_codes, column_name) {
  
  variant_data[[column_name]] <- toupper(variant_data[[column_name]])
  aa_codes <- toupper(aa_codes)
  names(reverse_aa_codes) <- toupper(names(reverse_aa_codes))
  # Convert codes and add new columns
  converted <- lapply(variant_data[[column_name]], function(x){check_aa(x, aa_codes, reverse_aa_codes)})
  variant_data[[paste0(column_name, "_1letter")]] <- sapply(converted, `[[`, "one")
  variant_data[[paste0(column_name, "_3letter")]] <- sapply(converted, `[[`, "three")
  variant_data[[paste0(column_name, "_3letter")]] <- stringr::str_to_title(variant_data[[paste0(column_name, "_3letter")]])
  
  # Count conversions
  one_letter_count <- sum(nchar(variant_data[[column_name]]) == 1 & !is.na(variant_data[[paste0(column_name, "_3letter")]]))
  three_letter_count <- sum(nchar(variant_data[[column_name]]) == 3 & !is.na(variant_data[[paste0(column_name, "_1letter")]]))
  invalid_count <- sum(is.na(variant_data[[paste0(column_name, "_1letter")]]) & 
                         is.na(variant_data[[paste0(column_name, "_3letter")]]))
  
  return(variant_data)
}

#' Function to parse variant data; 1-check for format, 2-check for correct column names, 
#' 3-check for missense variants, and 4-convert amino acid codes
#'
#' @param variant_data Variant data in data.frame format
#' @param aa_codes Vector to convert 1-letter amino acid codes to 3-letter amino acid codes
#' @param reverse_aa_codes Vector to convert 3-letter amino acid codes to 1-letter amino acid codes
#' @return variant_data: controlled for format and missense variants
#' @export
#' 
parse_variant_data <- function(variant_data, aa_codes, reverse_aa_codes) {
  
  if(!is.data.frame(variant_data)) {
    stop(paste(' ', error_sign, " Please provide the variant data in data.frame format", sep = ''))
  }
  tryCatch(
    {
      as.character(unlist(variant_data['REF']))
      as.character(unlist(variant_data['ALT']))
    },
    error = function(e) {
      stop(paste(' ', error_sign, " Please provide the variant data in data.frame format with columns 'REF' and 'ALT'", sep = ''))
    }
  )
  
  variant_data <- convert_aa_codes(variant_data, aa_codes, reverse_aa_codes, column_name = "REF")
  variant_data <- convert_aa_codes(variant_data, aa_codes, reverse_aa_codes, column_name = "ALT")
  
  native_oneLetter <- as.character(unlist(variant_data['REF_1letter']))
  mutant_oneLetter <- as.character(unlist(variant_data['ALT_1letter']))
  native_threeLetter <- as.character(unlist(variant_data['REF_3letter']))
  mutant_threeLetter <- as.character(unlist(variant_data['ALT_3letter']))
  
  natInd <- which(nchar(native_oneLetter) > 0)
  mutInd <- which(nchar(mutant_oneLetter) > 0)
  missenseIdn <- intersect(natInd, mutInd)
  missenseIdn <- unique(c(missenseIdn, which(native_oneLetter == mutant_oneLetter)))
  
  if(length(missenseIdn) == 0) {
    stop(paste(' ', error_sign, ' Please provide the variant data containing missense mutations with correct amino acid notations!', sep = ''))
  }
  else if(length(missenseIdn) < nrow(variant_data)) {
    variant_data <- variant_data[missenseIdn,]
    cat(paste(' ', warning_sign, ' Some variants are not missense => ', nrow(variant_data)-length(missenseIdn), '/', nrow(variant_data), ' ', warning_sign, '\n', sep = ''))
    cat(paste(' ', warning_sign, ' Removing non-missense variants from the data as only missense variants are supported ', warning_sign, '\n', sep = ''))
  }
  
  cat(paste(' Parsed ', nrow(variant_data), ' variants for missense variants and single character amino acids ', completion_sign, '\n', sep = ''))
  return(variant_data)
}

#' Calculate ScoreCons for given multiple sequence alignment
#'
#' @param fasta_file fasta file containing the protein sequences
#' @param variant_data Variant data in data.frame format
#' @return Write and return file containing ScoreCons data
#' @export
#' 
run_scorecons <- function(fasta_file, variant_data) {
  # Extract the base name of the file (remove .fasta extension)
  base_name <- tools::file_path_sans_ext(fasta_file)
  
  cat(paste(' ', calculation_sign, ' Running ScoreCons for: ', fasta_file, ' ', calculation_sign, '\n', sep = ''))
  # Construct the command
  command <- sprintf("./scorecons -a %s -m resources/PET91mod.mat2 -o %s.scorecons", 
                     fasta_file, base_name)
  
  # Run the command
  dops <- system(command, intern = TRUE)
  dops <- as.numeric(gsub("DOPS score: ", "", dops))
  if(dops < 70 & dops > 50){
    cat(paste(' ', warning_sign, ' DOPS score is below 70, which indicates that the alignment is not informative ', warning_sign, '\n', sep = ''))
    cat('\n\n')
  } else if(dops <= 50){
    stop(paste(' ', error_sign, ' DOPS score is below 50, which indicates that the alignment is very uninformative!', sep = ''))
  }
  if(dops >= 70){
    cat(paste(' DOPS score is above 70, which indicates that the alignment is informative. ', completion_sign, '\n', sep = ''))
  }
  
  # Clean up the scorecons output file and write the cleaned content to a new file
  path_to_scorecons_file <- paste0(base_name, ".scorecons")  
  parse_raw_scorecons(file_path = path_to_scorecons_file, fasta_file = fasta_file)
  cat(paste(' ', info_sign, ' ScoreCons output file: ', path_to_scorecons_file, ' ', info_sign, '\n', sep = ''))
  
  # Read scorecons output and add to variant data
  cleaned_content <- readLines(paste0(base_name, ".scorecons"))
  cleaned_content <- trimws(sapply(strsplit(cleaned_content, '#', fixed = T), '[', 1))
  cleaned_content <- as.numeric(cleaned_content)
  scorecons_scores <- cleaned_content[variant_data$MSALocation]
  variant_data['ScoreCons'] <- scorecons_scores

  return(variant_data)
}

#' Clean up the ScoreCons output file to make it ready for further processing. Write the cleaned content to a new file
#'
#' @param file_path ScoreCons output file obtained by running run_scorecons function or supplied externally
#' @param fasta_file fasta file containing the protein sequences matching the ScoreCons output dimensions
#' 
parse_raw_scorecons <- function(file_path, fasta_file) {
  # Read the content of the file
  content <- readLines(file_path, warn = FALSE)
  fasta_file <- Biostrings::readAAStringSet(fasta_file)

  # Remove empty lines
  content <- content[nzchar(content)]
  
  # Remove lines that were split incorrectly
  combined_content <- content[which(nchar(content) == median(nchar(content)))]
  
  check_content <- sapply(strsplit(combined_content, '#              ', fixed = T), '[', 2)
  check_content <- lapply(check_content, function(x) {strsplit(x, '', fixed = T)[[1]]})
  if(length(unlist(unique(lapply(check_content, length)))) == 1){
    for(i in 1:length(check_content)){
      scoreConsSeq <- paste(unlist(lapply(check_content, function(x){x[i]})), collapse = '')
      fastaSeq <- as.character(fasta_file[[i]])
      if(scoreConsSeq != fastaSeq){
        stop(paste(' ', error_sign, ' The ScoreCons output file is not formatted correctly. Please check the file and try again!', sep = ''))
      }
    }
  }
  else{
    stop(paste(' ', error_sign, ' The ScoreCons output file is not formatted correctly. Please check the file and try again!', sep = ''))
  }
  
  # Write the cleaned content back to the file
  writeLines(combined_content, file_path)
}

#' Reads the ScoreCons output file and annotates the variant data with the ScoreCons scores
#'
#' @param ScoreConsFile ScoreCons output file obtained by running run_scorecons function or supplied externally
#' @param variant_data Variant data containing variants, and some already added features in data.frame format
#' @return return variant data annotated with ScoreCons scores
#' @export
#' 
use_precalculated_scorecons <- function(ScoreConsFile, variant_data){
  
  # Read scorecons output and add to variant data
  cleaned_content <- readLines(ScoreConsFile)
  cleaned_content <- trimws(sapply(strsplit(cleaned_content, '#', fixed = T), '[', 1))
  cleaned_content <- as.numeric(cleaned_content)
  scorecons_scores <- cleaned_content[variant_data$MSALocation]
  variant_data['ScoreCons'] <- scorecons_scores
  
  return(variant_data)
}

#' Calculate BLOSUM score for given variant data
#'
#' @param variant_data Variant data in data.frame format
#' @param subsMatrix Substitution matrix to be used. Default is BLOSUM80. Either obtain from ftp or used from internal BLOSUM80 matrix
#' @return variant_data annotated with BLOSUM score
#' @export
#' 
get_blosum_score <- function(variant_data, subsMatrix = 'BLOSUM80') {
  
  blosumChoice <- subsMatrix
  
  if(subsMatrix == "BLOSUM80"){
    tryCatch(
      {
        file <- "ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM80"
        subsMatrix <- as.matrix(read.table(file, check.names=FALSE))
      },
      error = function(e) {
        subsMatrix <- readRDS('resources/internal_blosum80.rds')
      }
    )
  } else if(subsMatrix == "BLOSUM62"){
    tryCatch(
      {
        file <- "ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62"
        subsMatrix <- as.matrix(read.table(file, check.names=FALSE))
      },
      error = function(e) {
        subsMatrix <- readRDS('resources/internal_blosum62.rds')
      }
    )
  }
  
  variant_data_split <- split(variant_data[,c('REF_1letter', 'ALT_1letter')], seq(nrow(variant_data[,c('REF_1letter', 'ALT_1letter')])))
  BLOSUM_scores <- lapply(variant_data_split, function(x){subsMatrix[x$REF_1letter, x$ALT_1letter]})
  BLOSUM_scores <- as.numeric(unlist(BLOSUM_scores))
  cat(paste(' Annotated ', length(BLOSUM_scores), ' variants with ', blosumChoice, ' scores. ', completion_sign, '\n', sep = ''))
  
  variant_data[blosumChoice] <- BLOSUM_scores
  return(variant_data)
}

#' Calculate Grantham distance for given variant data
#'
#' @param variant_data Variant data in data.frame format
#' @return variant_data annotated with grantham distance
#' @export
#' 
get_grantham_distance <- function(variant_data) {
  
  grantham_scores <- grantham::grantham_distance(variant_data$REF_3letter, variant_data$ALT_3letter)
  grantham_scores <- as.data.frame(grantham_scores)
  grantham_scores <- grantham_scores$d
  cat(paste(' Annotated ', length(grantham_scores), ' variants with Grantham scores. ', completion_sign, '\n', sep = ''))
  
  variant_data['Grantham'] <- grantham_scores
  variant_data['REF_1letter'] <- NULL
  variant_data['ALT_1letter'] <- NULL
  variant_data['REF_3letter'] <- NULL
  variant_data['ALT_3letter'] <- NULL
  return(variant_data)
}

#' A wrapper function to calculate sequence-based features. Calls other functions to parse and check format of variant data object, calculate ScoreCons, BLOSUM and Grantham scores
#'
#' @param variants_inherited variant data containing variants in data.frame format, output of variant inheritance process obtained from inherit_variants_from_msa function
#' @param ScoreConsFile ScoreCons output file obtained by running run_scorecons function or supplied externally
#' @param alignmentFile alignment file in fasta format containing sequences, with ids matching the variant data
#' @param usePrecalculatedScoreCons logical value indicating whether to use precalculated ScoreCons scores or not
#' @param BLOSUM_Matrix Substitution matrix to be used. Default is BLOSUM80. Either obtain from ftp or used from internal BLOSUM80 matrix
#' @return return variant data annotated with ScoreCons scores
#' @export
#' 
get_sequence_based_features <- function(variants_inherited, ScoreConsFile, alignmentFile, usePrecalculatedScoreCons, BLOSUM_Matrix = 'BLOSUM80'){
  
  if(is.null(ScoreConsFile) & isTRUE(usePrecalculatedScoreCons)){
    stop(paste(' ', error_sign, ' Please either provide a ScoreCons file or allow a new ScoreCons run', sep = ''))
  }
  
  underscore_start("Calculating sequence-based features")
  
  variant_data_annotated <- parse_variant_data(variant_data = variants_inherited, aa_codes, reverse_aa_codes)
  
  if(isFALSE(usePrecalculatedScoreCons)){
    variant_data_annotated <- run_scorecons(fasta_file = alignmentFile, 
                                            variant_data = variant_data_annotated)
  }
  else{
    variant_data_annotated <- use_precalculated_scorecons(ScoreConsFile = ScoreConsFile, variant_data = variant_data_annotated)
    cat(paste(' ', info_sign, ' Using precalculated ScoreCons scores ', info_sign, '\n', sep = ''))
  }

  variant_data_annotated <- get_blosum_score(variant_data = variant_data_annotated, subsMatrix = BLOSUM_Matrix)
  variant_data_annotated <- get_grantham_distance(variant_data = variant_data_annotated)
  variant_data_annotated <- as.data.table(variant_data_annotated)
  variant_data_annotated <- convert_aa_codes(variant_data_annotated, aa_codes, reverse_aa_codes, column_name = "REF")
  variant_data_annotated <- convert_aa_codes(variant_data_annotated, aa_codes, reverse_aa_codes, column_name = "ALT")
  variant_data_annotated[, Accession := sapply(strsplit(Variant_Seq, "_", fixed = T), '[', 1)]
  variant_data_annotated[, target_id := paste(Accession, '_', REF, Variant_Location, ALT, sep = '')]
  variant_data_annotated[, protVar_identifiers := paste(Accession, ' ', REF_3letter, Variant_Location, ALT_3letter, sep = '')]
  variant_data_annotated[, Accession := NULL]
  
  underscore_finish("Sequence-based features obtained successfully")
  
  return(variant_data_annotated)
}
