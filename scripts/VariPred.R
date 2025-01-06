#' Generate a file ready to be used in VariPred
#'
#' @param variant_data a data frame containing the variant data
#' @param fullLengthSequences a fasta file containing the full length protein sequences related to variants
#' @export
#' 
generateVariPredInput <- function(variant_data, fullLengthSequences){
  
  variant_data <- unique(variant_data[,.(target_id, Variant_Location, REF, ALT)])
  target_id <- variant_data$target_id
  aa_index <- variant_data$Variant_Location
  wt_aa <- variant_data$REF
  mt_aa <- variant_data$ALT
  
  accessions <- sapply(strsplit(target_id, '_', fixed = T), '[', 1)
  sequence <- Biostrings::readAAStringSet(fullLengthSequences)
  sequence <- as.character(as.character(sequence[accessions]))
  
  varipred_df <- data.table(target_id = target_id, aa_index = aa_index, wt_aa = wt_aa, mt_aa = mt_aa, wt_seq = sequence)

  varipred_df$mt_seq <- sapply(1:nrow(varipred_df), function(i) {
    new_seq <- varipred_df$wt_seq[i]
    substr(new_seq, varipred_df$aa_index[i], varipred_df$aa_index[i]) <- varipred_df$mt_aa[i]
    new_seq
  })
  
  write.csv(varipred_df, 'VariPred_files/varipred_input.csv', row.names = F)
}

#' Run VariPred on the variant data or use precalculated output
#'
#' @param variant_data a data frame containing the variant data
#' @param fullLengthSequences a fasta file containing the full length protein sequences related to variants
#' @param use_existing_results a boolean indicating whether VariPred is installed locally
#' @param precalculated_varipred a character string indicating the path to the precalculated VariPred output
#' @export
#' 
varipred_main <- function(use_existing_results = TRUE, variant_data, fullLengthSequences, precalculated_varipred){

  underscore_start("Obtaining VariPred output")
  
  if(isTRUE(use_existing_results)){
    if(!file.exists(precalculated_varipred)){
      stop(paste(' ', error_sign, ' If local varipred installation is not used, please provide the path to the precalculated_varipred file', sep = ''))
    }
    else{
      cat(paste(' ', info_sign, ' Using precalculated VariPred output ', info_sign, '\n', sep = ''))
      precalculated_varipred <- fread(precalculated_varipred, sep = '\t', header = T)
      input_check <- merge(variant_data, precalculated_varipred, by = 'target_id')
      if(nrow(input_check) != nrow(variant_data)){
        stop(paste(error_sign, ' The variant data frame does not match the precalculated varipred output Please check the target_id column', sep = ''))
      }
      else{
        cat(paste(' Annotated variants with VariPred classifications ', completion_sign, '\n', sep = ''))
        variant_data <- input_check
        underscore_finish("VariPred output obtained")
        return(variant_data)
      }
    }
  }
  else{
    cat(paste(' ', calculation_sign, ' Generating VariPred input file ', calculation_sign, '\n', sep = ''))
    generateVariPredInput(variant_data = variant_data, fullLengthSequences = fullLengthSequences)
    if(!file.exists('VariPred_files/varipred_input.csv')){
      stop(paste(' ', error_sign, ' VariPred input file not generated! Please check the input data', sep = ''))
    }
    else{
      cat(paste(' VariPred file generated ', completion_sign, '\n', sep = ''))
    }
    
    if(!dir.exists('VariPred')){
      stop(paste(' ', error_sign, ' The VariPred directory not found! If VariPred installed in your system, please provide the path to the directory', sep = ''))
    }
    setwd('VariPred/VariPred')
    cat(paste(' ', calculation_sign, ' Running VariPred ', calculation_sign, '\n', sep = ''))
    system(paste('python3 main.py -p ../../VariPred_files/ -i varipred_input -o varipred_output', sep = ''))
    setwd('../../')
    
    if(!file.exists('VariPred/example/output_results/varipred_output.txt')){
      stop(paste(' ', error_sign, ' VariPred output not generated! Please check the input data', sep = ''))
    }
    else{
      file.move('VariPred/example/output_results/varipred_output.txt', 'new_run/varipred_output.txt')
      cat(paste(' VariPred output generated and saved in new_run/varipred_output.txt ', completion_sign, '\n', sep = ''))
    }
    varipred_output <- fread('new_run/varipred_output.txt', sep = '\t', header = T)
    colnames(varipred_output)[2] <- 'VariPred_Class'
    
    input_check <- merge(variant_data, varipred_output, by = 'target_id')
    if(nrow(input_check) != nrow(variant_data)){
      cat(paste(' ', warning_sign, ' Dimensions of the variant data does not match the varipred output. ', warning_sign, '\n', sep = ''))
    }
    else{
      cat(paste(' Annotated variants with VariPred classifications ', completion_sign, '\n', sep = ''))
    }
    variant_data <- input_check
    underscore_finish("VariPred output obtained")
    return(variant_data)
  }
}
