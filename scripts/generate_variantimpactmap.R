if (!requireNamespace("Biostrings", quietly = TRUE)) {
  install.packages('Biostrings')   
}
library(BiocManager)
if (!requireNamespace("dplyr", quietly = TRUE)) {
  BiocManager::install('dplyr')
}
library(data.table)
if (!requireNamespace("tidyr", quietly = TRUE)) {
  BiocManager::install('tidyr')
}

#' Function to generate all possible single amino acid substitutions from input sequence
#'
#' @param sequence Amino acid sequence of the target protein
#' @param start Start position of the target protein sequence to start generating substitutions
#' @param end End position of the target protein sequence to stop generating substitutions
#' @return A data frame containing all possible single amino acid substitutions
#' @export
#' 
generate_substitutions_from_sequence <- function(sequence, start, end) {
  aa_list <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  
  positions <- start:end
  original_aa <- strsplit(substr(sequence, start, end), "")[[1]]
  
  substitutions <- data.frame(
    position = rep(positions, each = length(aa_list)),
    original = rep(original_aa, each = length(aa_list)),
    substitution = rep(aa_list, times = length(positions))
  )
  
  substitutions <- substitutions[substitutions$original != substitutions$substitution, ]
  
  substitutions$wt_seq <- sequence
  substitutions$mt_seq <- sapply(1:nrow(substitutions), function(i) {
    new_seq <- sequence
    substr(new_seq, substitutions$position[i], substitutions$position[i]) <- substitutions$substitution[i]
    new_seq
  })
  
  return(substitutions)
}

#' Function to convert single letter amino acid code to three letter code
#' @param aa Amino acid in single letter code
#' @return Amino acid in three letter code
#' @export
#' 
aa_to_three_letter <- function(aa) {
  aa_dict <- c(
    "A" = "Ala", "R" = "Arg", "N" = "Asn", "D" = "Asp", "C" = "Cys",
    "Q" = "Gln", "E" = "Glu", "G" = "Gly", "H" = "His", "I" = "Ile",
    "L" = "Leu", "K" = "Lys", "M" = "Met", "F" = "Phe", "P" = "Pro",
    "S" = "Ser", "T" = "Thr", "W" = "Trp", "Y" = "Tyr", "V" = "Val"
  )
  return(aa_dict[aa])
}

#' Main script to generate variant impact map - writes the output to a txt file
#' @export
#' 
generate_variant_impact_map <- function(domain_boundaries_file = 'resources/domain_boundaries.csv', 
                                        domain_sequences_fasta = 'resources/Human_DHs_prots_full_length.fasta', 
                                        sequence_id_annotation_file = 'resources/Human_DH_Annotations.txt', 
                                        output_file = 'resources/variant_impact_map.txt') {
  
  if(!file.exists(domain_boundaries_file)){
    stop(paste(' ', error_sign, ' path must contain domain boundaries file: ', domain_boundaries_file, ' !', sep = ''))
  }
  # Read input file containing UniProt IDs and ranges
  domain_boundaries_data <- read.csv(domain_boundaries_file, header = TRUE, stringsAsFactors = FALSE)
  
  if(!file.exists(domain_sequences_fasta)){
    stop(paste(' ', error_sign, ' path must contain domain sequences file: ', domain_sequences_fasta, ' !', sep = ''))
  }
  # Read FASTA file
  domain_sequences_fasta <- readAAStringSet("Human_DHs_prots_full_length.fasta")
  
  if(!file.exists(sequence_id_annotation_file)){
    stop(paste(' ', error_sign, ' path must contain sequences id annotations file: ', sequence_id_annotation_file, ' !', sep = ''))
  }
  # Read sequence ID annotation file
  sequence_id_annotation_file <- fread('Human_DH_Annotations.txt', header = TRUE, sep = '\t')
  
  cat(paste(' ', calculation_sign, ' Generating all possible substitutions for each given sequence ', calculation_sign, '\n', sep = ''))
  # Process each UniProt ID and range
  output <- lapply(1:nrow(domain_boundaries_data), function(i) {
    
    track_progress(i, nrow(domain_boundaries_data))
    
    gene_symbol <- domain_boundaries_data$gene_symbol[i]
    uniprot_id <- domain_boundaries_data$uniprot_id[i]
    start <- domain_boundaries_data$start[i]
    end <- domain_boundaries_data$end[i]
    
    # Find the sequence for the current UniProt ID
    seq_index <- which(names(domain_sequences_fasta) == uniprot_id)
    
    if (length(seq_index) == 0) {
      cat("Sequence not found for UniProt ID:", uniprot_id, "\n")
      return(NULL)
    }
    
    sequence <- as.character(domain_sequences_fasta[[seq_index]])
    
    # Generate substitutions
    substitutions <- generate_substitutions_from_sequence(sequence, start, end)
    
    # Add additional columns
    substitutions %>%
      mutate(
        Accession = uniprot_id,
        GeneSymbol = gene_symbol,
        target_id = paste0(Accession, "_", original, position, substitution),
        AA1_VarDB = original,
        AA_VarDB_v1 = paste0(original, position, substitution),
        AA2_VarDB = substitution,
        protLoc_VarDB = position,
        protLoc_UniProt = position,
        AA_VarDB_v2 = paste(position, original, substitution),
        mapToDomain = "DH",
        AA_VarDB_v3 = paste0(aa_to_three_letter(original), position, aa_to_three_letter(substitution)),
        begin_dh_uniprot_canonical = start,
        end_dh_uniprot_canonical = end,
        begin_dh_VarDB = start,
        end_dh_VarDB = end,
        AA_UniProt_v1 = paste0(original, position, substitution),
        AA_UniProt_v2 = paste(position, original, substitution),
        AA_UniProt_v3 = paste0(aa_to_three_letter(original), position, aa_to_three_letter(substitution)),
        AA1_UniProt = original,
        AA2_UniProt = substitution,
        varDomainID = paste0(Accession, "_", start, "_", end),
        Source = 'Synthetic',
        Class = '-',
        rsID = '-',
        Name = '-',
        Phenotype = '-',
        Reference = '-',
        Allele_Count = '-',
        Allele_Frequency = '-',
        Nuc_Change = '-',
        LastEvaluated = '-',
        ChrNucChange = '-',
        Chr_v1_GRCh38 = '-',
        Chr_v2_GRCh38 = '-'
      )
  }) %>%
    bind_rows()
  
  output <- merge(output, sequence_id_annotation_file, by = "GeneSymbol")
  output <- as.data.table(output)
    
  cat(paste(' ', calculation_sign, ' Generating varipred input file for all possible substitutions ', calculation_sign, '\n', sep = ''))
  
  varipred_file <- output[,.(target_id, position, original, substitution, wt_seq, mt_seq)]
  colnames(varipred_file) <- c('target_id', 'aa_index', 'wt_aa', 'mt_aa', 'wt_seq', 'mt_seq')
  write.csv(varipred_file, file = 'varipred_file.csv', row.names = F)
  cat(paste(' Varipred file is written to "varipred_file.csv" ', completion_sign, '\n', sep = ''))
  
  output[, wt_seq := NULL]
  output[, mt_seq := NULL]
  output[, position := NULL]
  output[, original := NULL]
  output[, substitution := NULL]
  
  # Write output to file
  fwrite(output, file = paste("resources/", output_file, sep = ''), sep = '\t')
  cat(paste(' Impact map has been written to ', paste("resources/", output_file, sep = ''), completion_sign, '\n', sep = ''))
}

#' Run the main function
#' @export
generate_variant_impact_map()