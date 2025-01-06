#' Parse supplied variant data, and multiple sequence alignment to inherit variants to the corresponding sequences
#'
#' @param alignmentFile in fasta format - sequence identifiers must match with the variant data
#' @param variantFile in tsv format - sequence identifiers must match with the alignment file
#' @param ScoreConsFile in tsv format - sequence identifiers must match with the alignment file
#' @return inherited variants data frame
#' @export
inherit_variants_using_alignment <- function(alignmentFile, variantFile, inheritingSeqIDs, saveDetailsToFile){
  
  underscore_start("Beginning Variant Inheritance Process")
  
  cat(paste(' ', warning_sign, ' Be careful not to supply trimmed alignment file ', warning_sign, '\n'), sep = '')
  if(grepl('trim', alignmentFile)){stop(paste(error_sign, ' Potential trimmed alignment file detected. (By the file name...)', sep = ''))}
  
  # ScoreConsFile <- fread(ScoreConsFile)
  alignmentObj <- Biostrings::readAAStringSet(alignmentFile)
  alignment_matrix <- as.matrix(alignmentObj)
  
  variants_df <- read.table(variantFile, sep = '\t', header = TRUE)
  
  cat(paste(' ', calculation_sign, ' Checking if the supplied sequences to inherit to are present in the alignment file ', calculation_sign, '\n', sep = ''))
  if(length(setdiff(inheritingSeqIDs, names(alignmentObj)))){
    cat(paste(' ', warning_sign, ' The following sequences were not present in the alignment file ', warning_sign, '\n', sep = ''))
    print(setdiff(inheritingSeqIDs, names(alignmentObj)))
    cat(paste(' ', warning_sign, 'Number of the remaining sequences to inherit variants to: ', warning_sign, length(intersect(names(alignmentObj), inheritingSeqIDs)), '\n', sep = ''))
  }
  inheritingSeqIDs <- intersect(names(alignmentObj), inheritingSeqIDs)
  if(length(inheritingSeqIDs) == 0){stop('Error: No sequences to inherit variants to. Check the supplied sequence IDs match the alignment file!')}
  
  cat('\n')
  inherited_variants_list <- list()
  for(i in 1:length(inheritingSeqIDs)){
    
    notInAlignment <- t(data.table(c(NA, NA, NA)))
    colnames(notInAlignment) <- c('rowNumber', 'Variant_Identifier', 'Variant_Seq')
    gapAtAlignment <- t(data.table(c(NA, NA, NA, NA, NA)))
    colnames(gapAtAlignment) <- c('rowNumber', 'Variant_Identifier', 'Variant_Seq', 'Location_on_Alignment', 'Inheriting_Seq')
    nativeAAmismatch <- t(data.table(c(NA, NA, NA, NA, NA, NA, NA)))
    colnames(nativeAAmismatch) <- c('rowNumber', 'Variant_Identifier', 'Variant_Sequence', 'Variant_Seq_Native_AA', 'Location_on_Alignment', 'Inheriting_Seq', 'Inherited_Seq_Native_AA')
    
    cat(paste(' ', calculation_sign, ' Inheriting variants to sequence ', i, ' of ', length(inheritingSeqIDs), ' ', calculation_sign, '\n', sep = ''))
    cat(paste(' ', info_sign, ' Sequence ID: ', inheritingSeqIDs[i], ' ', info_sign, '\n', sep = ''))
    
    inheritingSeqID <- inheritingSeqIDs[i]
    inheritingSeq <- alignment_matrix[inheritingSeqID, ]
    inheritingSeqInd <- numeric(length(inheritingSeq))
    inheritIdxModifier <- sapply(strsplit(inheritingSeqID, '_'), function(x) as.integer(x[2]))
    
    count <- 0
    for (j in 1:length(inheritingSeq)) {
      if (inheritingSeq[j] != "-") {
        count <- count + 1
        inheritingSeqInd[j] <- count
      }
      else{
        inheritingSeqInd[j] <- 0
      }
    }
    
    # inherited_variants <- data.frame(Variant_Seq = character(), Variant_Identifier = character(), UniProtNative = character(), UniProtMutant = character(), Variant_Location = integer(), Inherited_Seq = character(), Inherited_Location = integer(), ScoreCons = integer(), nativeSeqCheck = character(), Source = character(), Class = character(), stringsAsFactors = FALSE)
    inherited_variants <- data.frame(Variant_Seq = character(), Variant_Identifier = character(), UniProtNative = character(), UniProtMutant = character(), Variant_Location = integer(), Inherited_Seq = character(), Inherited_Location = integer(), nativeSeqCheck = character(), Source = character(), Class = character(), stringsAsFactors = FALSE)
    
    for (j in 1:nrow(variants_df)) {
      
      variant_source <- variants_df$Source[j]
      if(variants_df$ChrNucChange[j] != '-'){
        CHROM_POS_REF_ALT <- variants_df$ChrNucChange[j]
        CHROM_POS_REF_ALT <- gsub('chr', '', CHROM_POS_REF_ALT, fixed = T)
        CHROM_POS_REF_ALT <- gsub(':', ' ', CHROM_POS_REF_ALT, fixed = T)
        CHROM_POS_REF_ALT <- gsub('_', ' ', CHROM_POS_REF_ALT, fixed = T)
        CHROM_POS_REF_ALT <- gsub('/', ' ', CHROM_POS_REF_ALT, fixed = T)
      }
      else{
        CHROM_POS_REF_ALT <- '-'
      }
      variant_sequence_name <- variants_df$varDomainID[j]
      variant_location <- as.integer(variants_df$protLoc_UniProt[j])
      domain_start <- as.integer(variants_df$begin_dh_uniprot_canonical[j])
      variant_id <- variants_df$Name[j]
      variant_aa_uniprot <- variants_df$AA1_UniProt[j]
      mutant_aa_uniprot <- variants_df$AA2_UniProt[j]
      variant_class <- variants_df$Class[j]
      
      if(!variants_df$varDomainID[j] %in% rownames(alignment_matrix)){
        tmpDT <- data.table(j, variant_id, variant_sequence_name)
        colnames(tmpDT) <- colnames(notInAlignment)
        notInAlignment <- rbind(notInAlignment, tmpDT)
        next
      }
      
      variantSeq <- alignment_matrix[variant_sequence_name,]
      variantSeqVecInd <- numeric(length(variantSeq))
      count <- 0
      for (k in 1:length(variantSeq)) {
        if (variantSeq[k] != "-") {
          count <- count + 1
          variantSeqVecInd[k] <- count
        }
        else{
          variantSeqVecInd[k] <- 0
        }
      }
      
      variant_loc_onMSA <- variant_location - domain_start + 1
      variant_loc_onMSA <- which(variant_loc_onMSA == variantSeqVecInd)
      variant_aa_onMSA <- variantSeq[variant_loc_onMSA]
      inheritingSeq_aa_onMSA <- inheritingSeq[variant_loc_onMSA]
      # ScoreCons_Score <- ScoreConsFile$V1[variant_loc_onMSA]
      
      if(inheritingSeq_aa_onMSA != '-'){
        if(inheritingSeq_aa_onMSA != variant_aa_uniprot){
          tmpDT <- data.table(j, variant_id, variant_sequence_name, variant_aa_uniprot, variant_loc_onMSA, inheritingSeqID, inheritingSeq_aa_onMSA)
          colnames(tmpDT) <- colnames(nativeAAmismatch)
          nativeAAmismatch <- rbind(nativeAAmismatch, tmpDT)
          next
        }
      }
      if(variant_aa_onMSA == '-'){
        stop(paste(' ', error_sign, ' Query sequence has a gap at this location on MSA. This is an unrecoverable error!', sep = ''))
      }
      if(inheritingSeq_aa_onMSA == '-'){
        tmpDT <- data.table(j, variant_id, variant_sequence_name, variant_loc_onMSA, inheritingSeqID)
        colnames(tmpDT) <- colnames(gapAtAlignment)
        gapAtAlignment <- rbind(gapAtAlignment, tmpDT)
        next
      }
      else{
        inheritInd <- inheritingSeqInd[variant_loc_onMSA]
        variantSeq_native_inheritingSeq_native <- paste(variant_aa_onMSA, inheritingSeq_aa_onMSA, sep = '->')
        # inherited_variants <- rbind(inherited_variants, data.frame(Variant_Seq = variant_sequence_name, CHROM_POS_REF_ALT = CHROM_POS_REF_ALT, Variant_Identifier = variant_id, UniProtNative = variant_aa_uniprot, UniProtMutant = mutant_aa_uniprot, Variant_Location = variant_location, Inherited_Seq = inheritingSeqID, Inherited_Location = inheritInd+inheritIdxModifier-1, ScoreCons = ScoreCons_Score, MSALocation = variant_loc_onMSA, nativeSeqCheck = variantSeq_native_inheritingSeq_native, Source = variant_source, Class = variant_class))
        inherited_variants <- rbind(inherited_variants, data.frame(Variant_Seq = variant_sequence_name, CHROM_POS_REF_ALT = CHROM_POS_REF_ALT, Variant_Identifier = variant_id, UniProtNative = variant_aa_uniprot, UniProtMutant = mutant_aa_uniprot, Variant_Location = variant_location, Inherited_Seq = inheritingSeqID, Inherited_Location = inheritInd+inheritIdxModifier-1, MSALocation = variant_loc_onMSA, nativeSeqCheck = variantSeq_native_inheritingSeq_native, Source = variant_source, Class = variant_class))
        inherited_variants <- as.data.table(inherited_variants)
      }
      
    }
    notInAlignment <- notInAlignment[-1,]
    gapAtAlignment <- gapAtAlignment[-1,]
    nativeAAmismatch <- nativeAAmismatch[-1,]
    inherited_variants_list[[i]] <- list(inherited_variants, notInAlignment, gapAtAlignment, nativeAAmismatch)
    names(inherited_variants_list[[i]]) <- c('Inherited_Variants', 'Variant_sequence_not_in_alignment', 'Variant_coincides_with_gap_for_inheriting_sequence', 'Native_amino_acid_mismatch_at_the_given_position')
    
    cat(paste(' ', completion_sign, ' Variant inheritance to sequence: ', inheritingSeqID, ' ', completion_sign, '\n', sep = ''))
    cat(paste(' ', info_sign, ' Total variants inherited to this sequence: ', nrow(inherited_variants), info_sign, '\n', sep = ''))
    cat(paste(' ', info_sign, ' Variants removed because sequence associated with variant not found in the alignment: ', nrow(notInAlignment), info_sign, '\n', sep = ''))
    cat(paste(' ', info_sign, ' Variants coinciding with a gap in the inheriting sequence: ', nrow(gapAtAlignment), info_sign, '\n', sep = ''))
    cat(paste(' ', info_sign, ' Variants with different native amino acids in variant sequence and inheriting sequence: ', nrow(nativeAAmismatch), info_sign, '\n', sep = ''))

  }
  names(inherited_variants_list) <- inheritingSeqIDs
  variants_inherited_dt <- lapply(inherited_variants_list, function(x){x[[1]]})
  variants_inherited_dt <- do.call(rbind, variants_inherited_dt)
  variants_inherited_dt <- as.data.frame(variants_inherited_dt)
  colnames(variants_inherited_dt)[4:5] <- c('REF', 'ALT')
  
  variants_inherited_dt$duplicateChecker <- paste(variants_inherited_dt$Inherited_Seq, '_', variants_inherited_dt$REF, variants_inherited_dt$MSALocation, variants_inherited_dt$ALT, sep = '')
  variants_inherited_duplicated <- variants_inherited_dt[duplicated(variants_inherited_dt$duplicateChecker),]
  variants_inherited_dt <- variants_inherited_dt[!duplicated(variants_inherited_dt$duplicateChecker),]
  variants_inherited_dt$duplicateChecker <- NULL
  cat(paste(' ', info_sign, ' Variants that were not inherited because the same variant exists in the inheriting sequence: ', nrow(variants_inherited_duplicated), info_sign, '\n', sep = ''))
  
  if(isTRUE(saveDetailsToFile)){
    
    notInAlignment <- lapply(inherited_variants_list, function(x){x[[2]]})
    notInAlignment <- as.data.table(unique(do.call(rbind, notInAlignment)))
    foundGapAtAlignment <- lapply(inherited_variants_list, function(x){x[[3]]})
    foundGapAtAlignment <- as.data.table(unique(do.call(rbind, foundGapAtAlignment)))
    nativeMismatch <- lapply(inherited_variants_list, function(x){x[[4]]})
    nativeMismatch <- as.data.table(unique(do.call(rbind, nativeMismatch)))
    
    fwrite(variants_inherited_dt, file = 'new_run/Inherited_Variants.txt', sep = '\t')
    fwrite(notInAlignment, file = 'new_run/Variant_sequence_not_in_alignment.txt', sep = '\t')
    fwrite(foundGapAtAlignment, file = 'new_run/Variant_coincides_with_gap_for_inheriting_sequence.txt', sep = '\t')
    fwrite(nativeMismatch, file = 'new_run/Native_amino_acid_mismatch_at_the_given_position.txt', sep = '\t')
    
    cat('\n')
    cat(paste(' ', completion_sign, ' Detailed variant inheritance data is saved to files ', completion_sign, '\n', sep = ''))
  }
  
  underscore_finish("Variant Inheritance Process Finished")
  
  return(variants_inherited_dt)
}
