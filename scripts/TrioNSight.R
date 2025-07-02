#' @import BiocManager
#' @import data.table
#' @import httr
#' @import jsonlite
#' @import rjson
#' @import Biostrings
#' @import stringr
#' @import dplyr
#' @import stats
#' @import bio3d
#' @import grantham
#' @import adabag
#' @import caret
#' @import flextable
#' @import magrittr
#' 
TrioNSight <- function(variantFile, data_overlap, alignmentFile, seqs2Inherit, InheritanceDetails2File = FALSE,
                       BLOSUM_Matrix = 'BLOSUM80', ScoreConsFile, PrecalculatedScoreCons = TRUE,
                       PrecalculatedVariPred = TRUE, dirVariPred, precalculated_varipred_path, VariPred_output_path, 
                       fullLengthSequences, 
                       PrecalculatedSAAPdap = TRUE, SAAP3_path, 
                       PrecalculatedClusteringResults = FALSE, PrecalculatedClusters = TRUE) {
  
  source('scripts/helper.R')
  
  underscore_start("Starting TrioNsight workflow")
  
  suppressWarnings(variants_inherited <- inherit_variants_using_alignment(alignmentFile = alignmentFile,
                                                                          variantFile = variantFile,
                                                                          inheritingSeqIDs = seqs2Inherit,
                                                                          saveDetailsToFile = InheritanceDetails2File))
  
  suppressWarnings(variant_data_annotated <- get_sequence_based_features(variants_inherited = variants_inherited, 
                                                                         ScoreConsFile = ScoreConsFile,
                                                                         usePrecalculatedScoreCons = PrecalculatedScoreCons, 
                                                                         alignmentFile = alignmentFile,
                                                                         BLOSUM_Matrix = BLOSUM_Matrix))
  
  variant_data_annotated <- protVar_features_wrapper(variant_data = variant_data_annotated, 
                                                     variants = unique(variant_data_annotated$protVar_identifiers), 
                                                     varFormat = 'UniProt', outputFolder = NULL, 
                                                     chunk_size = 100, max_attempts = 10)
  
  
  variant_data_annotated <- varipred_main(use_existing_results = PrecalculatedVariPred, variant_data = variant_data_annotated, 
                                          fullLengthSequences = fullLengthSequences, dirVariPred = dirVariPred, 
                                          precalculated_varipred = precalculated_varipred_path)
  
  variant_data_annotated <- merge(variant_data_annotated, PHACT_TrioN_like[,.(target_id, PHACT_Score, PHACTboost_Score)], by = 'target_id')
  
  suppressWarnings(variant_data_annotated <- saapdap_main(use_existing_results = PrecalculatedSAAPdap, SAAP3_path,
                                                          variant_data = variant_data_annotated, variants_inherited = variants_inherited, 
                                                          Prots2PDBsAnnot = Prots2PDBsAnnot, 
                                                          aa_codes = aa_codes, reverse_aa_codes = reverse_aa_codes, 
                                                          DHDomains2PDBsAnnot = DHDomains2PDBsAnnot, check_models = FALSE))
  
  variant_data_annotated <- variant_clusters_main(use_existing_results = PrecalculatedClusteringResults, use_existing_clusters_only = PrecalculatedClusters, 
                                                  modelsFolder = 'models/', fixed_structure_file = 'O75962_1284_1470_af.pdb', 
                                                  aligned_str_folder = 'models_aligned/',
                                                  reference_clustering = reference_clustering, model_domainboundaries = model_domainboundaries,
                                                  max_aln_attempts = 10, max_aln_cycles = 10, 
                                                  aa_codes = aa_codes, reverse_aa_codes = reverse_aa_codes, 
                                                  strucclus_res_folderpath = 'strucclus_res/', 
                                                  distance_reference_pdb = 'TrioN_DHPH_1290_1594_A_1nty', best_n_reference = 'n5', 
                                                  variant_data_annotated = variant_data_annotated,
                                                  data_overlap = data_overlap)
  
  fwrite(variant_data_annotated, 'new_run/TrioNsight_features_dataset.txt', sep = '\t', quote = FALSE, row.names = FALSE)
  
  underscore_finish("TrioNsight workflow complete. Results saved to 'new_run/TrioNsight_features.txt' file.")
  
  return(variant_data_annotated)
}


TrioNsight_features <- TrioNSight(variantFile = 'resources/TrioNsight_variant_dataset.txt',
                                  alignmentFile = 'resources/Human_DHs_mTMalign.fasta',
                                  seqs2Inherit = c("A1L390_80_293", "O15068_621_815", "O43307_94_309", "O60229_1273_1459", "O60229_1917_2119",
                                                   "O75962_1284_1470", "O75962_1958_2161", "P10911_488_679", "Q58EX7_724_924", "Q86YR7_615_823",
                                                   "Q96PX9_795_1000", "Q9H7P9_101_302", "Q9ULL1_91_312"),
                                  InheritanceDetails2File = FALSE,
                                  BLOSUM_Matrix = 'BLOSUM80', ScoreConsFile = 'resources/Human_DHs_mTMalign.scorecons', PrecalculatedScoreCons = TRUE, 
                                  PrecalculatedVariPred = TRUE, dirVariPred = 'varipred_files/', precalculated_varipred_path = 'stable_run/varipred_output.txt', 
                                  fullLengthSequences = 'resources/Human_DHs_prots_full_length.fasta',
                                  PrecalculatedSAAPdap = TRUE, SAAP3_path = 'structure_analysis/SAAP3/',
                                  data_overlap = trainingSet_Overlaps, PrecalculatedClusteringResults = FALSE, PrecalculatedClusters = FALSE)
