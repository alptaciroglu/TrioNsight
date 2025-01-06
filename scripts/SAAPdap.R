#' Function to check the structural models for SAAPdap run.
#' 
#' The function checks the residue numbers and chains in the PDB files and corrects them if necessary.
#' @export
#' 
check_models_for_SAAPdap <- function(){
  
  if(!dir.exists('structure_analysis/')){
    stop(paste(' ', error_sign, ' Working directory must contain "structure_analysis" directory!', sep = ''))
  }
  setwd('structure_analysis/')
  if(!dir.exists('models/')){
    stop(paste(' ', error_sign, ' structure_analysis folder must contain "models/" directory!', sep = ''))
  }
  setwd('models/')
  if(!dir.exists('DH/')){
    stop(paste(' ', error_sign, ' models directory must contain "DH/", "DHPH/" and "Experimental/" folders separately. Folder "DH/" not found!', sep = ''))
  }
  if(!dir.exists('DHPH/')){
    stop(paste(' ', error_sign, ' models directory must contain "DH/", "DHPH/" and "Experimental/" folders separately. Folder "DHPH/" not found!', sep = ''))
  }
  if(!dir.exists('Experimental/')){
    stop(paste(' ', error_sign, ' Directory to models must contain "DH/", "DHPH/" and "Experimental/" folders separately. Folder "Experimental/" not found!', sep = ''))
  }
  
  cat(paste(' ', info_sign, ' Residue numbers in the PDB files must start with 1 ', info_sign, '\n', sep = ''))
  cat(paste(' ', info_sign, ' PDB files must contain only chain "A" ', info_sign, '\n', sep = ''))
  
  setwd('DH/')
  DH_Models <- list.files()
  cat(paste(' ', calculation_sign, ' Checking the DH models for residue numbers and chains ', calculation_sign, '\n', sep = ''))
  
  for(i in 1:length(DH_Models)){
    
    pdb_start <- as.numeric(as.character(sapply(strsplit(DH_Models[i], '_', fixed = T), '[', 2)))
    tmpModel <- bio3d::read.pdb(DH_Models[i])
    
    if(tmpModel$atom$resno[1] != 1){
      tmpModel$atom$resno <- tmpModel$atom$resno + (1 - tmpModel$atom$resno[1])
    }
    tmpModel$atom$eleno <- 1:nrow(tmpModel$atom)
    
    bio3d::write.pdb(tmpModel, file = DH_Models[i])
  }
  
  cat(paste(' PDB files checked for DH models ', attempts, ' ', completion_sign, '\n', sep = ''))
  
  setwd('../DHPH/')
  DHPH_Models <- list.files()
  cat(paste(' ', calculation_sign, ' Checking the DHPH models for residue numbers and chains ', calculation_sign, '\n', sep = ''))
  
  for(i in 1:length(DHPH_Models)){
    
    pdb_start <- as.numeric(as.character(sapply(strsplit(DHPH_Models[i], '_', fixed = T), '[', 2)))
    tmpModel <- bio3d::read.pdb(DHPH_Models[i])
    
    if(tmpModel$atom$resno[1] != 1){
      tmpModel$atom$resno <- tmpModel$atom$resno + (1 - tmpModel$atom$resno[1])
    }
    tmpModel$atom$eleno <- 1:nrow(tmpModel$atom)
    
    bio3d::write.pdb(tmpModel, file = DHPH_Models[i])
  }
  
  cat(paste(' PDB files checked for DHPH models ', attempts, ' ', completion_sign, '\n', sep = ''))
  
  setwd('../Experimental/')
  Experimental_Models <- list.files()
  cat(paste(' ', calculation_sign, ' Checking the Experimental models for residue numbers and chains ', calculation_sign, '\n', sep = ''))
  
  for(i in 1:length(Experimental_Models)){
    
    pdb_start <- as.numeric(unlist(lapply(strsplit(Experimental_Models[i], '_', fixed = T), function(x){x[length(x)-3]})))
    tmpModel <- bio3d::read.pdb(Experimental_Models[i])
    
    protChain <- unlist(lapply(strsplit(Experimental_Models[i], '_', fixed = T), function(x){x[length(x)-1]}))
    chainsTotal <- unique(tmpModel$atom$chain)
    
    if(length(chainsTotal) > 1){
      tmpModel$atom <- tmpModel$atom[tmpModel$atom$chain == protChain,]
    }
    
    if(tmpModel$atom$resno[1] != 1){
      tmpModel$atom$resno <- tmpModel$atom$resno + (1 - tmpModel$atom$resno[1])
    }
    tmpModel$atom$eleno <- 1:nrow(tmpModel$atom)
    
    bio3d::write.pdb(tmpModel, file = Experimental_Models[i])
  }
  
  cat(paste(' PDB files checked for Experimental models ', attempts, ' ', completion_sign, '\n', sep = ''))
  setwd('../../../')
  
}

#' Function to convert inherited variant data to SAAPdap format
#' 
#' @param variants_inherited_list Variant data inherited to other proteins in the MSA in data.frame format
#' @param Prots2PDBsAnnot data.frame that holds information about the PDB files for each protein in the variants_inherited_list
#' @return variant_data_in_SAAPdap_format: list of data.tables that hold the variant data in SAAPdap format
#' @export
#' 
## Map variants to different models based on their location in the PDB file
## There are some models that have different domain start locations than what should be the default
## Those models are listed below in the script but I don't know why they are behaving this way.
inherited_variants_2_SAAPdap_format <- function(variants_inherited_list, Prots2PDBsAnnot){
  
  cat(paste(' ', calculation_sign, ' Converting variant data to SAAPdap data format ', calculation_sign, '\n', sep = ''))
  cat(paste(' ', info_sign, ' The resulting files will have the following columns in order ', info_sign, '\n', sep = ''))
  cat(paste(' ', info_sign, ' "UniProt Accession" "Native Residue" "Residue Number (overall protein sequence)" "Mutant Residue" "PDB file"  "Residue Number (PDB file)" ', info_sign, '\n', sep = ''))
  
  variant_data_in_SAAPdap_format <- list()
  
  for(i in 1:length(variants_inherited_list)){
    
    variants_inherited_nextProt <- variants_inherited_list[[i]]
    
    SAAPdap_format_df <- as.data.table(matrix(ncol = 6, nrow = nrow(variants_inherited_nextProt), data = ''))
    colnames(SAAPdap_format_df) <- c('UniProtAC', 'native', 'resnum', 'mutant', 'pdbfile', 'resid')
    
    UniAC <- sapply(strsplit(variants_inherited_nextProt$Variant_Seq, '_', fixed = T), '[', 1)
    SAAPdap_format_df[, UniProtAC := UniAC]
    
    resNum <- variants_inherited_nextProt$Variant_Location
    firstAA <- variants_inherited_nextProt$wt_3letter
    secondAA <- variants_inherited_nextProt$mt_3letter
    variant_by_inheritance <- variants_inherited_nextProt$variant_by_inheritance
    
    SAAPdap_format_df[, native := firstAA]
    SAAPdap_format_df[, resnum := resNum]
    SAAPdap_format_df[, mutant := secondAA]
    SAAPdap_format_df[, variant_by_inheritance := variant_by_inheritance]
    
    nextPDB <- paste(names(variants_inherited_list)[i], '_af.pdb', sep = '')
    allPDB_fromProt <- as.character(Prots2PDBsAnnot[AF_DH == nextPDB][,.(AF_DH, AF_DHPH, AF_Full, Experimental)])
    allPDB_fromProt <- allPDB_fromProt[allPDB_fromProt != '']
    SAAPdap_format_template <- SAAPdap_format_df
    
    for(j in 1:length(allPDB_fromProt)){
      
      if(j == 4){ # Goes into experimental pdb files
        pdb_exp_iter <- strsplit(allPDB_fromProt[j], ', ', fixed = T)[[1]]
        for(k in 1:length(pdb_exp_iter)){
          pdb_exp_next <- pdb_exp_iter[k]
          SAAPdap_format_df[, pdbfile := pdb_exp_next]
          pdb_start <- lapply(strsplit(pdb_exp_next, '_', fixed = T), function(x){x[length(x)-3]})[[1]]
          pdb_start <- as.numeric(pdb_start)
          
          SAAPdap_format_df[, resid := variants_inherited_nextProt$Inherited_Location - pdb_start + 1]
          SAAPdap_format_df[, resid := paste('A', resid, sep = '')]
          SAAPdap_format_df <- merge(SAAPdap_format_df, unique(Prots2PDBsAnnot[,.(UniProtAC, UniProtID)]), by = 'UniProtAC')
          
          pdbName <- gsub('.pdb', '', pdb_exp_next, fixed = T)
          variant_data_in_SAAPdap_format[[pdbName]] <- SAAPdap_format_df
          SAAPdap_format_df <- SAAPdap_format_template
        }
      }
      else{
        SAAPdap_format_df[, pdbfile := allPDB_fromProt[j]]
        pdb_start <- as.character(sapply(strsplit(allPDB_fromProt[j], '_', fixed = T), '[', 2))
        if(pdb_start == 'full'){
          pdb_start <- 1
        }
        else{
          pdb_start <- as.numeric(pdb_start)
        }
        if(allPDB_fromProt[j] %in% c("Q58EX7_724_924_af.pdb", "Q86YR7_615_823_af.pdb", "Q86YR7_615_974_af.pdb", "Q9H7P9_101_302_af.pdb")){
          SAAPdap_format_df[, resid := variants_inherited_nextProt$Inherited_Location - pdb_start]
        }
        else{
          SAAPdap_format_df[, resid := variants_inherited_nextProt$Inherited_Location - pdb_start + 1]
        }
        SAAPdap_format_df[, resid := paste('A', resid, sep = '')]
        SAAPdap_format_df <- merge(SAAPdap_format_df, unique(Prots2PDBsAnnot[,.(UniProtAC, UniProtID)]), by = 'UniProtAC')
        
        pdbName <- gsub('_af.pdb', '', allPDB_fromProt[j], fixed = T)
        pdbName <- gsub('.pdb', '', pdbName, fixed = T)
        pdbName <- paste(pdbName, '_af', sep = '')
        variant_data_in_SAAPdap_format[[pdbName]] <- SAAPdap_format_df
        SAAPdap_format_df <- SAAPdap_format_template
      }
    }
  }
  cat(paste(' Variants converted to SAAPdap data format ', completion_sign, '\n', sep = ''))
  
  return(variant_data_in_SAAPdap_format)
}

#' Wrapper function to prepare the data for SAAPdap pipeline
#' 
#' The files created this way can be supplied to SAAPdap pipeline and used to predict the effect of variants on protein stability
#' @param variants_inherited Variant data inherited to other proteins in the MSA in data.frame format
#' @param Prots2PDBsAnnot data.frame that holds information about the PDB files for each protein in the variants_inherited_list
#' @param aa_codes Vector to convert 1-letter amino acid codes to 3-letter amino acid codes
#' @param reverse_aa_codes Vector to convert 3-letter amino acid codes to 1-letter amino acid codes
#' @export
#' 
prepare_data_for_SAAPdap <- function(variants_inherited, Prots2PDBsAnnot, aa_codes, reverse_aa_codes){
  
  if(!dir.exists('structure_analysis/')){
    stop(paste(' ', error_sign, ' Working directory must contain "structure_analysis" directory!', sep = ''))
  }
  setwd('structure_analysis/')
  if(!dir.exists('models/')){
    stop(paste(' ', error_sign, ' structure_analysis folder must contain "Models" directory!', sep = ''))
  }
  setwd('models/')
  if(!dir.exists('DH/')){
    stop(paste(' ', error_sign, ' Models directory must contain "DH", "DHPH" and "Experimental" folders separately. Folder "DH" not found!', sep = ''))
  }
  if(!dir.exists('DHPH/')){
    stop(paste(' ', error_sign, ' Models directory must contain "DH", "DHPH" and "Experimental" folders separately. Folder "DHPH" not found!', sep = ''))
  }
  if(!dir.exists('Full/')){
    stop(paste(' ', error_sign, ' Models directory must contain "DH", "DHPH" and "Experimental" folders separately. Folder "Full" not found!', sep = ''))
  }
  if(!dir.exists('Experimental/')){
    stop(paste(' ', error_sign, ' Directory to models must contain "DH", "DHPH" and "Experimental" folders separately. Folder "Experimental" not found!', sep = ''))
  }
  if(!dir.exists('../saapdap_files/')){
    dir.create('../saapdap_files/')
    cat(paste(' ', calculation_sign, ' Creating "saapdap_files" directory to store all structural models ', calculation_sign, '\n', sep = ''))  
  }
  
  ## Read PDB files into a list
  pdbList <- list()
  
  modelTypes <- c('DH', 'DHPH', 'Full', 'Experimental')
  
  cat(paste(' ', calculation_sign, ' Reading structural models ', calculation_sign, '\n', sep = ''))
  cat(paste(' ', info_sign, ' This might take a few minutes ', info_sign, '\n', sep = ''))
  
  warningBool <- TRUE
  for(i in 1:length(modelTypes)){
    
    setwd(modelTypes[i])
    tmpModels <- list.files(pattern = '.pdb')
    for(j in 1:length(tmpModels)){
      
      pdbList[[tmpModels[j]]] <- bio3d::read.pdb(tmpModels[j])
      
    }
    current_models <- list.files(pattern = '.pdb')
    all_models <- list.files(pattern = '.pdb', path = '../../saapdap_files/')
    if(length(setdiff(current_models, all_models)) > 0){
      file.copy(from = current_models, to = '../../saapdap_files/', 
                overwrite = TRUE, recursive = FALSE, 
                copy.mode = TRUE)
      if(isTRUE(warningBool)){
        cat(paste(' ', info_sign, ' Copying models to "saapdap_files" directory ', info_sign, '\n', sep = ''))
        cat(paste(' ', info_sign, ' This will be only done for the missing models in the "saapdap_files" directory ', info_sign, '\n', sep = ''))
        warningBool <- FALSE
      }
    }
    setwd('../')
  }
  
  setwd('../')
  cat(paste(' PDB files read ', completion_sign, '\n', sep = ''))
  
  ## Control if the residues are the same in the PDB file and the SAAPdap file
  ## Then write the SAAPdap files with correctly projected to native residues to a directory
  
  cat(paste(' ', calculation_sign, ' Mapping inherited & reported variants to the structural models ', calculation_sign, '\n', sep = ''))
  cat(paste(' ', info_sign, ' Residue in model must match the residue in the protein sequence ', info_sign, '\n', sep = ''))
  
  ## variants_inherited <- do.call(rbind, variants_inherited)
  colnames(variants_inherited)[4:5] <- c('wt', 'mt')
  variants_inherited <- convert_aa_codes(variants_inherited, aa_codes, reverse_aa_codes, column_name = "wt")
  variants_inherited <- convert_aa_codes(variants_inherited, aa_codes, reverse_aa_codes, column_name = "mt")
  
  variants_inherited <- variants_inherited %>% group_by(Inherited_Seq) %>% group_split()
  variants_inherited <- lapply(variants_inherited, function(x){as.data.frame(x)})
  names(variants_inherited) <- unlist(lapply(variants_inherited, function(x){unique(x$Inherited_Seq)}))
  for(i in 1:length(variants_inherited)){
    variants_inherited[[i]]$variant_by_inheritance <- ''
    for(j in 1:nrow(variants_inherited[[i]])){
      wt <- variants_inherited[[i]]$wt[j]
      inherited_loc <- variants_inherited[[i]]$Inherited_Location[j]
      mt <- variants_inherited[[i]]$mt[j]
      inheriting_seq <- sapply(strsplit(names(variants_inherited)[i], '_', fixed = T), '[', 1)
      variant_by_inheritance <- paste(inheriting_seq, '_', wt, inherited_loc, mt, sep = '')
      variants_inherited[[i]]$variant_by_inheritance[j] <- variant_by_inheritance
    }
  }
  
  cat(paste(' Variants mapped to structural models ', completion_sign, '\n', sep = ''))
  
  variant_data_in_SAAPdap_format <- inherited_variants_2_SAAPdap_format(variants_inherited_list = variants_inherited, Prots2PDBsAnnot = Prots2PDBsAnnot)
  
  cat(paste(' ', calculation_sign, ' Writing input files necessary for the SAAPdap run ', calculation_sign, '\n', sep = ''))
  cat(paste(' ', info_sign, ' One file per model will be written ', info_sign, '\n', sep = ''))
  
  variant_data_as_SAAPdap_files <- list()
  for(i in 1:length(variant_data_in_SAAPdap_format)){
    
    nextPDB <- paste(names(variant_data_in_SAAPdap_format)[i], '.pdb', sep = '')
    tmpAccession <- sapply(strsplit(names(variant_data_in_SAAPdap_format)[i], '_', fixed = T), '[', 1)
    tmpVarObjFit <- variant_data_in_SAAPdap_format[[i]]
    tmpVarObjFit[, resid := as.numeric(gsub('A', '', resid, fixed = T))]
    tmpVarObjFit[, native := toupper(native)]
    
    pdbFile <- pdbList[[nextPDB]]$atom
    pdbFile <- as.data.table(pdbFile)
    pdbFile <- pdbFile[,.(resno, resid)]
    pdbFile <- unique(pdbFile)
    colnames(pdbFile) <- c('resid', 'PDB_AA')
    
    tmpVarObjFit <- merge(tmpVarObjFit, pdbFile, by = 'resid', all.x = T)
    if(nrow(tmpVarObjFit) > 0){
      tmpVarObjFit <- tmpVarObjFit[native == PDB_AA,]
      if(nrow(tmpVarObjFit) > 0){
        variant_data_as_SAAPdap_files[[nextPDB]] <- tmpVarObjFit
      }
    }
  }
  
  setwd('saapdap_files/')
  for(i in 1:length(variant_data_as_SAAPdap_files)){
    
    modelName <- gsub('.pdb', '', names(variant_data_as_SAAPdap_files)[i], fixed = T)
    fileName <- paste(modelName, '.txt', sep = '')
    variant_data_as_SAAPdap_files[[i]] <- variant_data_as_SAAPdap_files[[i]][, .(UniProtAC, native, resnum, mutant, pdbfile, resid, UniProtID)]
    if(!grepl('A', variant_data_as_SAAPdap_files[[i]]$resid[1], fixed = T)){
      variant_data_as_SAAPdap_files[[i]][, resid := paste('A', resid, sep = '')]
    }
    fwrite(variant_data_as_SAAPdap_files[[i]], fileName, sep = '\t', quote = F)
    
  }
  cat(paste(' Data is written to files! You can check these input files in "saapdap_files/" folder', completion_sign, '\n', sep = ''))
  
  setwd('../')
  if(dir.exists('saapdap_res')){
    cat(paste(' ', warning_sign, ' Output "saapdap_res" folder from a previous SAAPdap run found! Removing this folder to create a new run. ', warning_sign, '\n', sep = ''))
    unlink("saapdap_res", recursive = TRUE)
  }
  dir.create('saapdap_res')
  setwd('../')
}

#' Function to parse SAAPdap results
#' 
#' The files created this way can be supplied to SAAPdap pipeline and used to predict the effect of variants on protein stability
#' @param aa_codes Vector to convert 1-letter amino acid codes to 3-letter amino acid codes
#' @return SAAPdapResults.dt: a data.table that contains the results of the SAAPdap pipeline
#' @export
#' 
parse_SAAPdap_results <- function(use_existing_results, aa_codes, DHDomains2PDBsAnnot){
  
  cat(paste(' ', calculation_sign, ' Parsing SAAPdap results ', calculation_sign, '\n', sep = ''))
  
  if(isTRUE(use_existing_results)){
    cat(paste(' ', calculation_sign, ' Using results from the stable run ', calculation_sign, '\n', sep = ''))
    setwd("stable_run/saapdap_res/")
  } else {
    cat(paste(' ', calculation_sign, ' Using results from a new run ', calculation_sign, '\n', sep = ''))
    setwd("new_run/saapdap_res/")
  }
  
  modelsFolders <- list.files()
  for(i in 1:length(modelsFolders)){
    setwd(modelsFolders[i])
    file.remove('UniProtAC_native_resnum_mutant.json')
    setwd('../')
  }
  
  aa_codes <- toupper(aa_codes)
  
  ExternalProgramFailed <- ''
  GlycinesNotOK <- ''
  ProlinesNotOK <- ''
  
  for(i in 36:length(modelsFolders)){
    
    setwd(modelsFolders[i])
    
    variantsFiles <- list.files(pattern = '.json$')
    for(j in 1:length(variantsFiles)){
      
      tmpJSON <-  rjson::fromJSON(file = variantsFiles[j], simplify = T)
      tmpRes <- tmpJSON$SAAPS$pdbs[[1]]$SAAP$results
      
      if(is.null(tmpRes)){
        tmpModel <- modelsFolders[i]
        tmpVariant <- gsub('.json', '', variantsFiles[j], fixed = T)
        tmpRow <- c(tmpModel, tmpVariant, rep(NA, 48))
        tmpRow <- t(as.data.frame(tmpRow))
        SAAPdapResults.dt <- rbind(SAAPdapResults.dt, tmpRow)
        next
      }
      Model <- tmpJSON$SAAPS$pdbs[[1]]$SAAP$file
      Model <- gsub('.pdb', '', Model, fixed = T)
      
      Accession <- tmpJSON$SAAPS$uniprotac
      Native <- tmpJSON$SAAPS$native
      Native <- names(aa_codes[aa_codes == Native])
      ResNum <- tmpJSON$SAAPS$resnum
      Mutant <- tmpJSON$SAAPS$mutant
      Mutant <- names(aa_codes[aa_codes == Mutant])
      target_id <- paste(Accession, '_', Native, ResNum, Mutant, sep = '')
      
      Binding <- as.character(tmpRes$Binding)
      
      SProtFTs <- as.character(c(tmpRes$SProtFT, tmpRes$SprotFTdist))
      if(sum(grepl('External program failed', SProtFTs)) > 0){
        ExternalProgramFailed <- c(ExternalProgramFailed, paste(modelsFolders[i], ' -> ', variantsFiles[j], sep = ''))
        SProtFTs <- as.character(t(matrix(ncol = 1, nrow = 14, NA)))
      }
      else{
        SProtFT <- tmpRes$SProtFT
        SprotFTdist <- tmpRes$SprotFTdist
        if(!is.null(SProtFT$`SProtFT-ERROR`)){
          SProtFTs <- c(NA, NA, NA, SprotFTdist$`SprotFTdist-BOOL`, 
                        SprotFTdist$`SprotFTdist-ACT_SITE`, SprotFTdist$`SprotFTdist-BINDING`, SprotFTdist$`SprotFTdist-CA_BIND`, 
                        SprotFTdist$`SprotFTdist-DNA_BIND`, SprotFTdist$`SprotFTdist-NP_BIND`, SprotFTdist$`SprotFTdist-METAL`, 
                        SprotFTdist$`SprotFTdist-MOD_RES`, SprotFTdist$`SprotFTdist-CARBOHYD`, SprotFTdist$`SprotFTdist-MOTIF`, 
                        SprotFTdist$`SprotFTdist-LIPID`)
        }
        else{
          SProtFTs <- c(SProtFT$`SProtFT-FEATURES`, SProtFT$`SProtFT-BOOL`, SProtFT$`SProtFT-NAMES`, SprotFTdist$`SprotFTdist-BOOL`, 
                        SprotFTdist$`SprotFTdist-ACT_SITE`, SprotFTdist$`SprotFTdist-BINDING`, SprotFTdist$`SprotFTdist-CA_BIND`, 
                        SprotFTdist$`SprotFTdist-DNA_BIND`, SprotFTdist$`SprotFTdist-NP_BIND`, SprotFTdist$`SprotFTdist-METAL`, 
                        SprotFTdist$`SprotFTdist-MOD_RES`, SprotFTdist$`SprotFTdist-CARBOHYD`, SprotFTdist$`SprotFTdist-MOTIF`, 
                        SprotFTdist$`SprotFTdist-LIPID`)
        }
      }
      
      Interface <- tmpRes$Interface
      Interface <- as.numeric(Interface$`Interface-RELACCESS-MOL`) - as.numeric(Interface$`Interface-RELACCESS`)
      
      Impact <- tmpRes$Impact
      if(is.null(Impact$`Impact-ERROR`)){
        Impact <- Impact$`Impact-BOOL`
        if(is.null(Impact)){
          Impact <- NA
        }
      }
      else {
        Impact <- Impact$`Impact-ERROR`
      }
      
      
      HBonds <- tmpRes$HBonds
      if(HBonds$`HBonds-BOOL` == 'OK'){
        HBonds <- 20
      }
      else{
        HBonds <- -20
      }
      if(tmpRes$HBonds$`HBonds-ENERGY` != "NULL"){
        HBonds <- as.numeric(tmpRes$HBonds$`HBonds-ENERGY`)
      }
      
      CisPro <- tmpRes$CisPro
      CisPro <- CisPro$`CisPro-BOOL`
      
      Clash <- tmpRes$Clash
      Clash <- as.numeric(Clash$`Clash-ENERGY`)
      
      Glycine <- tmpRes$Glycine
      if(!is.null(Glycine$`Glycine-ERROR`)){
        Glycine <- -100
      }
      else {
        if(Glycine$`Glycine-BOOL` != 'OK'){
          GlycinesNotOK <- c(GlycinesNotOK, paste(modelsFolders[i], ' -> ', variantsFiles[j], sep = ''))
        }
        if(Glycine$`Glycine-NATIVE-ENERGY` != ''){
          Glycine <- as.numeric(Glycine$`Glycine-NATIVE-ENERGY`) - as.numeric(Glycine$`Glycine-MUTANT-ENERGY`)
        }
        else{
          Glycine <- -100
        }
      }
      
      Proline <- tmpRes$Proline
      if(!is.null(Proline$`Proline-ERROR`)){
        Proline <- -100
      }
      else {
        if(Proline$`Proline-BOOL` != 'OK'){
          ProlinesNotOK <- c(ProlinesNotOK, paste(modelsFolders[i], ' -> ', variantsFiles[j], sep = ''))
        }
        if(Proline$`Proline-NATIVE-ENERGY` != ''){
          Proline <- as.numeric(Proline$`Proline-NATIVE-ENERGY`) - as.numeric(Proline$`Proline-MUTANT-ENERGY`)
        }
        else{
          Proline <- -100
        }
      }
      
      SSGeom <- as.character(tmpRes$SSGeom)
      
      SurfacePhobic <- tmpRes$SurfacePhobic
      SurfacePhobic <- as.numeric(SurfacePhobic$`SurfacePhobic-NATIVE-HPHOB`) - as.numeric(SurfacePhobic$`SurfacePhobic-MUTANT-HPHOB`)
      
      CorePhilic <- tmpRes$CorePhilic
      CorePhilic <- as.numeric(CorePhilic$`CorePhilic-NATIVE-HPHOB`) - as.numeric(CorePhilic$`CorePhilic-MUTANT-HPHOB`)
      
      BuriedCharge <- tmpRes$BuriedCharge
      BuriedCharge <- as.numeric(BuriedCharge$`BuriedCharge-NATIVE-CHARGE`) - as.numeric(BuriedCharge$`BuriedCharge-MUTANT-CHARGE`)
      
      RelAccess <- as.numeric(tmpRes$BuriedCharge$`BuriedCharge-RELACCESS`)
      
      Voids <- tmpRes$Voids
      VoidsDiff <- as.numeric(Voids$`Voids-NATIVE-LARGEST`) - as.numeric(Voids$`Voids-MUTANT-LARGEST`)
      VoidsMut <- as.numeric(Voids$`Voids-MUTANT`)
      VoidsNat <- as.numeric(Voids$`Voids-NATIVE`)
      
      if(i == 1 & j == 1){
        SAAPdapResults.dt <- as.data.table(t(c(Model, target_id, Binding, SProtFTs, Interface, RelAccess, Impact, HBonds, SurfacePhobic, CorePhilic, 
                                               BuriedCharge, SSGeom, VoidsDiff, VoidsMut, VoidsNat, Clash, Glycine, Proline, CisPro)))
      }
      else{
        tmpRow <- as.data.table(t(c(Model, target_id, Binding, SProtFTs, Interface, RelAccess, Impact, HBonds, SurfacePhobic, CorePhilic, 
                                    BuriedCharge, SSGeom, VoidsDiff, VoidsMut, VoidsNat, Clash, Glycine, Proline, CisPro)))
        SAAPdapResults.dt <- rbind(SAAPdapResults.dt, tmpRow)
      }
    }
    setwd('../')
    track_progress(i, length(modelsFolders))
  }
  setwd('../../')
  
  if(length(ExternalProgramFailed) > 1){
    ExternalProgramFailed <- ExternalProgramFailed[-1]
    cat(paste(' ', warning_sign, ' Received external program failed warning for the following models and variants ', warning_sign, '\n', sep = ''))
    cat(paste('\t', ExternalProgramFailed, sep = ''), sep = '\n')
  }
  if(length(GlycinesNotOK) > 1){
    GlycinesNotOK <- GlycinesNotOK[-1]
    cat(paste(' ', warning_sign, ' Received "GlycineBOOL Not OK" warning for the following models and variants ', warning_sign, '\n', sep = ''))
    cat(paste('\t', GlycinesNotOK, sep = ''), sep = '\n')
  }
  if(length(ProlinesNotOK) > 1){
    ProlinesNotOK <- ProlinesNotOK[-1]
    cat(paste(' ', warning_sign, ' Received "ProlineBOOL Not OK" warning for the following models and variants ', warning_sign, '\n', sep = ''))
    cat(paste('\t', ProlinesNotOK, sep = ''), sep = '\n')
  }
  
  colnames(SAAPdapResults.dt) <- c("Model", "target_id", "Binding", "SProtFT_NAMES", "SProtFT_BOOL", "SProtFT_FEATURES", "SProtFTDist_BOOL", 
                                   "SProtFT_ACT_SITE", "SProtFT_BINDING", "SProtFT_CA_BIND", "SProtFT_DNA_BIND", "SProtFT_NP_BIND", 
                                   "SProtFT_METAL", "SProtFT_MOD_RES", "SProtFT_CARBOHYD", "SProtFT_MOTIF", "SProtFT_LIPID", 
                                   "Interface", "Relaccess", "Impact",	"HBonds",	"SPhobic", "CPhilic", "BCharge", "SSGeom", "Voids",	"MLargest1",
                                   "MLargest2", "MLargest3",	"MLargest4", "MLargest5", "MLargest6", "MLargest7", "MLargest8", "MLargest9", 
                                   "MLargest10", "NLargest1", "NLargest2", "NLargest3", "NLargest4", "NLargest5", "NLargest6", "NLargest7", 
                                   "NLargest8", "NLargest9", "NLargest10",	"Clash", "Glycine", "Proline", "CisPro")
  
  SAAPdapResults.dt[SSGeom == 'OK', SSGeom := 0]
  SAAPdapResults.dt[CisPro == 'OK', CisPro := 0]
  SAAPdapResults.dt[CisPro == 'BAD', CisPro := -1]
  SAAPdapResults.dt[SProtFT_BOOL == 'OK', SProtFT_BOOL := 0]
  SAAPdapResults.dt[SProtFTDist_BOOL == 'OK', SProtFTDist_BOOL := 0]
  SAAPdapResults.dt[SProtFTDist_BOOL == 'BAD', SProtFTDist_BOOL := -1]
  SAAPdapResults.dt[Binding == 'OK', Binding := 0]
  
  # Go through each column and replace NAs if column has only one unique non-NA value
  SAAPdapResults.dt[, names(SAAPdapResults.dt) := lapply(.SD, function(x) {
    unique_vals <- unique(x[!is.na(x)])
    
    if(length(unique_vals) == 1) {
      replace(x, is.na(x), unique_vals)
    } else {
      x
    }
  })]
  # Replace NA with column means, converting character columns to numeric when possible
  SAAPdapResults.dt[, names(SAAPdapResults.dt) := lapply(.SD, function(x) {
    if(is.character(x)) {
      numeric_col <- suppressWarnings(as.numeric(x))
      if(!any(is.na(numeric_col[!is.na(x)]))) {
        x <- numeric_col
      }
    }
    
    if(is.numeric(x)) {
      replace(x, is.na(x), mean(x, na.rm = TRUE))
    } else {
      x
    }
  })]
  
  check_zero_variation <- function(col) {
    length(unique(col)) == 1
  }
  zero_variation_columns <- names(SAAPdapResults.dt)[sapply(SAAPdapResults.dt, check_zero_variation)]
  if(length(zero_variation_columns) > 0){
    cat(paste(' ', warning_sign, ' The following columns have zero variation and will be removed ', warning_sign, '\n', sep = ''))
    cat(paste('\t', zero_variation_columns, sep = ''), sep = '\n')
  }
  
  SAAPdapResults.dt[, (zero_variation_columns) := NULL]
  SAAPdapResults.dt[, Impact := NULL]
  SAAPdapResults.dt <- SAAPdapResults.dt[complete.cases(SAAPdapResults.dt)]
  
  cols_to_convert <- names(SAAPdapResults.dt)[-c(1, 2)]  # Get names of columns to convert
  SAAPdapResults.dt[, (cols_to_convert) := lapply(.SD, as.numeric), .SDcols = cols_to_convert]
  SAAPdapResults.dt <- merge(DHDomains2PDBsAnnot, SAAPdapResults.dt, by = 'Model')
  SAAPdapResults.dt <- as.data.table(SAAPdapResults.dt)
  
  if(isFALSE(use_existing_results)){
    saveRDS(SAAPdapResults.dt, file = 'new_run/saapdap_results_.rds')
    cat(paste(' Parsed SAAPdap results are saved as: ../new_run/saapdap_results_.rds ', completion_sign, '\n', sep = ''))
  }
  
  return(SAAPdapResults.dt)
}

#' Main function to run the SAAPdap pipeline
#' 
#' @param variant_data A data.table containing the variant data
#' @param variants_inherited A data.table containing the output of variant inheritance process
#' @param Prots2PDBsAnnot A data.table containing the protein to PDB mapping
#' @param aa_codes Vector to convert 1-letter amino acid codes to 3-letter amino acid codes
#' @param reverse_aa_codes Vector to convert 3-letter amino acid codes to 1-letter amino acid codes
#' @return A data.table containing the variant data with SAAPdap results
#' @export
#' 
saapdap_main <- function(use_existing_results = TRUE, SAAP3_path, variant_data, variants_inherited = NULL, Prots2PDBsAnnot = NULL, aa_codes, reverse_aa_codes = NULL, DHDomains2PDBsAnnot, check_models = FALSE){
  
  underscore_start("Starting SAAPdap pipeline")
  
  if(!use_existing_results){
    cat(paste(' ', info_sign, ' "use_existing_results = FALSE" is selected. Starting a new SAAPdap run. ', info_sign, '\n', sep = ''))
    
    if(is.null(variants_inherited) | is.null(Prots2PDBsAnnot) | is.null(reverse_aa_codes)){
      stop(paste(' ', error_sign, ' variants_inherited and Prots2PDBsAnnot must be provided if use_existing_results = FALSE', sep = ''))
    }
    if(check_models){
      check_models_for_SAAPdap()
    }
    prepare_data_for_SAAPdap(variants_inherited, Prots2PDBsAnnot, aa_codes, reverse_aa_codes)
    
    input_folder="new_run/saapdap_files/"
    input_folder <- normalizePath(input_folder)
    input_files <- list.files(input_folder, pattern = '.*\\.txt$', full.names = TRUE)
    input_files <- gsub('//', '/', input_files, fixed = T)
    output_folder_template="new_run/saapdap_res/"
    output_folder_template <- normalizePath(output_folder_template)
    SAAP3_path <- paste(SAAP3_path, 'bin/multiUniprotPipeline', sep = '')
    setwd(input_folder)
    
    for(i in 1:length(input_files)){
      input_file <- input_files[i]
      output_folder <- strsplit(input_file, '/', fixed = T)[[1]]
      output_folder <- output_folder[length(output_folder)]
      output_folder <- paste(output_folder_template, '/', output_folder, sep = '')
      output_folder <- gsub('.txt', '', output_folder, fixed = T)
      dir.create(output_folder, showWarnings = FALSE)
      system(paste(SAAP3_path, ' -v -m -f ', output_folder, ' ', input_file, sep = ''))
    }
  }
  else {
    cat(paste(' ', warning_sign, ' Using SAAPdap results from the stable run! use_existing_results = FALSE to run SAAPdap pipeline from the start ', warning_sign, '\n', sep = ''))
    if(!dir.exists('stable_run/saapdap_res/')){
      stop(paste(' ', error_sign, ' Output "saapdap_res/" folder not found! Try "use_existing_results = FALSE" to run SAAPdap pipeline from the start', sep = ''))
    }
  }
  
  SAAPdapResults <- parse_SAAPdap_results(use_existing_results, aa_codes, DHDomains2PDBsAnnot)
  variant_data <- merge(SAAPdapResults, variant_data, by = c('target_id', 'Inherited_Seq'))
  
  underscore_finish('SAAPdap pipeline finished')
  return(variant_data)
}