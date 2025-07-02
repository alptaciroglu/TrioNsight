#' Aligns the structures in the Models folder to the fixed structure and saves the aligned structures in the models_aligned folder
#'
#' @param modelsFolder the folder containing the structures to be aligned
#' @param fixed_structure_file the fixed structure to which the other structures will be aligned
#' @param max_aln_attempts the maximum number of attempts to connect to the EBI server for each alignment
#' @param max_aln_cycles the maximum number of cycles for each alignment
#' @export
#'
aln_structures <- function(modelsFolder = 'models/', fixed_structure_file = 'O75962_1284_1470_af.pdb', aligned_str_folder = 'models_aligned/', max_aln_attempts = 10, max_aln_cycles = 10){
  
  if(!dir.exists(modelsFolder)){
    stop(paste(' ', error_sign, ' "structure_analysis" directory must contain "models/" directory!', sep = ''))
  }
  setwd(modelsFolder)
  if(!file.exists(fixed_structure_file)){
    stop(paste(' ', error_sign, ' models/ folder must contain the specified fixed structure -> ' , fixed_structure_file, ' !', sep = ''))
  }
  aligned_str_folder <- paste('../', aligned_str_folder, sep = '')
  dir.create(aligned_str_folder, showWarnings = FALSE)
  
  mobile_structures <- list.files(pattern = ".pdb")
  mobile_structures <- setdiff(mobile_structures, fixed_structure_file)
  fixed_structure <- read.pdb(fixed_structure_file)
  
  cat(paste(' ', info_sign, ' Pairwise aligning ', length(mobile_structures), ' to the fixed structure -> ', fixed_structure_file, ' ', info_sign, '\n', sep = ''))
  cat(paste(' ', info_sign, ' EBI servers will be used for the alignment. Maximum of ', max_aln_attempts, ' will be made to connect to the EBI server for each alignment', info_sign, '\n', sep = ''))
  cat(paste(' ', info_sign, ' Aligned structures will be saved in the ', gsub('../', '', aligned_str_folder, fixed = T), ' folder in the working directory ', info_sign, '\n', sep = ''))
  
  for(i in 1:length(mobile_structures)){
    mobile_structure <- read.pdb(mobile_structures[i])
    cat(paste(' ', calculation_sign, ' Aligning fixed structure to ', mobile_structures[i], calculation_sign, '\n', sep = ''))
    
    attempts <- 0
    while (attempts < max_aln_attempts) {
      attempts <- attempts + 1
      tryCatch({
        
        invisible(capture.output(result <- bio3d::struct.aln(fixed = fixed_structure, mobile = mobile_structure, web.args=list(email="trionsight@gmail.com"), max.cycles = max_aln_cycles)))
        
        file.remove('aln.fa')
        setwd('fitlsq')
        if(!file.exists(paste('../', aligned_str_folder, fixed_structure_file, sep = ''))){
          file.copy('fixed.pdb', 'O75962_1284_1470_af.pdb', overwrite = FALSE)
        }
        fitted_structures <- list.files(pattern = 'mobile_')
        fitted_structures <- gsub('mobile_', '', fitted_structures, fixed = T)
        fitted_structures <- gsub('.pdb', '', fitted_structures, fixed = T)
        fitted_structures <- as.numeric(fitted_structures)
        fitted_structures <- max(fitted_structures)
        best_fit <- paste('mobile_', fitted_structures, '.pdb', sep = '')
        file.copy(best_fit, paste('../', aligned_str_folder, mobile_structures[i], sep = ''), overwrite = FALSE)
        setwd('..')
        dir.remove('fitlsq')
        
        cat(paste(' Alignment successfully completed on attempt ', attempts, ' ', completion_sign, '\n', sep = ''))
        cat(paste(' Best aligned structure is saved ', completion_sign, '\n', sep = ''))

        next        
      }, error = function(e) {
        if (attempts < max_aln_attempts) {
          cat(paste(' ', calculation_sign, ' Attempt ', attempts, ' failed. Retrying ', calculation_sign, '\n', sep = ''))
        } else {
          cat(paste(' ', warning_sign, ' All ', max_aln_attempts ,' attempts failed for this alignment. Moving on to the next alignment ', warning_sign, '\n', sep = ''))
        }
      })
    }
  }
  setwd('../')
}

#' Creates structure clustering files required by the strucclus algorithm
#' 
#' @param variant_data data.frame holding the variant data
#' @param aligned_str_folder the folder containing the PDB files
#' @param aa_codes the amino acid codes
#' @param reverse_aa_codes the reverse amino acid codes
#' @return the structure clustering files
#' @export
#' 
create_strucclus_dat <- function(variant_data, aligned_str_folder, aa_codes, reverse_aa_codes){
  
  cat(paste(' ', info_sign, ' Creating files containing coordinates to guide clustering ', info_sign, '\n', sep = ''))
  
  setwd(aligned_str_folder)
  pdb_files <- list.files(pattern = ".pdb", full.names = TRUE)
  pdb_list <- list()
  
  cat(paste(' ', calculation_sign, ' Reading pdb files in the folder ', calculation_sign, '\n', sep = ''))
  for(i in 1:length(pdb_files)){
    
    invisible(capture.output(pdb_list[[i]] <- bio3d::read.pdb(pdb_files[i], rm.alt=TRUE, verbose = FALSE)))
    names(pdb_list)[i] <- gsub(".pdb", "", basename(pdb_files[i]), fixed = T)
    
    track_progress(i, length(pdb_files))
  }
  
  variants2cluster <- variant_data[,c("target_id", "Variant_Seq", "Class")]
  variants2cluster <- unique(variants2cluster)
  colnames(variants2cluster)[3] <- 'class'
  
  var_spatial_loc_list <- list()
  for(i in 1:nrow(variants2cluster)){
    next_variant <- variants2cluster$target_id[i]
    variant_loc <- sapply(strsplit(next_variant, "_", fixed = T), function(x) x[2])
    variant_loc <- extractNumbersFromString(variant_loc)
    
    variant_Seq <- variants2cluster$Variant_Seq[i]
    variant_Seq <- paste(variant_Seq, '_af', sep = '')
    
    pdb_start <- sapply(strsplit(variant_Seq, "_", fixed = T), function(x) x[2])
    pdb_start <- as.numeric(pdb_start)
    variant_loc_onpdb <- variant_loc - pdb_start + 1
    
    pdb_variant <- pdb_list[[variant_Seq]]
    ca_atoms <- atom.select(pdb_variant, "calpha")
    ca_coords <- pdb_variant$atom[ca_atoms$atom, c("resno", "x", "y", "z")]
    ca_coords <- as.data.table(ca_coords)
    
    var_spatial_loc <- ca_coords[resno == variant_loc_onpdb]
    x <- var_spatial_loc$x
    y <- var_spatial_loc$y
    z <- var_spatial_loc$z
    var_spatial_loc_list[[next_variant]] <- c(x, y, z)
  }
  
  cat(paste(' ', calculation_sign, ' Mapping the variants to supplied structural models based on closest residues coordinates ', calculation_sign, '\n', sep = ''))
  closest_residue_list <- list()
  for(i in 1:length(pdb_list)){
    next_pdb <- pdb_list[[i]]
    closest_residue_list[[names(pdb_list)[i]]] <- list()
    for(j in 1:length(var_spatial_loc_list)){
      target_id <- names(var_spatial_loc_list)[j]
      next_var <- var_spatial_loc_list[[j]]
      closest_residue <- find_closest_residue(pdb = next_pdb, var_coords = next_var, var_id = target_id, aa_codes = aa_codes, reverse_aa_codes = reverse_aa_codes)
      if(j == 1){
        closest_residue_df <- as.data.frame(closest_residue)
      }
      else{
        closest_residue <- as.data.frame(closest_residue)
        closest_residue_df <- rbind(closest_residue_df, closest_residue)
      }
    }
    closest_residue_df$resid <- str_to_title(closest_residue_df$resid)
    closest_residue_list[[names(pdb_list)[i]]] <- as.data.table(unique(closest_residue_df))
    
    track_progress(i, length(pdb_list))
  }
  
  closest_residue_list <- lapply(closest_residue_list, function(x){convert_aa_codes(x, aa_codes, reverse_aa_codes, column_name = "resid")})
  closest_residue_list <- lapply(closest_residue_list, function(x){x[,resid := NULL]; x[,resid_3letter := NULL]; x[, resno := as.numeric(resno)]; colnames(x)[4] <- 'resid'; x <- unique(x); x[,.(chain, resno, resid, target_id)]})
  
  cat(paste(' ', calculation_sign, ' Writing mapped variant files to .dat files ', calculation_sign, '\n', sep = ''))
  for(i in 1:length(closest_residue_list)){
    fileName <- paste(names(closest_residue_list)[i], '.dat', sep = '')
    fileName <- gsub(pattern = '//', replacement = '/', fileName, fixed = T)
    closest_residue_list[[i]] <- merge(closest_residue_list[[i]], variants2cluster[,c("target_id", 'class')], by = 'target_id')
    closest_residue_list[[i]] <- closest_residue_list[[i]][,.(chain, resno, resid, target_id, class)]
    fwrite(closest_residue_list[[i]], file = fileName, col.names = F, quote = F, sep = ' ')
  }
  setwd('../')
  return(closest_residue_list)
}

#' Parse a single clustering guide file created by create_strucclus_dat function
#' 
#' @param file_path Path to the clustering guide file
#' @return A data.frame of clusters, each containing the cluster number, coordinates, residue numbers, PDB indices, and target IDs
#' @export
#' 
parse_single_cluster_file <- function(file_path) {
  
  # Read the file
  lines <- readLines(file_path)
  
  # Find the index where "Cluster 1" starts
  start_index <- which(str_detect(lines, "^Cluster 1"))[1]
  
  # Extract only the relevant lines
  data_lines <- lines[start_index:length(lines)]
  
  # Initialize variables
  clusters <- list()
  current_cluster <- NULL
  
  # Parse the lines
  for (line in data_lines) {
    if (str_detect(line, "^Cluster")) {
      if (!is.null(current_cluster)) {
        clusters[[length(clusters) + 1]] <- current_cluster
      }
      current_cluster <- list(
        cluster_num = as.integer(str_extract(line, "\\d+")),
        data = data.frame(Coordinates = character(), Residue = integer(), PDB_Idx = character(), target_id = character(), Class = character(), stringsAsFactors = FALSE)
      )
    } else if (str_detect(line, "#")) {
      parts <- str_split(line, "#")[[1]]
      coordinates <- str_trim(parts[1])
      info <- str_trim(parts[2])
      
      residue <- as.integer(str_extract(info, "\\d+"))
      pdb_idx <- sapply(strsplit(info, ' ', fixed = T), function(x){paste(x[1], x[2], x[3], sep = ' ')})
      target_id <- sapply(strsplit(info, ' ', fixed = T), function(x){x[length(x)-1]})
      variant_class <- sapply(strsplit(info, ' ', fixed = T), function(x){x[length(x)]})
      
      current_cluster$data <- rbind(current_cluster$data, 
                                    data.frame(Coordinates = coordinates, Residue = residue, PDB_Idx = pdb_idx, target_id = target_id, Class = variant_class, stringsAsFactors = FALSE))
    }
  }
  # Add the last cluster
  if (!is.null(current_cluster)) {
    clusters[[length(clusters) + 1]] <- current_cluster
  }
  
  # Combine all clusters into a single dataframe
  result <- do.call(rbind, lapply(clusters, function(cluster) {
    data.frame(
      Residue = cluster$data$Residue,
      Coordinates = cluster$data$Coordinates,
      Cluster = cluster$cluster_num,
      PDB_Idx = cluster$data$PDB_Idx,
      target_id = cluster$data$target_id,
      Class = cluster$data$Class,
      stringsAsFactors = FALSE
    )
  }))
  return(result)
}

#' Process all cluster files in a folder calling parse_single_cluster_file on each
#' 
#' @param folder_path Path to the folder containing the cluster files
#' @return A list of data.frames, each containing the clusters from a single file
#' @export
#' 
process_all_cluster_files <- function(folder_path) {
  
  # Check if the folder exists
  if (!dir.exists(folder_path)) {
    stop(paste(' ', error_sign, ' The specified folder does not exist: ', folder_path, sep = ''))
  }
  
  cat(paste(' ', calculation_sign, ' Processing strucclus results ', calculation_sign, '\n', sep = ''))
  
  # Get all txt files in the folder
  all_files <- list.files(path = folder_path, pattern = "*.txt", full.names = TRUE)
  
  # Check if there are any files
  if (length(all_files) == 0) {
    stop(paste(' ', error_sign, ' No .txt files found in the specified folder: ', folder_path, sep = ''))
  }
  
  # Get all txt files in the folder
  all_files <- list.files(path = folder_path, pattern = "*.txt", full.names = TRUE)
  
  # Extract base names (without n1, n2, etc.)
  base_names <- unique(str_replace(basename(all_files), "_n\\d+\\.txt$", ""))
  
  # Initialize the result list
  result_list <- list()
  
  # Process files for each base name
  i <- 1
  folder_count <- length(base_names)
  for (base_name in base_names) {
    # Get all files for this base name
    base_files <- all_files[str_detect(all_files, paste0(base_name, "_n\\d+\\.txt$"))]
    
    # Sort the files numerically by the n value
    base_files <- base_files[order(as.numeric(str_extract(base_files, "(?<=_n)\\d+(?=\\.txt$)")))]
    
    # Process each file and store results
    base_results <- list()
    for (file in base_files) {
      n_value <- str_extract(file, "(?<=_n)\\d+(?=\\.txt$)")
      base_results[[paste0("n", n_value)]] <- parse_single_cluster_file(file)
    }
    
    # Add to the main result list
    result_list[[base_name]] <- base_results
    track_progress(i, folder_count)
    i <- i + 1
  }
  return(result_list)
}

#' Function to perform chi-square test and return p-value and statistic for given data.frame
#' 
#' @param df A data.frame containing the cluster results
#' @return A list containing the p-value and statistic of the chi-square test
#' @export
#' 
chi_square_test_single <- function(df) {
  # Check if there's more than one cluster
  if (length(unique(df$Cluster)) < 2) {
    return(NA)
  }
  
  contingency_table <- table(df$Cluster, df$Class)
  
  # Check if any expected frequencies are less than 5
  chi_test <- suppressWarnings(chisq.test(contingency_table))
  return(list(chi_test$p.value, chi_test$statistic))
}

#' Function to determine the best n for each cluster based on chi-square test
#' 
#' @param base_results A list of data.frames containing the cluster results for each n
#' @return A list containing p-values and chi-square statistics for each n
#' @export
#' 
chi_square_test_multi <- function(base_results) {
  chisq_res <- sapply(base_results, function(df) {
    tryCatch({
      chi_square_test_single(df)
    }, error = function(e) {
      cat(paste(' ', warning_sign, ' Error in chi-square test: ', e$message, ' ', warning_sign, '\n', sep = ''))
      return(NA)
    })
  })
  
  p_values <- sapply(chisq_res, function(x){x[[1]]})
  x_squared_vals <- sapply(chisq_res, function(x){if(length(x) > 1){x[[2]]} else{NA}})
  names(x_squared_vals) <- names(p_values)
  
  # Remove NA values
  p_values <- p_values[!is.na(p_values)]
  x_squared_vals <- x_squared_vals[!is.na(x_squared_vals)]
  
  if (length(p_values) == 0) {
    cat(paste(' ', warning_sign, ' Valid p-values could not be calculated ', warning_sign, '\n', sep = ''))
    return(list(best_n = NA, p_values = NA, x_squared = NA))
  }
  
  best_n <- names(which.min(p_values))
  return(list(best_n = best_n, p_values = p_values, x_squared = x_squared_vals))
}

#' Function to process all results and determine best n for each base name
#' 
#' @param results A list of data.frames containing the cluster results for each n
#' @return A data.frame containing the best n results corresponding metrics for each n
#' @export
#' 
determine_best_n_results <- function(results) {
  
  if (!is.list(results) || length(results) == 0) {
    stop(paste(' ', error_sign, ' Invalid input: results should be a non-empty list ', sep = ''))
  }
  cat(paste(' ', calculation_sign, ' Analysing best n = number of clusters ', calculation_sign, '\n', sep = ''))
  cat(paste(' ', calculation_sign, ' Running chi-square test on clusters variant composition ', calculation_sign, '\n', sep = ''))
  
  cat(paste(' ', warning_sign, ' Cannot perform chi-square test on results with fewer than 2 clusters ', warning_sign, '\n', sep = ''))
  cat(paste(' ', warning_sign, ' Chi-square approximation may be incorrect in some results due to low expected frequencies in determining the best clustering ', warning_sign, '\n', sep = ''))
  
  best_n_results <- lapply(results, chi_square_test_multi)
  
  # Create a summary dataframe
  best_n_summarised <- do.call(rbind, lapply(names(best_n_results), function(base_name) {
    data.frame(
      BaseName = base_name,
      BestN = best_n_results[[base_name]]$best_n,
      PValue = if(is.na(best_n_results[[base_name]]$best_n)) NA else best_n_results[[base_name]]$p_values[best_n_results[[base_name]]$best_n],
      chi_squared = if(is.na(best_n_results[[base_name]]$best_n)) NA else best_n_results[[base_name]]$x_squared[best_n_results[[base_name]]$best_n]
    )
  }))
  
  all_columns <- unique(unlist(lapply(best_n_results, function(x){names(x$p_value)})))
  # Create a data frame with all unique columns, filling missing values with NA
  best_n_list <- lapply(best_n_results, function(x) {
    tmp <- setNames(rep(NA, length(all_columns)), all_columns)
    tmp[names(x$p_value)] <- x$p_value
    return(tmp)
  })
  
  best_n_all_df <- do.call(rbind, best_n_list)
  best_n_all_df <- as.data.frame(best_n_all_df)
  best_n_all_df$id <- rownames(best_n_all_df)
  fwrite(best_n_all_df, file = "../new_run/best_n_all.txt", sep = '\t')
  
  rownames(best_n_summarised) <- 1:nrow(best_n_summarised)
  best_n_summarised <- as.data.table(best_n_summarised)
  fwrite(best_n_summarised, file = "../new_run/best_n_summary.txt", sep = '\t')
  cat(paste(' Results are written to "new_run/best_n_all.txt" & "new_run/best_n_summary.txt ', completion_sign, '\n', sep = ''))
  if(nrow(best_n_summarised[best_n_summarised$PValue > 0.05,]) > 0){
    cat(paste(' ', warning_sign, ' Some strucclus clusterings resulted in p-value > 0.05 ', warning_sign, '\n', sep = ''))
  }
  cat(paste(' Completed ', completion_sign, '\n', sep = ''))
  return(best_n_summarised)
}

#' Calculate distance metric - contains the mathematical formula only
#' 
#' @param residue object containing coordinates of the residue to measure distance
#' @param best_clustering data.frame containing the best clustering results
#' @return A data.frame containing the distance of the residue to each cluster
#' @export
#' 
calculate_distance_metric_formula <- function(residue, best_clustering) {
  
  # Extract coordinates and convert to numeric
  residue_coords <- as.numeric(str_split(residue$Coordinates, " ")[[1]])
  
  # Function to calculate Euclidean distance
  euclidean_distance <- function(coord1, coord2) {
    sqrt(sum((coord1 - coord2)^2))
  }
  
  # Calculate cluster centers, and composition
  clusters <- best_clustering %>% 
    group_by(Cluster) %>% 
    summarise(
      center_x = mean(as.numeric(sapply(str_split(Coordinates, " "), `[`, 1))),
      center_y = mean(as.numeric(sapply(str_split(Coordinates, " "), `[`, 2))),
      center_z = mean(as.numeric(sapply(str_split(Coordinates, " "), `[`, 3))),
      .groups = 'drop')
  
  # Calculate distances to all cluster centers
  cluster_distances <- clusters %>% 
    rowwise() %>%
    mutate(distance = euclidean_distance(residue_coords, c(center_x, center_y, center_z))) %>%
    ungroup()
  
  return(cluster_distances$distance)
}

#' Function to process the PDB file to calculate distance metrics
#' 
#' @param pdb_file Path to the PDB file
#' @param best_clustering data.frame containing the best clustering results
#' @return A data.frame containing the distance of each residue to each cluster
#' @export
#' 
calculate_distance_metrics_onepdb <- function(pdb_file, best_clustering) {
  
  # Read PDB file and extract residue information
  # This is a simplified version and may need to be adapted based on your PDB file format
  invisible(capture.output(pdb_data <- bio3d::read.pdb(pdb_file, rm.alt=TRUE, verbose = FALSE)))
  pdb_data <- pdb_data$atom
  
  # Process each residue
  residues <- pdb_data %>% 
    filter(type == "ATOM") %>%
    filter(elety == "CA") %>%
    group_by(resno, chain) %>%
    summarise(Coordinates = paste(mean(x), mean(y), mean(z), collapse = " "), .groups = 'drop') %>%
    ungroup()
  
  # Calculate distance metric for each residue
  distance_to_clusters <- sapply(1:nrow(residues), function(i) {
    calculate_distance_metric_formula(residues[i,], best_clustering)
  })
  distance_to_clusters <- t(distance_to_clusters)
  colnames(distance_to_clusters) <- paste0("distance_to_cluster_", 1:ncol(distance_to_clusters))
  residues <- cbind(residues, distance_to_clusters)
  
  return(residues)
}

#' Function to calculate distance metrics for all variants on all structures from scratch
#' 
#' @param variant_data_annotated data.frame containing the annotated variant data
#' @param best_n_analysis data.frame containing the best N analysis results
#' @param aligned_str_folder Path to the folder containing aligned structures
#' @param clustering_parsed data.frame containing the parsed clustering results obtained by running process_all_cluster_files function
#' @param distance_reference_pdb Name of the reference PDB file to calculate distance metrics
#' @param best_n_reference Reference clustering count found to be best to calculate distance metrics
#' @return A list containing distance metrics for all variants
#' @export
#' 
calculate_distance_metrics_fromstart <- function(variant_data_annotated, best_n_analysis, aligned_str_folder, clustering_parsed, distance_reference_pdb, best_n_reference){
  
  cat(paste(' ', calculation_sign, ' Calculating distance to variant clusters for each cluster on every structure ', calculation_sign, '\n', sep = ''))
  setwd(aligned_str_folder)
  distance_to_clusters_list <- list()
  for(i in 1:nrow(best_n_analysis)){
    
    inherited_model <- best_n_analysis$BaseName[i]
    best_n_inherited <- best_n_analysis$BestN[i] ## Just to pick the coordinates. Otherwise, we are using reference pdb clusters to calculate the distances
    
    if(is.null(variant_data_annotated[Model == inherited_model])){
      cat(paste(' ', warning_sign, ' No variants were inherited for the model ', inherited_model, ' ', warning_sign, '\n', sep = ''))
      cat(paste(' ', warning_sign, ' This is likely due to native residue mismatch between the model and the variant data ', warning_sign, '\n', sep = ''))
      cat(paste(' ', warning_sign, ' Skipping model ', inherited_model, ' ', warning_sign, '\n', sep = ''))
      next
    }
    
    pdb_file_next <- paste(best_n_analysis$BaseName[i], ".pdb", sep = "")
    pdb_file_next <- gsub('//' , '/', pdb_file_next, fixed = T)
    
    variant_coords <- clustering_parsed[[inherited_model]][[best_n_inherited]]
    variant_coords$Model <- inherited_model
    variant_coords$Cluster <- NULL
    variant_coords$Residue <- NULL
    
    best_clustering_next <- clustering_parsed[[distance_reference_pdb]][[best_n_reference]]
    best_clustering_next$Model <- distance_reference_pdb
    
    distance_to_clusters_metrics <- calculate_distance_metrics_onepdb(pdb_file = pdb_file_next, best_clustering = best_clustering_next)
    distance_to_clusters_metrics <- merge(variant_coords, distance_to_clusters_metrics, by = 'Coordinates')
    distance_to_clusters_metrics <- merge(distance_to_clusters_metrics, unique(variant_data_annotated[Model == inherited_model][,.(target_id)]), by = 'target_id')
    
    if(nrow(unique(variant_data_annotated[,.(target_id, Model)][Model == inherited_model])) != nrow(distance_to_clusters_metrics)){
      stop(paste(' ', error_sign, ' The number of variants in the variant annotation list does not match the number of variants in the distance to clusters metrics ', sep = ''))
    }
    
    distance_colnames <- extractNumbersFromString(best_n_reference)
    distance_colnames <- seq(1, as.numeric(distance_colnames))
    distance_colnames <- paste0('distance_to_cluster_', distance_colnames)
    distance_colnames <- c('target_id', 'Model', distance_colnames)
    
    distance_to_clusters_metrics <- distance_to_clusters_metrics[,..distance_colnames]
    distance_to_clusters_list[[i]] <- distance_to_clusters_metrics
    track_progress(i, nrow(best_n_analysis))
  }
  setwd('../')
  return(distance_to_clusters_list)
  
}
  
#' Function to calculate distance metrics for all variants on all structures from known clusters
#' 
#' @param variant_data_annotated data.frame containing the annotated variant data
#' @param best_n_analysis data.frame containing the best N analysis results
#' @param aligned_str_folder Path to the folder containing aligned structures
#' @param reference_clustering a data.frame containing the reference clustering results
#' @param model_domainboundaries a data.frame containing the domain boundaries for each model
#' @param distance_reference_pdb Name of the reference PDB file to calculate distance metrics
#' @param best_n_reference Reference clustering count found to be best to calculate distance metrics
#' @return A list containing distance metrics for all variants
#' @export
#' 
calculate_distance_metrics_knownclusters <- function(variant_data_annotated, best_n_analysis, aligned_str_folder, reference_clustering, model_domainboundaries, distance_reference_pdb = 'TrioN_DHPH_1290_1594_A_1nty', best_n_reference = 'n5'){
  
  cat(paste(' ', calculation_sign, ' Calculating distance to variant clusters for known clusters ', calculation_sign, '\n', sep = ''))
  setwd(aligned_str_folder)
  distance_to_clusters_list <- list()
  
  for(i in 1:nrow(best_n_analysis)){
    
    inherited_model <- best_n_analysis$BaseName[i]
    best_n_inherited <- best_n_analysis$BestN[i] ## Just to pick the coordinates. Otherwise, we are using reference pdb clusters to calculate the distances
    
    if(is.null(variant_data_annotated[Model == inherited_model])){
      cat(paste(' ', warning_sign, ' No variants were inherited for the model ', inherited_model, ' ', warning_sign, '\n', sep = ''))
      cat(paste(' ', warning_sign, ' This is likely due to native residue mismatch between the model and the variant data ', warning_sign, '\n', sep = ''))
      cat(paste(' ', warning_sign, ' Skipping model ', inherited_model, ' ', warning_sign, '\n', sep = ''))
      next
    }
    
    pdb_file_next <- paste(best_n_analysis$BaseName[i], ".pdb", sep = "")
    pdb_file_next <- gsub('//' , '/', pdb_file_next, fixed = T)
    
    invisible(capture.output(variant_coords <- bio3d::read.pdb(pdb_file_next, rm.alt=TRUE, verbose = FALSE)))
    variant_coords <- variant_coords$atom
    # Process each residue
    variant_coords <- variant_coords %>% 
      filter(type == "ATOM") %>%
      filter(elety == "CA") %>%
      group_by(resno, chain, resid) %>%
      summarise(Coordinates = paste(mean(x), mean(y), mean(z), collapse = " "), .groups = 'drop') %>%
      ungroup()
    
    variant_coords <- as.data.table(variant_coords)
    variant_coords$resid <- str_to_title(variant_coords$resid)
    variant_coords$Model <- inherited_model
    
    variant_coords <- convert_aa_codes(variant_coords, aa_codes = aa_codes, reverse_aa_codes = reverse_aa_codes, column_name = 'resid')
    variant_coords[, PDB_Idx := paste(chain, resno, resid_1letter, sep = ' ')]
    variant_coords[, resid := NULL]
    variant_coords[, resid_3letter := NULL]
    variant_coords[, chain := NULL]
    variant_coords$Class <- '-'
    colnames(variant_coords)[1] <- 'Residue'
    colnames(variant_coords)[4] <- 'resid'
    
    substitutions <- generate_substitutions_from_df(df = variant_coords, position_field = 'Residue', residue_field = 'resid')
    
    the_model <- model_domainboundaries[Model == unique(variant_coords$Model)]
    accession <- the_model$Accession
    
    pos_start <- as.numeric(the_model$pos_start)
    pos_end <- pos_start + nrow(variant_coords) - 1
    positions <- sort(rep(seq(pos_start, pos_end, 1), 19))
    
    substitutions[, accession := accession]
    substitutions[, position := positions]
    substitutions[, target_id := paste(accession, '_', original, position, substitution, sep = '')]
    substitutions[, original := NULL]
    substitutions[, substitution := NULL]
    substitutions[, accession := NULL]
    colnames(substitutions)[1] <- 'Residue'
    variant_coords <- merge(variant_coords, substitutions, 'Residue')
    
    variant_coords[, Residue := NULL]
    variant_coords[, resid := NULL]
    variant_coords <- variant_coords[position >= as.numeric(as.character(the_model$domain_start)),]
    variant_coords <- variant_coords[position <= as.numeric(as.character(the_model$domain_end)),]
    variant_coords <- variant_coords[,.(Coordinates, PDB_Idx, target_id, Class, Model)]
    
    distance_to_clusters_metrics <- calculate_distance_metrics_onepdb(pdb_file = pdb_file_next, best_clustering = reference_clustering)
    distance_to_clusters_metrics <- merge(variant_coords, distance_to_clusters_metrics, by = 'Coordinates')
    distance_to_clusters_metrics <- merge(distance_to_clusters_metrics, unique(variant_data_annotated[Model == inherited_model][,.(target_id)]), by = 'target_id')
    
    if(nrow(unique(variant_data_annotated[,.(target_id, Model)][Model == inherited_model])) != nrow(distance_to_clusters_metrics)){
      stop(paste(' ', error_sign, ' The number of variants in the variant annotation list does not match the number of variants in the distance to clusters metrics ', sep = ''))
    }
    
    distance_colnames <- extractNumbersFromString(best_n_reference)
    distance_colnames <- seq(1, as.numeric(distance_colnames))
    distance_colnames <- paste0('distance_to_cluster_', distance_colnames)
    distance_colnames <- c('target_id', 'Model', distance_colnames)
    
    distance_to_clusters_metrics <- distance_to_clusters_metrics[,..distance_colnames]
    distance_to_clusters_list[[i]] <- distance_to_clusters_metrics
    track_progress(i, nrow(best_n_analysis))
  }
  setwd('../')
  return(distance_to_clusters_list)
  
}
  
#' Calculates the closest residue to the variant in the PDB structure
#' Variants that are reported in a structure and variants that are inherited to a structure are both covered
#' 
#' @param pdb the PDB structure
#' @param var_coords the coordinates of the variant
#' @param var_id the ID of the variant
#' @param aa_codes the amino acid codes 
#' @param reverse_aa_codes the reverse amino acid codes
#' @return the closest residue to the variant based on input pdb file
#' @export
#' 
find_closest_residue <- function(pdb, var_coords, var_id, aa_codes, reverse_aa_codes) {
  closest_residue <- NULL
  min_distance <- Inf
  
  # Extract CA atoms
  ca_atoms <- atom.select(pdb, "calpha")
  ca_coords <- pdb$atom[ca_atoms$atom, c("x", "y", "z")]
  
  # Calculate distances
  ca_coords$x <- ca_coords$x - var_coords[1]
  ca_coords$y <- ca_coords$y - var_coords[2]
  ca_coords$z <- ca_coords$z - var_coords[3]
  distances <- sqrt(rowSums(ca_coords^2))
  
  # Find the minimum distance and corresponding residue
  min_index <- which.min(distances)
  min_distance <- distances[min_index]
  closest_residue <- pdb$atom[ca_atoms$atom,][min_index, ]
  
  chain <- closest_residue$chain
  resno <- as.numeric(closest_residue$resno)
  resid <- closest_residue$resid
  closest_residue <- t(as.data.frame(c(chain, resno, resid, var_id)))
  colnames(closest_residue) <- c("chain", "resno", "resid", "target_id")
  rownames(closest_residue) <- NULL
  
  return(closest_residue)
}

#' Main function to calculate distance metrics for all PDB files
#' 
#' @param use_existing_results Boolean value for whether to re-use existing results as whole
#' @param use_existing_clusters_only Boolean value for whether to calculate new distances based on existing clusters
#' @param modelsFolder path to the folder containing the structural models
#' @param fixed_structure_file path to the fixed structure file - used for structural alignment
#' @param aligned_str_folder path to the folder containing the aligned structural models
#' @param max_aln_attempts maximum number of attempts for structural alignment
#' @param max_aln_cycles maximum number of cycles for structural alignment
#' @param aa_codes the amino acid codes
#' @param reverse_aa_codes the reverse amino acid codes
#' @param strucclus_res_folderpath path to the folder containing the results of the strucclus clustering
#' @param distance_reference_pdb the reference PDB file for distance calculation
#' @param best_n_reference the best N reference value - previously calculated from the strucclus clustering
#' @param variant_data_annotated the annotated variant data
#' @return a data.table containing the variant data annotated with the distance metrics
#' @export
#' 
variant_clusters_main <- function(use_existing_results = FALSE, use_existing_clusters_only = TRUE, modelsFolder, fixed_structure_file, aligned_str_folder, reference_clustering, model_domainboundaries, max_aln_attempts, max_aln_cycles, aa_codes, reverse_aa_codes, strucclus_res_folderpath, distance_reference_pdb = 'TrioN_DHPH_1290_1594_A_1nty', best_n_reference, variant_data_annotated, data_overlap = NULL){
  
  underscore_start("Starting variant clustering analysis")
  
  if(!dir.exists('structure_analysis')){
    stop(paste(' ', error_sign, ' Working directory must contain "structure_analysis/" directory!', sep = ''))
  }
  if(!dir.exists('structure_analysis/strucclus_V1.0/')){
    stop(paste(' ', error_sign, ' Working directory must contain "structure_analysis/strucclus_V1.0/" directory!', sep = ''))
  }
  if(isFALSE(use_existing_results)){
    setwd('structure_analysis/')
    aln_structures(modelsFolder = modelsFolder, fixed_structure_file = fixed_structure_file, 
                   aligned_str_folder = aligned_str_folder, max_aln_attempts = max_aln_attempts, max_aln_cycles = max_aln_cycles)
    
    if(!is.null(data_overlap)){
      cat(paste(' ', warning_sign, ' Removing variants that overlap with TrioNsight data ', warning_sign, '\n', sep = ''))
      variant_data_annotated <- as.data.table(variant_data_annotated)
      number_of_variants_before <- nrow(variant_data_annotated)
      variant_data_annotated <- variant_data_annotated[!target_id %in% data_overlap$TrioNsight]
      number_of_variants_after <- nrow(variant_data_annotated)
      cat(paste(' ', warning_sign, ' Removed ', number_of_variants_before - number_of_variants_after, ' variants that overlap with TrioNsight data ', warning_sign, '\n', sep = ''))
    }
    
    closest_residue_list <- create_strucclus_dat(variant_data = variant_data_annotated, 
                                                 aligned_str_folder = aligned_str_folder, 
                                                 aa_codes = aa_codes, reverse_aa_codes = reverse_aa_codes)
    
    setwd('strucclus_V1.0/')
    cat(paste(' ', calculation_sign, ' Running strucclus algorithm ', calculation_sign, '\n', sep = ''))
    system("chmod +x strucclus_multi.sh")
    system('./strucclus_multi.sh')
    cat(paste(' strucclus run completed ', completion_sign, '\n', sep = ''))
    setwd('../')
    
    clustering_parsed <- process_all_cluster_files(folder_path = strucclus_res_folderpath)
    best_n_analysis <- determine_best_n_results(clustering_parsed)
    distance_to_clusters <- calculate_distance_metrics_fromstart(variant_data_annotated, best_n_analysis, aligned_str_folder, clustering_parsed, distance_reference_pdb, best_n_reference)
    distance_to_clusters <- do.call(rbind, distance_to_clusters)
    fwrite(distance_to_clusters, '../new_run/distance_to_clusters_results.txt', sep = '\t')
    setwd('../')
  }
  else if(isFALSE(use_existing_results) & isTRUE(use_existing_clusters_only)){
    best_n_analysis <- fread('stable_run/best_n_summary.txt', sep = '\t', header = T)
    setwd('structure_analysis/')
    distance_to_clusters <- calculate_distance_metrics_knownclusters(variant_data_annotated, best_n_analysis, aligned_str_folder, reference_clustering, model_domainboundaries, distance_reference_pdb = 'TrioN_DHPH_1290_1594_A_1nty', best_n_reference = 'n5')
    distance_to_clusters <- do.call(rbind, distance_to_clusters)
    fwrite(distance_to_clusters, '../new_run/distance_to_clusters_results.txt', sep = '\t')
    setwd('../')
  }
  else if(isTRUE(use_existing_results)){
    distance_to_clusters <- fread('stable_run/distance_to_clusters_results.txt', sep = '\t')
  }
  variant_data_annotated <- merge(variant_data_annotated, distance_to_clusters, by = c('target_id', 'Model'))
  
  underscore_finish('Variant clustering analysis finished')
  
  return(variant_data_annotated)
}

