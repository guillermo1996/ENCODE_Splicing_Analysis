#' Title
#'
#' @param metadata
#' @param main_samples_path
#' @param annotated_SR_details
#' @param all_reads_combined
#' @param num_cores Number of multiprocessing cores to use. Memory requirements significantly increase with the number of cores.
#' @param rw_disk Whether to store the results in disk. By default, TRUE.
#' @param overwrite Whether to overwrite previously generated results from the
#'   function. If set to FALSE and 'rw_disk' is set to TRUE, the function looks
#'   for the files in memory and loads them if possible. By default, FALSE.
#'
#' @return
#' @export
juctionPairing <- function(metadata,
                           main_samples_path,
                           annotated_SR_details,
                           all_reads_combined,
                           num_cores = 4,
                           rw_disk = T,
                           overwrite = F){
  logger::log_info("\t Obtaining the raw distances for all samples (num_cores = ", num_cores, ").")
  
  ## If rw_disk argument is set to TRUE and overwrite is set to FALSE, it looks
  ## for the results of the function in disk. If they have already been
  ## generated, the function is not executed and the files are returned.
  if (rw_disk & !overwrite & file.exists(paste0(main_samples_path, "all_distances_raw.rds"))) {
    logger::log_info("\t\t Ignoring the pairing process.") 
    logger::log_info("\t\t File ", main_samples_path, "all_distances_raw.rds already exists.")
    logger::log_info("\t\t Loading file from disk.")
    return(readRDS(paste0(main_samples_path, "all_distances_raw.rds")))
  }
  
  ## Multiprocessing generation. Add argument "output.file" to
  ## parallel::makeCluster() if you want to output information about the
  ## parallel execution. 
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  
  ## ######### Calculation of the distances ############
  ##
  ## For every sample, we look for overlaps between the novel junctions and the
  ## annotated (or reference) junctions.
  ##
  ## 1. Divide the junctions in annotated, donor and acceptor and their strand
  ## (reverse or forward).
  ## 2. Look for overlaps between the novel junctions in both directions and the
  ## annotated junctions. A total of 4 comparison are made (donor and acceptor,
  ## either forward or reverse).
  ## 3. Join the results for each sample into a single dataframe.
  
  all_distances_raw <- foreach(i = 1:nrow(metadata), .packages = c("GenomicRanges", "tidyverse"), .export = "getDistancesDataFrame") %dopar%{
    row <- metadata[i, ]
    sample <- row$sample_id
    if(!(sample %in% names(all_reads_combined))) return(tibble::tibble())
    
    sample_reads <- all_reads_combined %>% 
      dplyr::select(junID, all_of(sample)) %>%
      tidyr::drop_na()
    
    sample_annotated_SR_details <- annotated_SR_details %>%
      dplyr::inner_join(sample_reads,
                        by = c("junID"),
                        sort = T) %>%
      dplyr::rename(counts = all_of(sample)) %>%
      GenomicRanges::GRanges()
    
    ## Annotated junctions
    all_annotated <- sample_annotated_SR_details[sample_annotated_SR_details$type == "annotated"]
    all_annotated_forward <- all_annotated[all_annotated@strand == "+"]
    all_annotated_reverse <- all_annotated[all_annotated@strand == "-"]
    
    ## Novel donors
    all_donor <- sample_annotated_SR_details[sample_annotated_SR_details$type == "novel_donor"]
    all_donor_forward <- all_donor[all_donor@strand == "+"]
    all_donor_reverse <- all_donor[all_donor@strand == "-"]
    
    ## Novel acceptors
    all_acceptor <- sample_annotated_SR_details[sample_annotated_SR_details$type == "novel_acceptor"]
    all_acceptor_forward <- all_acceptor[all_acceptor@strand == "+"]
    all_acceptor_reverse <- all_acceptor[all_acceptor@strand == "-"]
    
    ## Generation of the distance dataframes
    df_nd_f <- getDistancesDataFrame(all_donor_forward, all_annotated_forward, "end", sample, "novel_donor")
    df_nd_r <- getDistancesDataFrame(all_donor_reverse, all_annotated_reverse, "start", sample, "novel_donor")
    df_na_f <- getDistancesDataFrame(all_acceptor_forward, all_annotated_forward, "start", sample, "novel_acceptor")
    df_na_r <- getDistancesDataFrame(all_acceptor_reverse, all_annotated_reverse, "end", sample, "novel_acceptor")
    
    sample_distances <- rbind(df_nd_f, df_nd_r, df_na_f, df_na_r)
    return(sample_distances)
  }
  parallel::stopCluster(cl)
  
  ## Combine the resulted list into a single dataframe
  all_distances_raw <- all_distances_raw %>%
    dplyr::bind_rows() %>% 
    dplyr::distinct()
  invisible(gc())
  
  ## Save the output
  ## Recommended name: all_distances_raw.rds
  if(rw_disk){
    logger::log_info("\t\t Saving output to ", main_samples_path, "all_distances_raw.rds")
    all_distances_raw %>% saveRDS(paste0(main_samples_path, "all_distances_raw.rds"))
  }
  
  logger::log_info("\t\t Distance calculation executed successfully.")
  return(all_distances_raw)
}

#' Generates the distance dataframe for a particular novel junction type and
#' direction
#'
#' More information about the parameters in \link[GenomicRanges]{findOverlaps}.
#'
#' @param query GRanges object with junctions to find the overlaps.
#' @param subject GRanges object with junctions to find the overlaps.
#' @param type The type of overlap to look for.
#' @param sample The sample ID being studied.
#' @param junc_type The type of junction being studied.
#'
#' @return Dataframe with the novel junctions and the associated reference
#'   junctions for a given sample and type of novel donor and direction.
#' @export
getDistancesDataFrame <- function(query,
                                  subject,
                                  type,
                                  sample,
                                  junc_type) {
  ## Find the overlaps between the query and the subject GRanges
  overlaps <- GenomicRanges::findOverlaps(
    query = query,
    subject = subject,
    ignore.strand = FALSE,
    type = type
  )
  novel_junctions <- query[S4Vectors::queryHits(overlaps), ]
  ref_junctions <- subject[S4Vectors::subjectHits(overlaps), ]
  
  ## Calculation of the distances according to the direction of the junctions
  if (type == "start") {
    distance <- end(novel_junctions) - end(ref_junctions)
  } else {
    distance <- start(ref_junctions) - start(novel_junctions)
  }
  
  ## Dataframe with all the relevant information
  df <- tibble::tibble(
    sample = sample,
    type = junc_type,
    distance = distance,
    novel_junID = novel_junctions$junID,
    novel_counts = novel_junctions$counts,
    novel_seq = novel_junctions %>% seqnames() %>% as.character(),
    novel_start = novel_junctions %>% start(),
    novel_end = novel_junctions %>% end(),
    novel_strand = novel_junctions %>% strand() %>% as.character(),
    novel_width = novel_junctions %>% width(),
    #novel_ss5score = novel_junctions$ss5score,
    #novel_ss3score = novel_junctions$ss3score,
    ref_junID = ref_junctions$junID,
    ref_counts = ref_junctions$counts,
    ref_seq = ref_junctions %>% seqnames() %>% as.character(),
    ref_start = ref_junctions %>% start(),
    ref_end = ref_junctions %>% end(),
    ref_strand = ref_junctions %>% strand() %>% as.character(),
    ref_width = ref_junctions %>% width(),
    #ref_ss5score = ref_junctions$ss5score,
    #ref_ss3score = ref_junctions$ss3score
  )
  
  return(df)
}

#' Title
#'
#' @param all_distances_raw 
#' @param main_samples_path 
#' @param rw_disk 
#' @param overwrite 
#'
#' @return
#' @export
removeAmbiguousPairing <- function(all_distances_raw,
                                   main_samples_path,
                                   rw_disk = T,
                                   overwrite = F){
  logger::log_info("\t Starting the ambiguous junction pairing removal.")
  ## If rw_disk argument is set to TRUE and overwrite is set to FALSE, it looks
  ## for the results of the function in disk. If they have already been
  ## generated, the function is not executed and the files are returned.
  if (rw_disk & !overwrite & file.exists(paste0(main_samples_path, "all_distances_pruned.rds"))) {
    logger::log_info("\t\t Ignoring ambiguous pairing removal.") 
    logger::log_info("\t\t File ", main_samples_path, "all_distances_pruned.rds already exists.")
    logger::log_info("\t\t Loading file from disk.")
    return(readRDS(paste0(main_samples_path, "all_distances_pruned.rds")))
  }
  
  ## If the same novel junction overlaps two or more reference junctions, we
  ## keep one with the highest reads or the minimum distance. This is done per
  ## sample.
  logger::log_info("\t\t Per sample, filtering pairings that share a novel junction.")
  all_distances_pruned <- all_distances_raw %>%
    dplyr::group_by(sample, novel_junID) %>%
    dplyr::filter(ref_counts == max(ref_counts)) %>%
    dplyr::filter(abs(distance) == min(abs(distance))) %>%
    dplyr::ungroup()
  
  ## To remove ambiguous junctions, we remove the junctions in which the
  ## standard deviation of the distance is different than 0
  logger::log_info("\t\t For all samples, removing pairings that share a novel junction by .")
  all_distances_tidy <- all_distances_pruned %>%
    dplyr::group_by(novel_junID) %>%
    dplyr::mutate(
      distances_sd = distance %>% sd(),
      distances_sd = ifelse(is.na(distances_sd), 0, distances_sd)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(distances_sd == 0) %>%
    dplyr::select(-c(distances_sd))
  
  ## Save the output
  ## Recommended name: all_distances_pruned.rds
  if (rw_disk) {
    logger::log_info("\t\t Saving output to ", main_samples_path, "all_distances_pruned.rds")
    all_distances_tidy %>% saveRDS(paste0(main_samples_path, "all_distances_pruned.rds"))
  }
  
  logger::log_info("\t\t Successfully removed ambiguous junctions between all given projects")
  logger::log_info("")
  return(all_distances_tidy)
}

addCounts <- function(samples,
                      cluster_distances_pruned,
                      cluster_split_reads) {
  logger::log_info("\t\t\t Adding reads information to novel and reference junctions.")

  ## Generate the reads dataframes for novel junctions and reference junctions
  cluster_novel_counts <- cluster_split_reads %>%
    dplyr::filter(junID %in% unique(cluster_distances_pruned$novel_junID)) %>%
    dplyr::mutate(
      novel_n_individuals = rowSums(!is.na(across(all_of(samples)))),
      novel_sum_counts = rowSums(across(all_of(samples)), na.rm = T),
      novel_mean_counts = rowMeans(across(all_of(samples)), na.rm = T)
    ) %>%
    dplyr::select(-all_of(samples))

  cluster_ref_counts <- cluster_split_reads %>%
    dplyr::filter(junID %in% unique(cluster_distances_pruned$ref_junID)) %>%
    dplyr::mutate(
      ref_n_individuals = rowSums(!is.na(across(all_of(samples)))),
      ref_sum_counts = rowSums(across(all_of(samples)), na.rm = T),
      ref_mean_counts = rowMeans(across(all_of(samples)), na.rm = T)
    ) %>%
    dplyr::select(-all_of(samples))

  ## Merge both results in the previous distance dataframe
  cluster_distances_tidy <- cluster_distances_pruned %>%
    dplyr::left_join(cluster_novel_counts, by = c("novel_junID" = "junID")) %>%
    dplyr::left_join(cluster_ref_counts, by = c("ref_junID" = "junID"))

  return(cluster_distances_tidy)
}