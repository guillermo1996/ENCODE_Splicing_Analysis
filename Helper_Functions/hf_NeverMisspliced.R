#' Adds the never misspliced junctions to distances dataframe
#'
#' @param samples Vector containing the samples of the particular cluster.
#' @param cluster_split_reads Dataframe containing all the reads for each
#'   junction in the cluster.
#' @param cluster_distances_raw Dataframe with the cluster raw novel junctions
#'   and the associated reference junction with information about the distance
#'   and other reads statistics.
#' @param cluster_annotated_SR_details Dataframe containing the cluster
#'   junctions with their relevant information annotated.
#' @param cluster_distances_tidy Dataframe with the cluster pruned novel
#'   junctions and the associated reference junction.
#' @param cluster_path Path to where to save the files related to the cluster.
#' @param cluster_name Name of the cluster.
#' @param rw_disk Whether to store the results in disk. By default, TRUE.
#' @param overwrite Whether to overwrite previously generated results from the
#'   function. If set to FALSE and 'rw_disk' is set to TRUE, the function looks
#'   for the files in memory and loads them if possible. By default, FALSE.
#'
#' @return Dataframe with all the novel junctions and the associated reference junction and all the never mis-spliced junctions for the given cluster of samples. It also includes read statistics.
#' @export
addNeverMissplicedJunction <- function(samples,
                                       cluster_split_reads,
                                       cluster_distances_raw,
                                       cluster_annotated_SR_details,
                                       cluster_distances_tidy,
                                       cluster_path,
                                       cluster_name,
                                       rw_disk = T,
                                       overwrite = F){
  logger::log_info("\t\t\t Starting the generation of never mis-spliced junctions.")
  
  ## If rw_disk argument is set to TRUE and overwrite is set to FALSE, it looks
  ## for the results of the function in disk. If they have already been
  ## generated, the function is not executed and the files are returned.
  if (rw_disk & !overwrite & file.exists(paste0(cluster_path, cluster_name, "_distances_tidy_all.rds"))) {
    logger::log_info('\t\t\t\t Ignoring junction extraction. File ', cluster_path, cluster_name, '_distances_tidy_all.rds already exists.')
    return(readRDS(paste0(cluster_path, cluster_name, "_distances_tidy_all.rds")))
  }
  
  ############ Adding never mis-spliced information ############
  ##
  ##
  ## 1. Obtain a list of never mis-spliced junctions
  ## 
  ## 2. Generate read statistics of never mis-spliced junctions.
  ##
  ## 3. Appends the relevant information of the never mis-spliced junctions to
  ## the distances dataframe, previously containing all the relevant information
  ## from the novel junctions and their reference junctions.
  
  cluster_never_misspliced <- getNeverMissplicedJunctions(cluster_distances_raw,
                                                          cluster_annotated_SR_details)
  
  cluster_split_reads_never <- getNeverMissplicedReads(samples,
                                                       cluster_split_reads,
                                                       cluster_never_misspliced)
  
  cluster_distances_tidy_all <- appendNeverMissplicedIntrons(cluster_split_reads_never,
                                                             cluster_annotated_SR_details,
                                                             cluster_distances_tidy)
  
  ## Save the output
  ## Recommended name: cluster_distances_tidy_all.rds
  if (rw_disk) {
    logger::log_info("\t\t Saving output to ", cluster_path, cluster_name, "_distances_tidy_all.rds")
    cluster_distances_tidy_all %>% saveRDS(paste0(cluster_path, cluster_name, "_distances_tidy_all.rds"))
  }
  
  logger::log_info("\t\t\t Mis-spliced junctions added successfully.")
  return(cluster_distances_tidy_all)
}

#' Generates a list of never mis-spliced junctions
#'
#' Given the raw distances dataframe and the annotated junctions details, it
#' identifies the annotated junctions that do not overlap with any novel
#' junction under any circumstance (even if was removed as a potential reference
#' junction).
#'
#' @param cluster_distances_raw Dataframe with the cluster raw novel junctions
#'   and the associated reference junction with information about the distance
#'   and other reads statistics.
#' @param cluster_annotated_SR_details Dataframe containing the cluster
#'   junctions with their relevant information annotated.
#'
#' @return List of never mis-spliced junctions IDs.
#' @export
getNeverMissplicedJunctions <- function(cluster_distances_raw,
                                        cluster_annotated_SR_details) {
  logger::log_info("\t\t\t\t Extracting never mis-spliced junctions.")
  ## List of mis-spliced junctions
  miss_spliced_junctions <- cluster_distances_raw %>%
    pull(ref_junID) %>%
    unique()
  
  ## List of all annotated junctions
  all_junctions <- cluster_annotated_SR_details %>%
    filter(type == "annotated" & strand != "*") %>%
    pull(junID) %>%
    unique()
  
  ## Executing the set difference between all the annotated junctions and the
  ## mis-spliced junctions
  never_miss_spliced_junctions <- setdiff(all_junctions, miss_spliced_junctions)
  
  return(never_miss_spliced_junctions)
}

#' Extracts the reads for never mis-spliced junctions
#'
#' Given the list of never mis-spliced junctions and the read counts of every
#' junction, three different statistics about the reads per junction are
#' calculated: number of samples in which the junction is found, average number
#' of reads per sample and the total number of reads across all samples.
#'
#' @param samples Vector containing the samples of the particular cluster.
#' @param cluster_split_reads Dataframe containing all the reads for each junction in the cluster.
#' @param cluster_never_misspliced List of never mis-spliced junctions IDs.
#'
#' @return Read statistics for all the never mis-spliced junctions.
#' @export
getNeverMissplicedReads <- function(samples,
                                    cluster_split_reads,
                                    cluster_never_misspliced) {
  logger::log_info("\t\t\t\t Obtaining never mis-spliced reads information.")
  
  ## Generate the reads dataframes for never mis-spliced junctions
  cluster_split_reads_never <- cluster_split_reads %>%
    dplyr::filter(junID %in% cluster_never_misspliced) %>%
    dplyr::mutate(
      ref_n_individuals = rowSums(!is.na(across(all_of(samples)))),
      ref_sum_counts = rowSums(across(all_of(samples)), na.rm = T),
      ref_mean_counts = rowMeans(across(all_of(samples)), na.rm = T)
    ) %>%
    dplyr::select(-all_of(samples))
  
  return(cluster_split_reads_never)
}

#' Merge the distance dataframe with the never mis-spliced junctions
#'
#' From the distance dataframe containing all the information about the novel
#' junctions and their associated reference junction, we append the information
#' about the never mis-spliced junctions.
#'
#' @param cluster_split_reads_never Dataframe containing all the reads for each
#'   never misspliced junction in the cluster.
#' @param cluster_annotated_SR_details Dataframe containing the cluster
#'   junctions with their relevant information annotated.
#' @param cluster_distances_tidy Dataframe with the cluster pruned novel
#'   junctions and the associated reference junction.
#'
#' @return Dataframe containing the information of all novel junctions (and
#'   their associated reference junctions) and all the never mis-spliced
#'   junctions for the given cluster of samples. It also includes read
#'   statistics.
#' @export
appendNeverMissplicedIntrons <- function(cluster_split_reads_never,
                                         cluster_annotated_SR_details,
                                         cluster_distances_tidy) {
  logger::log_info("\t\t\t\t Adding never mis-spliced to the distances data.")
  ## Merge the never mis-spliced reads dataframe with the annotated junction
  ## dataframe to extract all their relevant information.
  df_never <- cluster_split_reads_never %>%
    dplyr::left_join(cluster_annotated_SR_details %>% select(junID, seqnames, start, end, strand, width),
                     by = "junID") %>%
    dplyr::rename(ref_junID = junID,
                  ref_seq = seqnames,
                  ref_start = start,
                  ref_end = end,
                  ref_strand = strand,
                  ref_width = width)
  
  
  ## Merge the distances dataframe with the newly created never mis-spliced
  ## junctions information. Also, append "gene_id_junction",
  ## "gene_name_junction" and "tx_id_junction" information of all reference and
  ## never mis-spliced junctions.
  df_all <- dplyr::bind_rows(cluster_distances_tidy, df_never)
  
  df_all <- df_all %>%
    dplyr::left_join(cluster_annotated_SR_details %>% select(junID, gene_id_junction, gene_name_junction, tx_id_junction),
                     by = c("ref_junID" = "junID")) %>%
    dplyr::rename(gene_id = gene_id_junction,
                  gene_name = gene_name_junction)
  return(df_all)
}