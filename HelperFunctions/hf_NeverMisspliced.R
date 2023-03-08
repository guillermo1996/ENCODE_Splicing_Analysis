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

addNeverMissplicedJunction <- function(samples,
                                       cluster_split_reads,
                                       cluster_distances_raw,
                                       cluster_annotated_SR_details,
                                       cluster_distances_tidy,
                                       project_path,
                                       cluster_name,
                                       rw_disk = T,
                                       overwrite = F){
  logger::log_info("\t\t\t Starting the generation of never mis-spliced junctions.")
  
  ## If rw_disk argument is set to TRUE and overwrite is set to FALSE, it looks
  ## for the results of the function in disk. If they have already been
  ## generated, the function is not executed and the files are returned.
  if (rw_disk & !overwrite & file.exists(paste0(project_path, cluster_name, "_distances_tidy_all.rds"))) {
    logger::log_info('\t\t\t\t Ignoring junction extraction. File ', project_path, cluster_name, '_distances_tidy_all.rds already exists.')
    return(readRDS(paste0(project_path, cluster_name, "_distances_tidy_all.rds")))
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
    logger::log_info("\t\t Saving output to ", project_path, cluster_name, "_distances_tidy_all.rds")
    cluster_distances_tidy_all %>% saveRDS(paste0(project_path, cluster_name, "_distances_tidy_all.rds"))
  }
  
  logger::log_info("\t\t\t Mis-spliced junctions added successfully.")
  return(cluster_distances_tidy_all)
}

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