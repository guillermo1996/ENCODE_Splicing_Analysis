#' Pipeline to generate the novel and reference junctions database dataframes
#'
#' @param cluster_distances_tidy_all Dataframe with all the novel junctions and
#'   the associated reference junction and all the never mis-spliced junctions
#'   for the given cluster of samples. It also includes read statistics.
#' @param cluster_annotated_SR_details Dataframe containing the cluster
#'   junctions with their relevant information annotated.
#' @param cluster_path Path to where to save the files related to the cluster.
#' @param cluster_name Name of the cluster.
#' @param u12_introns_path Path to where the u12 (or minor) intron list is
#'   located in disk.
#' @param u2_introns_path Path to where the u2 (or major) intron list is located
#'   in disk.
#' @param rw_disk Whether to store the results in disk. By default, TRUE.
#' @param overwrite Whether to overwrite previously generated results from the
#'   function. If set to FALSE and 'rw_disk' is set to TRUE, the function looks
#'   for the files in memory and loads them if possible. By default, FALSE.
#'
#' @return NULL
#' @export
generateDB <- function(cluster_distances_tidy_all,
                       cluster_annotated_SR_details,
                       cluster_path,
                       cluster_name,
                       u12_introns_path = "/home/grocamora/RytenLab-Research/Additional_files/minor_introns_tidy.rds",
                       u2_introns_path = "/home/grocamora/RytenLab-Research/Additional_files/major_introns_tidy.rds",
                       rw_disk = T,
                       overwrite = F) {
  logger::log_info("\t\t\t Generation the DB files.")
  
  ## If rw_disk argument is set to TRUE and overwrite is set to FALSE, it looks
  ## for the results of the function in disk. If they have already been
  ## generated, the function is not executed and the files are returned.
  if (rw_disk & !overwrite & 
      file.exists(paste0(cluster_path, cluster_name, "_db_introns.rds")) &
      file.exists(paste0(cluster_path, cluster_name, "_db_novel.rds"))){
    logger::log_info("\t\t\t\t Ignoring the DB generation process.")
    logger::log_info("\t\t\t\t File ", cluster_path, cluster_name, "_db_introns.rds already exists.")
    logger::log_info("\t\t\t\t File ", cluster_path, cluster_name, "_db_novel.rds already exists.")
    return()
  }
  
  df_misspliced <- cluster_distances_tidy_all %>%
    dplyr::filter(!is.na(novel_junID))
  
  df_never <- cluster_distances_tidy_all %>%
    dplyr::filter(is.na(novel_junID))
  
  db_novel <- generateNovelDB(df_misspliced, cluster_annotated_SR_details)
  db_introns <- generateIntronDB(df_misspliced, df_never, cluster_annotated_SR_details)
  
  ## Add biotype percentage
  db_introns <- db_introns %>%
    dplyr::left_join(cluster_annotated_SR_details %>% select(junID, protein_coding, lncRNA),
                     by = c("ref_junID" = "junID"))
  
  ## Add intron type
  db_introns <- addIntronType(db_introns,
                              u12_introns_path = u12_introns_path,
                              u2_introns_path = u2_introns_path)
  
  if(rw_disk){
    logger::log_info("\t\t\t\t Saving output to ", cluster_path, cluster_name, "_db_novel.rds")
    db_novel %>% saveRDS(paste0(cluster_path, cluster_name, "_db_novel.rds"))
    logger::log_info("\t\t\t\t Saving output to ", cluster_path, cluster_name, "_db_introns.rds")
    db_introns %>% saveRDS(paste0(cluster_path, cluster_name, "_db_introns.rds"))
  }
  
  logger::log_info("\t\t\t\t DB generation process executed successfully.")
}

#' Generates the novel junctions DB
#'
#' @param df_misspliced Dataframe containing the information of all novel
#'   junctions (and their associated reference junctions) and their read
#'   statistics.
#' @param cluster_annotated_SR_details Dataframe containing the cluster
#'   junctions with their relevant information annotated.
#'
#' @return Dataframe with all the novel junctions and their most relevant
#'   information.
#' @export
generateNovelDB <- function(df_misspliced,
                            cluster_annotated_SR_details) {
  logger::log_info("\t\t\t\t Generation the db_novel file.")
  ## Select and rename the relevant columns
  db_novel <- df_misspliced %>%
    dplyr::mutate(novel_coordinates = paste0("chr", novel_seq, ":", novel_start, "-", novel_end, ":", novel_strand)) %>%
    dplyr::select(novel_junID,
                  ref_junID,
                  novel_coordinates,
                  novel_seq,
                  novel_start,
                  novel_end,
                  novel_strand,
                  novel_width,
                  novel_type = type,
                  novel_reads = novel_sum_counts,
                  distance)
  
  db_novel <- db_novel %>%
    dplyr::left_join(cluster_annotated_SR_details %>% select(junID, ss5score, ss3score),
                     by = c("novel_junID" = "junID")) %>%
    dplyr::rename(novel_ss5score = ss5score,
                  novel_ss3score = ss3score) %>%
    dplyr::left_join(cluster_annotated_SR_details %>% select(junID, ss5score, ss3score),
                     by = c("ref_junID" = "junID")) %>%
    dplyr::rename(ref_ss5score = ss5score,
                  ref_ss3score = ss3score)
  
  return(db_novel)
}

#' Generates the reference junctions DB
#'
#' @param df_misspliced Dataframe containing the information of all novel
#'   junctions (and their associated reference junctions) and their read
#'   statistics.
#' @param df_never Dataframe containing the information of all never mis-spliced
#'   junctions and their read statistics.
#'
#' @return Dataframe with all the reference junctions and never mis-spliced
#'   junctions with their most relevant information.
#' @export
generateIntronDB <- function(df_misspliced, 
                             df_never,
                             cluster_annotated_SR_details) {
  logger::log_info("\t\t\t\t Generation the db_intron file.")
  ## Calculate the mis-splicing ratio of every reference junction
  logger::log_info("\t\t\t\t\t Calculating the mis-splicing ratios.")
  df_ref <- df_misspliced %>%
    dplyr::group_by(ref_junID, type) %>%
    dplyr::mutate(MSR = sum(novel_sum_counts) / (sum(novel_sum_counts) + ref_sum_counts)) %>%
    dplyr::distinct(MSR, .keep_all = T) %>%
    dplyr::ungroup()
  
  ## Combine the information of novel donors and novel acceptors when
  ## associated to the same reference junction
  db_ref <- df_ref %>%
    tidyr::spread(key = type, value = MSR, fill = NA) %>%
    dplyr::group_by(ref_junID) %>%
    dplyr::mutate(
      MSR_Donor = mean(novel_donor, na.rm = T),
      MSR_Acceptor = mean(novel_acceptor, na.rm = T)
    ) %>%
    dplyr::distinct(ref_junID, .keep_all = T) %>%
    dplyr::ungroup()
  
  ## Remove all the novel junctions information (only keep the annotated
  ## junctions information)
  logger::log_info("\t\t\t\t\t Pruning the resulted data.")
  db_ref <- db_ref %>%  
    dplyr::select(-all_of(starts_with("novel")), -distance)
  
  ## Remove all the novel junctions information from the never mis-spliced
  ## junctions
  db_never <- df_never %>%
    dplyr::select(-all_of(starts_with("novel")), -type, -distance)
  
  ## Combine the dataframes for the reference junctions and the never
  ## mis-spliced junctions
  db_introns <- dplyr::bind_rows(db_ref, db_never)
  db_introns[is.na(db_introns$MSR_Donor), "MSR_Donor"] <- 0
  db_introns[is.na(db_introns$MSR_Acceptor), "MSR_Acceptor"] <- 0
  
  ## Better present the data
  db_introns <- db_introns %>%
    dplyr::mutate(ref_coordinates = paste0("chr", ref_seq, ":", ref_start, "-", ref_end, ":", ref_strand)) %>%
    dplyr::select(-ref_n_individuals, -ref_mean_counts) %>%
    dplyr::rename(ref_reads = ref_sum_counts) %>%
    dplyr::relocate(ref_coordinates, .after = ref_junID)
  
  ## Add MaxEntScore
  db_introns <- db_introns %>%
    dplyr::left_join(cluster_annotated_SR_details %>% select(junID, ss5score, ss3score),
                     by = c("ref_junID" = "junID")) %>%
    dplyr::rename(ref_ss5score = ss5score,
                  ref_ss3score = ss3score) %>%
    tidyr::unnest(c(gene_id, gene_name))
  
  ## Annotate the type of annotated junction as "both", "novel", "donor" or
  ## "never"
  logger::log_info("\t\t\t\t\t Annotate the introns in mis-splicing categories.")
  db_introns <- db_introns %>%
    dplyr::rowwise() %>%
    dplyr::mutate(ref_type = missplicingClass(MSR_Donor, MSR_Acceptor)) %>%
    dplyr::relocate(ref_type, .after = ref_coordinates)
  
  rm(df_ref, db_ref, db_never)
  invisible(gc())
  return(db_introns)
}

#' Categorize the reference junctions in "novel", "donor", "acceptor" or "never"
#'
#' @param novel_donor_ratio Mis-splicing ratio of the novel donors.
#' @param novel_acceptor_ratio Mis-splicing ratio of the novel acceptor.
#'
#' @return Reference junction category.
#' @export
missplicingClass <- function(novel_donor_ratio, novel_acceptor_ratio) {
  ref_junction_category = ""
  if (novel_donor_ratio > 0 & novel_acceptor_ratio > 0) {
    ref_junction_category = "both"
  } else if (novel_donor_ratio > 0 & novel_acceptor_ratio == 0) {
    ref_junction_category = "donor"
  } else if (novel_donor_ratio == 0 & novel_acceptor_ratio > 0) {
    ref_junction_category = "acceptor"
  } else {
    ref_junction_category = "never"
  }
  
  return(ref_junction_category)
}

#' Categorize the introns in u12 or u2
#'
#' @param db_introns Dataframe with all the reference junctions and never
#'   mis-spliced junctions with their most relevant information.
#' @param u12_introns_path Path to where the u12 (or minor) intron list is
#'   located in disk.
#' @param u2_introns_path Path to where the u2 (or major) intron list is located
#'   in disk.
#'
#' @return Dataframe with all the reference junctions and never mis-spliced
#'   junctions with their most relevant information.
#' @export
addIntronType <- function(db_introns,
                          u12_introns_path = "",
                          u2_introns_path = ""){
  logger::log_info("\t\t\t\t Adding the intron types (u2 or u12).")
  ## Load the variables
  u12_introns <- readRDS(u12_introns_path) %>% GenomicRanges::GRanges()
  u2_introns <- readRDS(u2_introns_path) %>% GenomicRanges::GRanges()
  
  GRdata <- db_introns %>%
    dplyr::select(ref_junID,
                  seqnames = ref_seq,
                  start = ref_start,
                  end = ref_end,
                  strand = ref_strand) %>%
    dplyr::mutate(u12_intron = F, u2_intron = F) %>%
    GenomicRanges::GRanges()
  
  
  ## Minor introns
  u12_overlaps <- GenomicRanges::findOverlaps(query = u12_introns,
                                              subject = GRdata,
                                              ignore.strand = F,
                                              type = "equal")
  ## Major introns
  u2_overlaps <- GenomicRanges::findOverlaps(query = u2_introns,
                                             subject = GRdata,
                                             ignore.strand = F,
                                             type = "equal")
  
  ## Adding the type
  GRdata <- GRdata %>% tibble::as_tibble()
  GRdata[S4Vectors::subjectHits(u12_overlaps), "u12_intron"] <- T
  GRdata[S4Vectors::subjectHits(u2_overlaps), "u2_intron"] <- T
  
  db_introns <- db_introns %>%
    left_join(GRdata %>% select(ref_junID, u12_intron, u2_intron),
              by = "ref_junID")
  
  return(db_introns)
}