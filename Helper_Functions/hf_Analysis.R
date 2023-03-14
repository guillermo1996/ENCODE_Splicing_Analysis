#' Get all introns from the databases
#'
#' @param target_RBPs Vector of the target genes to study.
#' @param RBP_path Path to where the RBP files are located.
#' @param required_clusters  Vector of required clusters to study. Defaults to
#'   'c("case", "control")'
#' @param prune_columns Whether to prune the columns in the databases.
#'   Recommended to reduce memory usage. Default to TRUE.
#' @param num_cores Number of cores to parallelize the process.
#' @param file_output Path to where to store the results. Leave empty for no
#'   storing.
#' @param overwrite Whether to overwrite the previously generated file. If the
#'   'file_output' path exists, it will load the variable instead of generating
#'   a new one. Defaults to FALSE.
#'
#' @return Dataframe containing all introns obtained from the analysis.
#' @export
getAllIntrons <- function(target_RBPs,
                          RBP_path,
                          required_clusters = c("case", "control"),
                          prune_columns = T,
                          num_cores = 4,
                          file_output = "",
                          overwrite = F){
  if(!overwrite & file.exists(file_output)){
    global_introns <- readRDS(file_output)
    return(global_introns)
  }
  
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  
  global_introns <- foreach(i = seq(target_RBPs), .packages = c("tidyverse", "foreach")) %dopar%{
    target_gene <- target_RBPs[i]
    
    cluster_introns <- foreach(j = seq(required_clusters)) %do%{
      cluster <- required_clusters[j]
      
      file_path <- paste0(RBP_path, target_gene, "/", cluster, "/", cluster, "_db_introns.rds")
      if(!file.exists(file_path)) return(tibble::tibble())
      
      db_introns <- readRDS(file_path) %>%
        dplyr::mutate(project = target_gene, cluster = cluster)
      
      if(prune_columns){
        db_introns <- db_introns %>%
          dplyr::select(ref_junID, ref_type, MSR_Donor, MSR_Acceptor, ref_reads, project, cluster)
      }
      
      return(db_introns)
    }
    
    return(cluster_introns)
  } %>% dplyr::bind_rows()
  parallel::stopCluster(cl)
  
  if(file_output != ""){
    global_introns %>% saveRDS(file_output)
  }
  
  return(global_introns)
}


#' Get the common introns across samples
#'
#' @param target_RBPs Vector of the target genes to study.
#' @param RBP_path Path to where the RBP files are located.
#' @param required_clusters  Vector of required clusters to study. Defaults to
#'   'c("case", "control")'
#' @param prune_columns Whether to prune the columns in the databases.
#'   Recommended to reduce memory usage. Default to TRUE.
#' @param num_cores Number of cores to parallelize the process.
#' @param file_output Path to where to store the results. Leave empty for no
#'   storing.
#' @param overwrite Whether to overwrite the previously generated file. If the
#'   'file_output' path exists, it will load the variable instead of generating
#'   a new one. Defaults to FALSE.
#'
#' @return Dataframe containing only the common introns across all samples.
#' @export
getCommonIntrons <- function(target_RBPs,
                             RBP_path,
                             required_clusters = c("case", "control"),
                             prune_columns = T,
                             num_cores = 4,
                             file_output = "",
                             overwrite = F){
  if(!overwrite & file.exists(file_output)){
    common_introns <- readRDS(file_output)
    return(common_introns)
  }
  
  global_introns <- getAllIntrons(target_RBPs = target_RBPs,
                                  RBP_path = RBP_path,
                                  required_clusters = required_clusters,
                                  prune_columns = T,
                                  num_cores = 8,
                                  file_output = global_introns_path,
                                  overwrite = F)
  
  numer_of_projects <- global_introns$project %>% unique() %>% length()
  number_of_clusters <- global_introns$cluster %>% unique() %>% length()
  common_introns <- global_introns %>%
    dplyr::group_by(ref_junID) %>%
    dplyr::filter(n() == (numer_of_projects*number_of_clusters)) %>%
    dplyr::ungroup()
  
  if(file_output != ""){
    common_introns %>% saveRDS(file_output)
  }
  
  return(common_introns)
}

#' #' Get all novel junctions from the databases
#'
#' @param target_RBPs Vector of the target genes to study.
#' @param common_introns Dataframe containing only the common introns across all
#'   samples.
#' @param RBP_path Path to where the RBP files are located.
#' @param required_clusters  Vector of required clusters to study. Defaults to
#'   'c("case", "control")'
#' @param prune_columns Whether to prune the columns in the databases.
#'   Recommended to reduce memory usage. Default to TRUE.
#' @param num_cores Number of cores to parallelize the process.
#' @param file_output Path to where to store the results. Leave empty for no
#'   storing.
#' @param overwrite Whether to overwrite the previously generated file. If the
#'   'file_output' path exists, it will load the variable instead of generating
#'   a new one. Defaults to FALSE.
#'
#' @return Dataframe containing all the novel junctions associated to the common
#'   introns.
#' @export
getCommonNovel <- function(target_RBPs,
                           common_introns,
                           RBP_path,
                           required_clusters = c("case", "control"),
                           prune_columns = T,
                           num_cores = 4,
                           file_output = "",
                           overwrite = F){
  if(!overwrite & file.exists(file_output)){
    global_novel <- readRDS(file_output)
    return(global_novel)
  }
  
  common_introns_list <- common_introns$ref_junID %>% unique()
  
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  
  global_novel <- foreach(i = seq(target_RBPs), .packages = "tidyverse") %do%{
    target_gene <- target_RBPs[i]
    
    cluster_novel <- foreach(j = seq(required_clusters)) %do%{
      cluster <- required_clusters[j]
      
      file_path <- paste0(RBP_path, target_gene, "/", cluster, "/", cluster, "_db_novel.rds")
      if(!file.exists(file_path)) return(tibble::tibble())
      
      db_novel <- readRDS(file_path) %>%
        dplyr::mutate(project = target_gene, cluster = cluster)
      
      if(prune_columns){
        db_novel <- db_novel %>%
          dplyr::select(novel_junID, ref_junID, novel_type, novel_reads, distance, project, cluster)
      }
      
      db_novel <- db_novel %>%
        dplyr::filter(ref_junID %in% common_introns_list)
      
      return(db_novel)
    }
    
    return(cluster_novel)
  } %>% dplyr::bind_rows()
  parallel::stopCluster(cl)
  
  if(file_output != ""){
    global_novel %>% saveRDS(file_output)
  }
  
  return(global_novel)
}

generateMSRtests <- function(target_RBPs,
                             cluster_case,
                             cluster_control,
                             MSR_Table,
                             num_cores = 4,
                             file_output = "",
                             overwrite = F){
  if(!overwrite & file.exists(file_output)){
    MSR_tests <- readRDS(file_output)
    return(MSR_tests)
  }
  
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  MSR_tests <- foreach(i = seq(length(target_RBPs)), .packages = c("tidyverse")) %dopar%{
    project <- target_RBPs[i]
    
    MSR_project <- MSR_Table %>%
      dplyr::select(ref_junID,
                    case = paste0(project, "_", cluster_case),
                    control = paste0(project, "_", cluster_control)) %>%
      tidyr::drop_na() %>%
      tidyr::pivot_longer(cols = c("case", "control"), names_to = "group", values_to = "MSR") %>%
      dplyr::arrange(ref_junID)
    
    #MSR_case <- MSR_project %>% dplyr::filter(group == "case") %>% dplyr::pull(MSR)
    #MSR_control <- MSR_project %>% dplyr::filter(group == "control") %>% dplyr::pull(MSR)
    
    wilcox_test <- rstatix::wilcox_test(MSR_project, MSR ~ group, paired = TRUE, alternative = "greater")
    wilcox_effsize <- rstatix::wilcox_effsize(MSR_project, MSR ~ group, paired = TRUE, alternative = "greater")
    
    tibble(target_gene = project,
           statistical_test = NA,
           H0 = NA,
           H1 = NA,
           p.value = wilcox_test$p,
           effect_size = wilcox_effsize$effsize,
           magnitude = wilcox_effsize$magnitude)
  } %>% dplyr::bind_rows()
  parallel::stopCluster(cl)
  
  if(file_output != ""){
    MSR_tests %>% saveRDS(file_output)
  }
  
  return(MSR_tests)
}

#' Add the target gene functional category to the MSR test dataframe
#'
#' @param MSR_tests
#' @param metadata_RBPs
#'
#' @return
#' @export
#'
#' @examples
addMSRcategories <- function(MSR_tests, metadata_RBPs){
  MSR_tests <- MSR_tests %>%
    dplyr::left_join(metadata_RBPs %>% 
                       dplyr::select(target_gene, Category) %>% 
                       dplyr::distinct(),
                     by = "target_gene") %>%
    dplyr::relocate(Category, .after = target_gene)
  
  return(MSR_tests)
}

#' Applies the bonferroni correction to a MSR test dataframe
#'
#' @param MSR_tests
#'
#' @return
#' @export
#'
#' @examples
addBonferroniCorrection <- function(MSR_tests){
  MSR_tests <- MSR_tests %>% 
    dplyr::mutate(p.value.bonferroni = NA, .after = p.value)
  
  MSR_tests$p.value.bonferroni <- p.adjust(MSR_tests$p.value, method = "bonferroni")
  return(MSR_tests)
}