#' Creates the subfolders for a particular RBP
#'
#' @param RBP_metadata Dataframe containing all the metadata for the RBPs/NMDs.
#' @param RBP_path Path to the folder where the generated files will be stored.
#' @param RBP_clusters The different clusters found for the target gene.
#' @param generate_script Whether to generate a download script for the samples
#'   by cluster. By default, TRUE.
#'
#' @return NULL
#' @export
createSubFolders <- function(RBP_metadata,
                             RBP_path, 
                             RBP_clusters,
                             generate_script = T) {
  logger::log_info("\t\t Generating the folder structure in ", RBP_path)
  ## Loop through the clusters and generate a specific folder for each one.
  for (cluster in RBP_clusters) {
    cluster_path <- paste0(RBP_path, cluster, "/")
    cluster_metadata <- RBP_metadata %>% filter(experiment_type == cluster)
    
    dir.create(cluster_path, recursive = T, showWarnings = F)
    
    ## If generate_script is set to TRUE, a download script will be generated in
    ## the cluster folder. However, it is not recommended to use this script to
    ## download the files.
    if (generate_script == T) {
      generateDownloadScript(cluster_metadata, cluster_path)
    }
  }
}

#' Generates the download script
#'
#' Given the metadata for the current target gene and cluster (either case or
#' control), it generates a script to download the BAM files in the
#' corresponding location. The script can be executed to download the BAM files,
#' but it is not recommended.
#'
#' @param cluster_metadata Dataframe containing all the metadata for the current
#'   cluster and target gene.
#' @param cluster_path Path to the cluster folder where the generated files will
#'   be stored.
#'
#' @return NULL
#' @export
generateDownloadScript <- function(cluster_metadata,
                                   cluster_path) {
  logger::log_info("\t\t Generating the download scripts in ", cluster_path)
  
  ## Generate the script's text
  download_script <- "# Script to download the .bam files\n"
  for (i in 1:nrow(cluster_metadata)) {
    download_link <- getDownloadLinkMetadata(cluster_metadata[i, ])
    download_script <- paste0(download_script, "wget -c ", download_link, "\n")
  }
  
  ## Writing of the script.
  tryCatch(
    {
      file_path <- paste0(cluster_path, "/download.sh")
      file_conn <- file(file_path)
      writeLines(download_script, file_conn)
      close(file_conn)
      
      system2(command = "chmod", args = c("+x", file_path))
    },
    error = function(e) print(e)
  )
}

#' Generates the download link for a given sample
#'
#' @param cluster_metadata_row Row containing all the metadata for the current
#'   cluster and target gene.
#'
#' @return Download link for the sample.
#' @export
getDownloadLinkMetadata <- function(cluster_metadata_row){
  sample_id <- cluster_metadata_row %>% pull(sample_id)
  sample_format <- cluster_metadata_row %>% pull(file_format)
  
  return(paste0("https://www.encodeproject.org/files/", sample_id, "/@@download/", sample_id, ".", sample_format))
}

#' Download the BAM files and extracts them into JUNC file
#'
#' Given the metadata for the current cluster (either case or control), this
#' function reads the download link and automatically start the process of
#' downloading and extracting the BAM files. The process requires the path for
#' \href{http://www.htslib.org/}{samtools} and
#' \href{https://regtools.readthedocs.io/en/latest/}{regtools}, but they can be
#' left empty if they are already in the default PATH.
#'
#' It is possible to control certain parameters from samtools, like the number
#' of threads or the maximum memory per core. More information in the official
#' samtools \href{http://www.htslib.org/doc/samtools-sort.html}{documentation}.
#'
#' Tested on samtools 1.16.1 and regtools 0.5.2
#'
#' @param RBP_metadata Dataframe containing all the metadata for the RBPs/NMDs.
#'   It is required to be have the field "file_format", "experiment_type" and
#'   "sample_id" in every row.
#' @param RBP_path Path to the folder where the generated files will be stored.
#' @param num_cores Number of multiprocessing cores to use. Memory requirements
#'   significantly increase with the number of cores.
#' @param samtools_threads Number of threads to use in the samtools sort
#'   command.
#' @param samtools_memory Maximum memory per core to use in the samtools sort
#'   command.
#' @param samtools_path Path to the samtools executable. Can be left empty if
#'   samtools is in default PATH.
#' @param regtools_path Path to the regtools executable. Can be left empty if
#'   regtools is in default PATH.
#' @param overwrite Whether to overwrite previously generated results from the
#'   function. If set to FALSE and 'rw_disk' is set to TRUE, the function looks
#'   for the files in memory and loads them if possible. By default, FALSE.
#'
#' @return NULL
#' @export
downloadExtractBamFiles <- function(RBP_metadata,
                                    RBP_path,
                                    num_cores = 4,
                                    samtools_threads = 1,
                                    samtools_memory = "5G",
                                    samtools_path = "",
                                    regtools_path = "",
                                    overwrite = F) {
  logger::log_info("\t Starting the download and extraction process.")
  if(checkDownloadedFiles(RBP_metadata, RBP_path)){
    logger::log_info("\t\t Ignoring download and extraction. All junctions already extracted!")
    return()
  }
  
  ## Multiprocessing generation. Add argument "output.file" to
  ## parallel::makeCluster() if you want to output information about the parallel
  ## execution
  cl <- parallel::makeCluster(num_cores, outfile = "")
  doParallel::registerDoParallel(cl)
  
  ## Multiprocessing loop
  metrics <- foreach(i = 1:nrow(RBP_metadata), .export = "getDownloadLinkMetadata", .packages = "tidyverse") %dopar% {
    ## Definition of the variables
    sample_id <- RBP_metadata[i, ] %>% pull(sample_id)
    sample_cluster <- RBP_metadata[i, ] %>% pull(experiment_type)
    sample_target_gene = RBP_metadata[i, ] %>% pull(target_gene)
    download_link <- getDownloadLinkMetadata(RBP_metadata[i, ])
    file_path <- paste0(RBP_path, sample_cluster, "/", sample_id, ".bam")
    sort_path <- paste0(RBP_path, sample_cluster, "/", sample_id, ".bam.sort")
    junc_path <- paste0(RBP_path, sample_cluster, "/", sample_id, ".bam.sort.s0.junc")
    
    ## Check if file BAM file is already extracted
    if(!overwrite & file.exists(junc_path)) {
      return(paste0("\t\t Ignoring extraction and download. Junction file for sample ", sample_id, " (", sample_target_gene, " - ", sample_cluster, ") already found!"))
    }
    
    ## Download the BAM file
    #logger::log_info("Starting download of sample ", sample_id, " (", sample_target_gene, " - ", sample_cluster, ")")
    download.file(download_link, file_path, method = "wget", extra = "-c --quiet --no-check-certificate")
    
    ## Extract the BAM file into JUNC using samtools and regtools
    #logger::log_info("Starting the sorting process.")
    system2(command = paste0(samtools_path, "samtools"), args = c(
      "sort", file_path,
      "-o", sort_path,
      "--threads", samtools_threads,
      "-m", samtools_memory
    ))
    if(!file.exists(sort_path)) return("Sorting failed")
    
    #logger::log_info("\t Starting the indexing process.")
    system2(command = paste0(samtools_path, "samtools"), args = c(
      "index", sort_path,
      "-@", samtools_threads
    ))
    
    #logger::log_info("\t Starting the extraction process.")
    system2(command = paste0(regtools_path, "regtools"), args = c(
      "junctions extract", sort_path,
      "-m 25",
      "-M 1000000",
      "-s 0",
      "-o", junc_path
    ))
    
    ## Remove the files that are not necessary
    system2(command = "rm", args = c(file_path))
    system2(command = "rm", args = c(sort_path))
    system2(command = "rm", args = c(paste0(sort_path, ".bai")))
    if(!file.exists(junc_path)){
      return(c(
        paste0("\t\t Error extracting sample ", sample_id, " (", sample_target_gene, " - ", sample_cluster, ")", "."),
        paste0("\t\t Removing all intermediary files. Please run the analysis again.")
      ))
    }else{
      return(paste0("\t\t Successfully extracted the junctions. Removing all intermediary files (BAM included)."))
    }
  }
  ## Stop the parallel cluster
  parallel::stopCluster(cl)
  
  ## Print the metrics
  invisible(lapply(metrics, function(x) logger::log_info(x[1])))
  invisible(lapply(metrics, function(x) if (!is.na(x[2])) logger::log_info(x[2])))
}

#' Checks if all the JUNC files in the input dataframe already exists
#'
#' @param RBP_metadata Dataframe containing all the metadata for the RBPs/NMDs.
#' @param RBP_path Path to the folder where the generated files should be
#'   stored.
#'
#' @return Whether the JUNC files for the input metadata exists or not.
#' @export
checkDownloadedFiles <- function(RBP_metadata,
                                 RBP_path){
  file_names <- apply(RBP_metadata, 1, function(x) {
    paste0(RBP_path, x["experiment_type"], "/", x["sample_id"], ".bam.sort.s0.junc")
  })
  
  return(all(file.exists(file_names)))
}
