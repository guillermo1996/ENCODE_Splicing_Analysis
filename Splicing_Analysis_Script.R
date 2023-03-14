# 1. Load libraries and variables ----
## Required libraries ----
shhh <- suppressPackageStartupMessages
shhh(library(here, warn.conflicts = F))
shhh(library(logger, warn.conflicts = F))
shhh(library(GenomicRanges, warn.conflicts = F))
shhh(library(doParallel, warn.conflicts = F))
shhh(library(foreach, warn.conflicts = F))
shhh(library(tidyverse, warn.conflicts = F))
options(dplyr.summarise.inform = FALSE)
options(lifecycle_verbosity = "warning")

## Load helper functions ----
source(here::here("Helper_Functions/hf_DownloadExtract.R"))
source(here::here("Helper_Functions/hf_JunctionReadAnnotate.R"))
source(here::here("Helper_Functions/hf_Distances.R"))
source(here::here("Helper_Functions/hf_NeverMisspliced.R"))
source(here::here("Helper_Functions/hf_GenerateDB.R"))

## Logger options ----
log_file <- here::here("Splicing_Analysis.log")
logger::log_appender(logger::appender_tee(log_file, append = T))
logger_layout <- logger::layout_glue_generator(format = '[{time}] [{level}] {msg}')
logger::log_layout(logger_layout)

## Relevant Paths ----
main_samples_path <- here::here("RBPs/")
metadata_path <- here::here("Metadata/metadata_complete.tsv")

tools_path <- "/home/grocamora/tools/"
regtools_path <- paste0(tools_path, "regtools/build/")
samtools_path <- paste0(tools_path, "samtools/bin/")
bedtools_path <- paste0(tools_path, "bedtools/")
fordownload_path <- paste0(tools_path, "fordownload/")

additional_files_path <- "/home/grocamora/RytenLab-Research/Additional_files/"
fasta_path <- paste0(additional_files_path, "Homo_sapiens.GRCh38.dna.primary_assembly.fa")
gtf_path <- paste0(additional_files_path, "Homo_sapiens.GRCh38.105.chr.gtf")
blacklist_path <- paste0(additional_files_path, "hg38-blacklist.v2.bed")
u12_introns_path <- paste0(additional_files_path, "minor_introns_tidy.rds")
u2_introns_path <- paste0(additional_files_path, "major_introns_tidy.rds")

## Script parameters ----

### Whether to write files to disk and to overwrite previously created files. If
### rw_disk is set to TRUE, the files exists and overwrite is set to FALSE, the
### variables are loaded from disk.
rw_disk <- T      
overwrite <- F    

### Download options
skip_download = F
download_cores = 2

### Junction extraction options
samtools_threads = 2
samtools_memory = "5G"

### Multiprocessing
processing_cores = 8

## Load metadata ----
metadata <- readr::read_delim(metadata_path, show_col_types = F)

target_RBPs <- metadata %>%
  dplyr::filter(if_any(c(Splicing_regulation, Spliceosome, Exon_junction_complex, NMD), ~ . != 0)) %>%
  dplyr::pull(target_gene) %>%
  unique()
metadata_RBPs <- metadata %>% dplyr::filter(target_gene %in% target_RBPs)

# 2. Download the BAM files ----
if(skip_download){
  logger::log_info("Download and extraction of the BAM files is skipped. 
                   Please set 'skip_download' to FALSE if you want to execute the download and extraction.")
}else{
  logger::log_info("Starting the download BAM files process.")
  for(target_RBP in target_RBPs){
    logger::log_info("\t Downloading target RBP: ", target_RBP)
    
    ## Target RBP variables
    RBP_metadata <- metadata_RBPs %>% 
      dplyr::filter(target_gene == target_RBP) %>% 
      dplyr::arrange(experiment_type)
    RBP_path <- paste0(main_samples_path, target_RBP, "/")
    RBP_clusters <- RBP_metadata %>%
      dplyr::pull(experiment_type) %>% unique()
    
    ## If the junctions files are already found
    if(checkDownloadedFiles(RBP_metadata, RBP_path)){
      logger::log_info("\t\t Ignoring download and extraction. All junctions already extracted!")
      next
    }
    
    ## Create the subfolders
    createSubFolders(RBP_metadata,
                     RBP_path,
                     RBP_clusters,
                     generate_script = T)
    
    downloadExtractBamFiles(RBP_metadata = RBP_metadata,
                            RBP_path = RBP_path,
                            num_cores = download_cores,
                            samtools_threads = samtools_threads,
                            samtools_memory = samtools_memory,
                            samtools_path = samtools_path,
                            regtools_path = regtools_path,
                            overwrite = overwrite)
  }
}

# 3. Junction Reading and Annotation ----
logger::log_info("Starting splicing noise analysis for ENCODE samples.")
all_reads_combined <- junctionReading(metadata = metadata_RBPs, 
                                      main_samples_path = main_samples_path,
                                      num_cores = processing_cores,
                                      rw_disk = rw_disk,
                                      overwrite = overwrite)

annotated_SR_details <- junctionAnnotation(all_reads_combined = all_reads_combined, 
                                           main_samples_path = main_samples_path,
                                           blacklist_path = blacklist_path,
                                           gtf_path = gtf_path,
                                           bedtools_path = bedtools_path,
                                           fasta_path = fasta_path,
                                           fordownload_path = fordownload_path,
                                           rw_disk = rw_disk,
                                           overwrite = overwrite)

# 4. Pairing and distances calculation ----
all_distances_raw <- juctionPairing(metadata = metadata_RBPs,
                                    main_samples_path = main_samples_path,
                                    annotated_SR_details = annotated_SR_details,
                                    all_reads_combined = all_reads_combined,
                                    num_cores = 2,
                                    rw_disk = rw_disk,
                                    overwrite = overwrite)

all_distances_pruned <- removeAmbiguousPairing(all_distances_raw,
                                               main_samples_path,
                                               rw_disk = rw_disk,
                                               overwrite = overwrite)

# 5. Clustering of the sample ----
r <- foreach(i = seq(target_RBPs)) %do%{
  target_RBP <- target_RBPs[i]
  
  ## Target RBP variables
  RBP_metadata <- metadata_RBPs %>% 
    dplyr::filter(target_gene == target_RBP) %>% 
    dplyr::arrange(experiment_type)
  RBP_path <- paste0(main_samples_path, target_RBP, "/")
  RBP_clusters <- RBP_metadata %>%
    dplyr::pull(experiment_type) %>% unique()
  
  foreach(j = seq(RBP_clusters)) %do%{
    cluster_name <- RBP_clusters[j]
    cluster_path <- paste0(RBP_path, cluster_name, "/")
    
    cluster_samples <- RBP_metadata %>% 
      dplyr::filter(experiment_type == cluster_name) %>% 
      dplyr::pull(sample_id)
    
    ## Filter the general dataframes to contain only junctions from the sample
    ## reads
    logger::log_info("\t\t\t Loading Split reads and Annotated SR details for the cluster.")
    cluster_split_reads <- all_reads_combined %>%
      dplyr::select(junID, all_of(cluster_samples)) %>%
      dplyr::filter(if_any(all_of(cluster_samples),  ~ !is.na(.)))
    cluster_annotated_SR_details <- annotated_SR_details %>% 
      dplyr::filter(junID %in% (cluster_split_reads %>% dplyr::pull(junID)))
    cluster_distances_raw <- all_distances_raw %>%
      filter(sample %in% cluster_samples)
    
    ## Tidy the distances
    logger::log_info("\t\t\t Loading pruned distances for the cluster.")
    cluster_distances_pruned <- all_distances_pruned %>%
      dplyr::filter(sample %in% cluster_samples) %>%
      dplyr::distinct(novel_junID, ref_junID, .keep_all = T) %>%
      dplyr::select(-sample, -ref_counts, -novel_counts) %>%
      dplyr::relocate(novel_junID, ref_junID)
    
    cluster_distances_tidy <- addCounts(cluster_samples,
                                        cluster_distances_pruned,
                                        cluster_split_reads)
    
    ## Combine the never mis-spliced with the tidy distances
    cluster_distances_tidy_all <- addNeverMissplicedJunction(cluster_samples,
                                                             cluster_split_reads,
                                                             cluster_distances_raw,
                                                             cluster_annotated_SR_details,
                                                             cluster_distances_tidy,
                                                             cluster_path,
                                                             cluster_name,
                                                             rw_disk = F,
                                                             overwrite = overwrite)
    
    ## Generate the Database
    generateDB(cluster_distances_tidy_all,
               cluster_annotated_SR_details,
               cluster_path = cluster_path,
               cluster_name = cluster_name,
               u12_introns_path = u12_introns_path,
               u2_introns_path = u2_introns_path,
               rw_disk = rw_disk,
               overwrite = overwrite)
  }
}
