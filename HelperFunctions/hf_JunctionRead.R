#' Loads the reference genome into memory
#'
#' The different versions can be downloaded
#' \href{http://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/}{here}.
#'
#' @param gtf_path path to the reference genome .gtf file.
#'
#' @return the connection to the reference genome DB.
#' @export
loadEdb <- function(gtf_path) {
  if (!exists("edb")) {
    logger::log_info("\t\t Loading the v105 reference genome.")
    edb <<- ensembldb::ensDbFromGtf(gtf_path, outfile = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
    edb <<- ensembldb::EnsDb(x = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
  } else {
    logger::log_info("\t\t Variable 'edb' already loaded!")
  }
  
  return(edb)
}

#' Loads the ENCODE blacklisted regions into memory
#'
#' The different versions can be downloaded
#' \href{https://github.com/Boyle-Lab/Blacklist/tree/master/lists}{here}.
#'
#' @param blacklist_path path to the ENCODE blacklisted regions .bed file.
#'
#' @return the ENCODE blacklisted region in GRanges object.
#' @export
loadEncodeBlacklist <- function(blacklist_path) {
  if (!exists("encode_blacklist_hg38")) {
    logger::log_info("\t\t Loading the v2 ENCODE blacklisted regions.")
    encode_blacklist_hg38 <<- rtracklayer::import(blacklist_path) %>% diffloop::rmchr()
  } else {
    logger::log_info("\t\t Variable 'encode_blacklist_hg38' is already loaded!")
  }
  
  return(encode_blacklist_hg38)
}

junctionReading <- function(metadata,
                            main_samples_path,
                            num_cores = 4,
                            rw_disk = T,
                            overwrite = F){
  logger::log_info("\t Starting the junction reading process.")
  
  ## If rw_disk argument is set to TRUE and overwrite is set to FALSE, it looks
  ## for the results of the function in disk. If they have already been
  ## generated, the function is not executed and the files are returned.
  if (rw_disk & !overwrite & file.exists(paste0(main_samples_path, "all_reads_combined.rds"))) {
    logger::log_info("\t\t Ignoring junction reading.") 
    logger::log_info("\t\t File ", main_samples_path, "all_reads_combined.rds already exists.")
    logger::log_info("\t\t Loading file from disk.")
    return(readRDS(paste0(main_samples_path, "all_reads_combined.rds")))
  }
  
  ############ Read the junctions ############
  ##
  ## 1. Loops through every .junc file located in control and cases.
  ## 2. Reads the columns and apply transformations
  ## 3. Join all the information together
  
  ## Multiprocessing generation. Add argument "output.file" to
  ## parallel::makeCluster() if you want to output information about the
  ## parallel execution.
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  
  ## Multiprocessing loop
  logger::log_info("\t\t Reading all extracted BAM files (num_cores = ", num_cores, ").")
  all_junc <- foreach(i = 1:nrow(metadata), .packages = c("tidyverse")) %dopar% {
    ## Definition of the variables
    target_gene <- metadata[i, ] %>% pull(target_gene)
    sample_id <- metadata[i, ] %>% pull(sample_id)
    cluster <- metadata[i, ] %>% pull(experiment_type) %>% tolower()
    junc_path <- paste0(main_samples_path, target_gene, "/", cluster, "/", sample_id, ".bam.sort.s0.junc")
    
    if(!file.exists(junc_path)) return(tibble())
    
    ## Read the junction file into a tibble using readr::read_table()
    ## The "locale" argument is to read comma separated values (i.e. 25,12)
    junc <- readr::read_table(
      junc_path,
      col_names = F,
      show_col_types = F,
      progress = F,
      locale = readr::locale(grouping_mark = "")
    )
    
    ## Transformations of the junctions
    junc <- junc %>%
      dplyr::select(
        chr = X1,
        start = X2,
        stop = X3,
        junID = X4,
        reads = X5,
        strand = X6,
        blockSizes = X11
      ) %>%
      dplyr::mutate(strand = ifelse(strand == "?", "*", strand)) %>%
      dplyr::mutate(sampleID = sample_id) %>%
      tidyr::separate(col = blockSizes, sep = ",", c("blockSizesStart", "blockSizesEnd"), conver = T) %>%
      dplyr::mutate(start = start + blockSizesStart + 1, stop = stop - blockSizesEnd) %>%
      GenomicRanges::GRanges() %>%
      diffloop::rmchr() %>%
      GenomeInfoDb::keepSeqlevels(value = c(seq(1, 22) %>% as.character(), "X", "Y"), pruning.mode = "tidy") %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        seqnames = as.character(seqnames),
        strand = as.character(strand)
      ) %>%
      dplyr::select(-junID, -blockSizesStart, -blockSizesEnd)
    
    return(junc)
  }
  ## Stop the parallel cluster
  parallel::stopCluster(cl)
  
  ## Combine the resulted list into a single dataframe
  all_junc <- all_junc %>%
    dplyr::bind_rows() %>% 
    dplyr::distinct()
  invisible(gc())
  
  ############ Group and name the junctions ############
  ##
  ## 1. Name every junction by its seqname, start and end.
  ## 2. Generates a dataframe with the number of reads per junction and sample.
  ## 3. Store dataframe in disk.
  logger::log_info("\t\t Grouping the junctions across samples.")
  
  ## Generate the name of the junction by its seqname, start and end
  all_junc <- all_junc %>%
    dplyr::group_by(seqnames, start, end) %>%
    dplyr::mutate(junID = dplyr::cur_group_id(), .before = seqnames) %>%
    dplyr::ungroup()
  
  ## Generate the reads table
  all_reads_combined <- all_junc %>%
    tidyr::pivot_wider(values_from = reads, names_from = sampleID)
  
  ## Save to disk only if rw_disk is set to TRUE
  if (rw_disk) {
    logger::log_info("\t\t Saving output to ", main_samples_path, "all_reads_combined.rds")
    all_reads_combined %>% saveRDS(paste0(main_samples_path, "all_reads_combined.rds"))
  }
  
  logger::log_info("\t\t Reading process completed successfully.")
  return(all_reads_combined)
}

junctionAnnotation <- function(all_reads_combined,
                               main_samples_path,
                               blacklist_path = "/home/grocamora/RytenLab-Research/Additional_files/hg38-blacklist.v2.bed",
                               gtf_path = "/home/grocamora/RytenLab-Research/Additional_files/Homo_sapiens.GRCh38.105.chr.gtf",
                               bedtools_path = "/home/grocamora/tools/bedtools/",
                               fasta_path = "/home/grocamora/RytenLab-Research/Additional_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
                               fordownload_path = "/home/grocamora/tools/fordownload/",
                               rw_disk = T,
                               overwrite = F){
  logger::log_info("\t Starting the annotation process.")
  
  ## If rw_disk argument is set to TRUE and overwrite is set to FALSE, it looks
  ## for the results of the function in disk. If they have already been
  ## generated, the function is not executed and the files are returned.
  if (rw_disk & !overwrite & file.exists(paste0(main_samples_path, "annotated_SR_details.rds"))) {
    logger::log_info("\t\t Ignoring annotation process.") 
    logger::log_info("\t\t File ", main_samples_path, "annotated_SR_details.rds already exists.")
    logger::log_info("\t\t Loading file from disk.")
    return(readRDS(paste0(main_samples_path, "annotated_SR_details.rds")))
  }
  
  all_reads <- all_reads_combined %>%
    dplyr::select(junID, seqnames, start, end, width, strand) %>%
    GenomicRanges::GRanges()
  
  ############ Junction annotation ############
  ##
  ## Every step is divided into its own function, so that the process can be
  ## easily understood
  ##
  ## 1. Remove every junction from the ENCODE blacklisted regions.
  ## 2. Use the Dasper package to annotate and categorize every junction.
  ## 3. Remove junctions shorter than 25bp.
  ## 4. Remove junctions not categorized as novel_donor, novel_acceptor or 
  ##    annotated
  ## 5. Remove the junctions from ambiguous genes.
  ## 6. Generate the biotype percentage.
  ## 7. Generate the MaxEntScore.
  ## 8. Save the results.
  
  ## Remove ENCODE blacklisted regions
  encode_blacklist_hg38 <- loadEncodeBlacklist(blacklist_path)
  all_reads_pruned <- removeEncodeBlacklistRegions(all_reads, encode_blacklist_hg38)
  
  ## Annotate Dasper
  ## print(gtf_path)
  edb <- loadEdb(gtf_path)
  annotated_SR_details <- annotateDasper(all_reads_pruned, edb) %>% 
    tibble::as_tibble()
  
  ## Remove junctions shorter than 25bp
  annotated_SR_details <- removeShortJunctions(annotated_SR_details)
  
  ## Remove uncategorized junctions
  annotated_SR_details <- removeUncategorizedJunctions(annotated_SR_details)
  
  ## Remove junctions with ambiguous genes
  annotated_SR_details <- removeAmbiguousGenes(annotated_SR_details)
  
  ## Add biotype percentage
  annotated_SR_details <- generateBiotypePercentage(annotated_SR_details,
                                                    edb)
  ## Generate MaxEntScore
  annotated_SR_details <- generateMaxEntScore(annotated_SR_details,
                                              bedtools_path = bedtools_path,
                                              fasta_path = fasta_path,
                                              fordownload_path = fordownload_path)
  
  
  ## Save the results.
  ## Recommended name: annotated_SR_details.rds
  if (rw_disk) {
    logger::log_info("\t\t Saving output to ", main_samples_path, "annotated_SR_details.rds")
    annotated_SR_details %>% saveRDS(paste0(main_samples_path, "annotated_SR_details.rds"))
  }
  
  logger::log_info("\t\t Junction annotation executed successfully.")
  rm(all_reads, all_reads_pruned)
  invisible(gc())
  
  return(annotated_SR_details)
}



#' Remove the junctions from the ENCODE blacklisted regions
#'
#' @param GRdata GRanges class object with the relevant junctions.
#' @param encode_blacklist_hg38 GRanges class object with the blacklisted
#'   regions.
#'
#' @return junctions that do not overlap with the blacklisted regions.
#' @export
removeEncodeBlacklistRegions <- function(GRdata,
                                         encode_blacklist_hg38) {
  logger::log_info("\t\t Removing junctions from ENCODE blacklisted regions.")
  ## Look fot the overlaps between GRdata and the ENCODE blacklisted region
  overlaps <- GenomicRanges::findOverlaps(
    query = encode_blacklist_hg38,
    subject = GRdata,
    ignore.strand = F,
    type = "any"
  )
  
  idxs <- S4Vectors::subjectHits(overlaps)
  
  ## If an overlap is found, remove the junctions
  if (length(idxs) > 0) {
    #logger::log_info(paste0("A total of ", length(unique(idxs)), " junctions overlap with an ENCODE blacklisted region."))
    GRdata <- GRdata[-idxs, ]
    #logger::log_info(paste0("A total of ", length(GRdata), " are left after the removal."))
  } else {
    logger::log_info("No junctions overlapping with an ENCODE blacklist region.")
  }
  
  return(GRdata)
}

#' Annotates and categorize every junction
#'
#' @param GRdata GRanges class object with the relevant junctions.
#' @param edb The connection to the reference genome DB.
#'
#' @return Annotated input by dasper.
#' @export 
annotateDasper <- function(GRdata, edb) {
  logger::log_info("\t\t Annotating using dasper::junction_annot().")
  GRdata <- dasper::junction_annot(GRdata,
                                   ref = edb,
                                   ref_cols = c("gene_id", "gene_name", "symbol", "tx_id"),
                                   ref_cols_to_merge = c("gene_id", "gene_name", "tx_id")
  )
  
  return(GRdata)
}

#' Remove the junctions shorter than 25bp.
#'
#' @param input_SR_details Tibble (or data.frame) bject with the relevant
#'   junctions.
#'
#' @return Junctions bigger than 25bp.
#' @export
removeShortJunctions <- function(input_SR_details) {
  logger::log_info("\t\t Removing junctions smaller than 25bp.")
  output_SR_details <- input_SR_details %>%
    dplyr::filter(width >= 25)
  
  return(output_SR_details)
}

#' Remove the junctions with uncategorized classification.
#'
#' @param input_SR_details Tibble (or data.frame) bject with the relevant
#'   junctions.
#'
#' @return Junctions categorized as "annotated", "novel_donor" or
#'   "novel_acceptor".
#' @export
removeUncategorizedJunctions <- function(input_SR_details) {
  logger::log_info("\t\t Removing junctions not classified as novel_acceptor, novel_donor or annotated.")
  output_SR_details <- input_SR_details %>%
    dplyr::filter(type %in% c("novel_acceptor", "novel_donor", "annotated"))
  
  return(output_SR_details)
}

#' Remove the junctions from ambiguous genes.
#'
#' @param input_SR_details Tibble (or data.frame) bject with the relevant
#'   junctions.
#'
#' @return Junctions assigned to only one gene.
#' @export
removeAmbiguousGenes <- function(input_SR_details) {
  logger::log_info("\t\t Removing junctions associated to more than one gene.")
  output_SR_details <- input_SR_details[which(sapply(input_SR_details$gene_id_junction, length) == 1), ]
  
  return(output_SR_details)
}

#' Executes the MaxEntScan
#'
#' @param input_SR_details Tibble (or data.frame) bject with the relevant
#'   junctions.
#' @param bedtools_path path to the
#'   \href{https://bedtools.readthedocs.io/en/latest/}{bedtools} executable. Can
#'   be left empty if bedtools is in default PATH.
#' @param fasta_path path to the fasta .fa file for the reference genome.
#' @param fordownload_path path to the MaxEntScan pearl scripts. Can be
#'   downloaded from
#'   \href{http://hollywood.mit.edu/burgelab/maxent/download/}{here}.
#'
#' @return Junction dataframe with the ss5 and ss3 scores.
#' @export
generateMaxEntScore <- function(input_SR_details,
                                bedtools_path = "~/tools/bedtools/",
                                fasta_path = "~/RytenLab-Research/Additional_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
                                fordownload_path = "~/tools/fordownload/") {
  logger::log_info("\t\t Calculating the MaxEntScore.")
  current_wd <- getwd()
  
  ## Add extra positions to calculate the MaxEntScore
  GRdata <- input_SR_details %>% 
    dplyr::mutate(
      seqnames = seqnames %>% as.character(),
      donorSeqStart = ifelse(strand == "-", end - 6, start - 4),
      donorSeqStop = ifelse(strand == "-", end + 3, start + 5),
      acceptorSeqStart = ifelse(strand == "-", start - 4, end - 20),
      acceptorSeqStop = ifelse(strand == "-", start + 19, end + 3)
    )
  
  ## Generate the dataframes for both donor and acceptors.
  donor_df <- GRdata %>%
    dplyr::select(seqnames,
                  starts = donorSeqStart,
                  end = donorSeqStop,
                  names = junID,
                  strands = strand
    ) %>%
    dplyr::mutate(scores = ".", .before = strands)
  
  acceptor_df <- GRdata %>%
    dplyr::select(seqnames,
                  starts = acceptorSeqStart,
                  end = acceptorSeqStop,
                  names = junID,
                  strands = strand
    ) %>%
    dplyr::mutate(scores = ".", .before = strands)
  
  ## Generate temporary files to store the sequences
  tmp_file <- tempfile()
  tmp_file_seq <- tempfile()
  
  ############ MaxEntScore ############
  ##
  ## 1. For both donor and acceptor dataframes, extract the sequences of the
  ## junctions and add them to the main dataframe.
  ## 2. Calculate the ss5 and ss3 scores (donor and acceptor sequences).
  ## 3. Add the scores to the main dataframe.
  
  ## Donor sequence
  data.table::fwrite(donor_df, file = tmp_file, quote = F, sep = "\t", row.names = F, col.names = F, scipen = 50)
  rm(donor_df)
  system2(
    command = paste0(bedtools_path, "bedtools"),
    args = c(
      "getfasta -nameOnly -s",
      "-fi", fasta_path,
      "-bed", tmp_file,
      "-tab -fo", tmp_file_seq
    )
  )
  donor_sequences_input <- read.delim(tmp_file_seq, header = F)
  stopifnot(identical(
    gsub("\\(\\+\\)", "", gsub("\\(\\*\\)", "", gsub("\\(-\\)", "", as.character(donor_sequences_input$V1)))),
    GRdata$junID %>% as.character()
  ))
  
  ## Acceptor sequence
  data.table::fwrite(acceptor_df, file = tmp_file, quote = F, sep = "\t", row.names = F, col.names = F, scipen = 50)
  rm(acceptor_df)
  system2(
    command = paste0(bedtools_path, "bedtools"),
    args = c(
      "getfasta -nameOnly -s",
      "-fi", fasta_path,
      "-bed", tmp_file,
      "-tab -fo", tmp_file_seq
    )
  )
  acceptor_sequences_input <- read.delim(tmp_file_seq, header = F)
  stopifnot(identical(
    gsub("\\(\\+\\)", "", gsub("\\(\\*\\)", "", gsub("\\(-\\)", "", as.character(acceptor_sequences_input$V1)))),
    GRdata$junID %>% as.character()
  ))
  
  ## Add to dataframe
  GRdata <- GRdata %>%
    dplyr::mutate(donor_sequence = as.character(donor_sequences_input$V2)) %>%
    dplyr::mutate(acceptor_sequence = as.character(acceptor_sequences_input$V2))
  
  ## Generate MaxEntScore
  tmp_file <- tempfile()
  setwd(fordownload_path)
  
  data.table::fwrite(list(GRdata$donor_sequence), file = tmp_file, quote = F, row.names = F, col.names = F, scipen = 50)
  ss5score_vector <- read.delim(pipe(paste0("perl ", fordownload_path, "score5.pl ", tmp_file)), header = F)
  GRdata <- GRdata %>% dplyr::mutate(ss5score = ss5score_vector$V2) 
  rm(ss5score_vector)
  
  data.table::fwrite(list(GRdata$acceptor_sequence), file = tmp_file, quote = F, row.names = F, col.names = F, scipen = 50)
  ss3score_vector <- read.delim(pipe(paste0("perl ", fordownload_path, "score3.pl ", tmp_file)), header = F)
  GRdata <- GRdata %>% dplyr::mutate(ss3score = ss3score_vector$V2)
  rm(ss3score_vector)
  
  setwd(current_wd)
  rm(tmp_file, tmp_file_seq, donor_sequences_input, acceptor_sequences_input)
  invisible(gc())
  return(GRdata)
}

#' Calculates the biotype percentage
#'
#' @param input_SR_details Tibble (or data.frame) bject with the relevant
#'   junctions.
#' @param edb The connection to the reference genome DB.
#'
#' @return Annotated junctions with their protein coding and lncRNA biotypes
#'   percentage added.
#' @export
generateBiotypePercentage <- function(input_SR_details,
                                      edb){
  logger::log_info("\t\t Calculating the biotype percentage.")
  ## Extract the reference transcriptome
  transcript_v105 <- ensembldb::transcripts(edb) %>% 
    tibble::as_tibble() %>% 
    dplyr::select(tx_id, tx_biotype, tx_gene_id = gene_id)
  
  ## Merge the transcript biotype to the intron data.frame
  df_all_introns <- input_SR_details %>% 
    dplyr::select(junID, tx_id_junction, gene_id_junction) %>%
    tidyr::unnest(tx_id_junction) %>%
    dplyr::left_join(transcript_v105,
                     by = c("tx_id_junction" = "tx_id"))
  
  ## Calculate percentage of unambiguous junctions (no multiple genes)
  df_all_percentage <- df_all_introns %>%
    dplyr::group_by(junID) %>%
    dplyr::filter(n_distinct(tx_gene_id) == 1)  %>%
    dplyr::group_by(junID, tx_biotype) %>%
    dplyr::distinct(tx_id_junction, .keep_all = T) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(percent = (n/sum(n))*100) %>%
    dplyr::ungroup()
  
  ## Filter only protein coding and lncRNA
  df_all_percentage_merged <- df_all_percentage %>% 
    dplyr::group_by(junID) %>%
    tidyr::pivot_wider(junID, names_from = tx_biotype, values_from = percent) %>%
    dplyr::select(junID, percent_PC = protein_coding, percent_lncRNA = lncRNA) %>%
    replace(is.na(.), 0) %>%
    dplyr::ungroup()
  
  ## Merge biotype percentage with the annotated SR details
  output_SR_details <- input_SR_details %>%
    dplyr::left_join(y = df_all_percentage_merged,
                     by = "junID") %>%
    dplyr::rename(protein_coding = percent_PC,
                  lncRNA = percent_lncRNA)
  
  return(output_SR_details)
}
