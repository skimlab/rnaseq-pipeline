library("tidyverse")
library("rjson")
library("tximport")

collect_meta_data <- function(a_salmon_folder) {
  
  lib_format <-
    fromJSON(file = file.path(a_salmon_folder, "lib_format_counts.json"))
  
  lib_format <-
    lib_format[c(
      "num_frags_with_concordant_consistent_mappings",
      "num_frags_with_inconsistent_or_orphan_mappings",
      "strand_mapping_bias"
    )]
  
  meta_info <-
    fromJSON(file = file.path(a_salmon_folder, "aux_info", "meta_info.json"))
  
  meta_info <-
    meta_info[c(
      "salmon_version",
      "opt_type",
      "frag_dist_length",
      "frag_length_mean",
      "frag_length_sd",
      "seq_bias_correct",
      "gc_bias_correct",
      "num_bias_bins",
      "mapping_type",
      "keep_duplicates",
      "num_valid_targets",
      "num_decoy_targets",
      "num_eq_classes",
      "num_bootstraps",
      "num_processed",
      "num_mapped",
      "num_decoy_fragments",
      "num_dovetail_fragments",
      "num_fragments_filtered_vm",
      "num_alignments_below_threshold_for_mapped_fragments_vm",
      "percent_mapped",
      "call",
      "start_time",
      "end_time"
    )]
  
  c(lib_format, meta_info)
}

split_gtf_names <- function(gtf_names, indices_to_keep = 1:7) {
  gtf_annots_list <- strsplit(gtf_names, split = "|", fixed = TRUE)
  
  lapply(gtf_annots_list,
         function(x) {
           xn <-
             c(x[indices_to_keep], paste(x[-indices_to_keep], collapse = "|"))
           names(xn) <-
             c(
               "ENST_ID",
               "ENSG_ID",
               "HAVANA_GENE",
               "HAVANA_TRANSCRIPT",
               "TRANSCRIPT_ID",
               "GENE_ID",
               "TRANSCRIPT_LENGTH",
               "OTHERS"
             )[c(indices_to_keep, 8)]
           xn
         })
}


get_sample_folders <- function(base_folder) {
  sample_ids <- dir(base_folder)
  sample_folders <- file.path(base_folder, sample_ids)
  names(sample_folders) <- sample_ids
  
  sample_folders
}

get_sample_files <- function(sample_folders, filename = "quant.sf") {
  sample_files <- file.path(sample_folders, filename)
  if (!is.null(names(sample_files))) {
    names(sample_files) <- names(sample_folders)
  } else {
    names(sample_files) <- basename(sample_folders)
  }
  
  sample_files
}

aggregate_meta_data <- function(sample_folders) {
  lapply(sample_folders,
         collect_meta_data) %>%
    bind_rows(.id = "sample_ids")
}


aggregate_salmon_files <- function(sample_folders) {
  salmon_files <- get_sample_files(sample_folders, filename = "quant.sf")

  txi <- tximport(salmon_files, type = "salmon", txOut = TRUE)

  txi$rowData <-
    split_gtf_names(gtf_names = rownames(txi$abundance)) %>% bind_rows()

  for (n in names(txi)) {
    if ( ("matrix" %in% class(txi[[n]])) && (nrow(txi[[n]]) == nrow(txi$rowData))) {
      rownames(txi[[n]]) <- txi$rowData[[1]]
    }
  }

  txi$colData <-
    aggregate_meta_data(sample_folders)
  
  txi
}


summarizeToGene <- function(object, tx2gene) {
  object_gene <- 
    object %>% tximport::summarizeToGene(tx2gene = tx2gene)
  
  object_gene[["colData"]] <- object[["colData"]]
    
  rowData <- data.frame(rownames(object_gene[["abundance"]]))
  names(rowData)[1] <- names(tx2gene)[2]
  object_gene[["rowData"]] <- rowData
  
  object_gene
}

save_aggregated_ouptputs <- function(object, file_base = "salmon") {
  write_csv(as_tibble(object[["abundance"]], rownames = "GENE_ID"), file = sprintf("%s_abundance.csv", file_base))
  write_csv(as_tibble(object[["counts"]], rownames = "GENE_ID"), file = sprintf("%s_counts.csv", file_base))
  write_csv(object[["colData"]], file = sprintf("%s_sample_info.csv", file_base))
  write_csv(object[["rowData"]], file = sprintf("%s_annots.csv", file_base))
  write_csv(as_tibble(object[["length"]], rownames = "GENE_ID"), file = sprintf("%s_estimated_length.csv", file_base))
}


# folder_list <- c(
#   "FASTQ_Generation_2021-10-06_16_42_24Z-469231862/outs/",
#   "FASTQ_Generation_2021-10-15_19_59_52Z-473785312/outs"
# )

folder_list <- c(
  "FASTQ_Generation_2021-10-06_16_42_24Z-469231862/outs_bc/",
  "FASTQ_Generation_2021-10-15_19_59_52Z-473785312/outs_bc"
)

salmon_transcript_list <-
  lapply(folder_list,
         get_sample_folders) %>%
  unlist() %>%
  aggregate_salmon_files()


salmon_gene_list <-
  salmon_transcript_list %>% summarizeToGene(tx2gene = salmon_transcript_list$rowData[c("ENST_ID", "GENE_ID")])


save_aggregated_ouptputs(salmon_transcript_list, file_base = file.path("jhu_2021-10-26_outs_bc_aggregated", "jhu_2021-10-26_transcript"))

save_aggregated_ouptputs(salmon_gene_list, file_base = file.path("jhu_2021-10-26_outs_bc_aggregated", "jhu_2021-10-26_gene"))

