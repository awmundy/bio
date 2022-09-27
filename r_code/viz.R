suppressPackageStartupMessages({
  library(tidyverse)
  library(tximport)
  library(ensembldb)
  library(EnsDb.Mmusculus.v79)
  library(EnsDb.Hsapiens.v86)
  library(edgeR)
  library(matrixStats)
  library(cowplot)
  library(zeallot)
}) 

get_abundance_paths_old <- function(sra_accessions, abundance_root_dir) {
  # extract sra_accession values (which are the abundance subdir names)
  seq_folders <- unlist(study_design$sra_accession)
  
  # construct abundance paths
  abundance_paths <- c()
  for(i in 1:length(seq_folders))
    abundance_paths[i] <- paste(abundance_root_dir, 
                                seq_folders[i], 
                                '/abundance.tsv', 
                                sep='')
  
  stopifnot(all(file.exists(abundance_paths))) 
  return(abundance_paths)
}

get_abundance_paths <- function(abundance_root_dir) {
  abundance_dirs <- list.dirs(abundance_root_dir)[-1]
  
  abundance_paths <- c()
  for(i in 1:length(abundance_dirs))
    abundance_paths[i] <- paste(normalizePath(abundance_dirs[i]), 
                                '/abundance.tsv', 
                                sep='')
  
  stopifnot(all(file.exists(abundance_paths)))
  return(abundance_paths)
}

get_transcript_to_gene_df <- function(annotation_db){
  # query the ensembl db to get gene transcripts id: gene name map
  tx <- transcripts(annotation_db, columns=c("tx_id", "gene_name"))
  # convert to tibble and clean up
  tx <- as_tibble(tx)
  tx <- dplyr::rename(tx, target_id = tx_id)
  tx <- dplyr::select(tx, "target_id", "gene_name")
  return(tx)
}

add_row_descriptive_stats <- function(df, col_subset) {
  df$sd <- rowSds(as.matrix(abund[, col_subset]))
  df$mean <- rowMeans(as.matrix(abund[, col_subset]))
  df$median <- rowMedians(as.matrix(abund[, col_subset]))
  return(df)
}

convert_tx_gene_mtx_to_tibble <- function(mtx, sample_labels) {
  # add column names that are the names of the samples
  stopif(is.null(sample_labels))
  colnames(mtx) <- sample_labels
  df <- as_tibble(mtx, rownames = 'gene')
  return(df)
}

build_digital_gene_expression_list <- function(gene_counts, sample_labels) {
  # @return A list containing a matrix of counts (same as what was passed in) 
  #         and a dataframe of samples containing the total count and a 
  #         normalization factor of 1
  
  # subset to sample labels cols only, and convert to matrix
  dge_list <- DGEList(as.matrix(gene_counts[, sample_labels]))  
  return(dge_list)
}

build_log_cpm_df_old <- function(dge_list, long) {
  # Construct a logged counts per million df
  
  log_cpm <- cpm(dge_list, log=TRUE)
  log_cpm <- as_tibble(log_cpm, rownames = "gene_id")
  # convert gene_id to int for later libraries that require it
  log_cpm <- as_tibble(transform(log_cpm, gene_id=as.numeric(gene_id)))
  if (long == TRUE) {
    log_cpm <- pivot_longer(log_cpm,
                            cols = HS01:CL13,
                            names_to = "samples",
                            values_to = "expression")
  }
  return(log_cpm)
}

build_log_cpm_df <- function(dge_list, long) {
  # Construct a logged counts per million df
  
  log_cpm <- cpm(dge_list, log=TRUE)
  log_cpm <- as_tibble(log_cpm, rownames = "gene_id")
  # convert gene_id to int for later libraries that require it
  log_cpm <- as_tibble(transform(log_cpm, gene_id=as.numeric(gene_id)))
  if (long == TRUE) {
    sample_cols <- colnames(log_cpm)
    sample_cols <- all_of(sample_cols[-1])
    log_cpm <- pivot_longer(log_cpm,
                            cols = sample_cols,
                            names_to = "sample_label",
                            values_to = "expression")
  }
  return(log_cpm)
}

build_log_cpm_plot_old <- function(cpm_df, subtitle) {
  plt <- ggplot(cpm_df) +
    aes(x=samples, y=expression, fill=samples) +
    geom_violin(trim = FALSE, show.legend = FALSE) +
    stat_summary(fun = "median", 
                 geom = "point", 
                 shape = 95, 
                 size = 10, 
                 color = "black", 
                 show.legend = FALSE) +
    labs(y="log2 expression", x = "sample",
         title="Log2 Counts per Million (CPM)",
         subtitle=subtitle,
         caption=paste0("produced on ", Sys.time())) +
    theme_bw()
  return(plt)
}

build_log_cpm_plot <- function(cpm_df, subtitle) {
  stopifnot('sample_label' %in% colnames(cpm_df))
  
  plt <- ggplot(cpm_df) +
    aes(x=sample_label, y=expression, fill=sample_label) +
    geom_violin(trim = FALSE, show.legend = FALSE) +
    stat_summary(fun = "median", 
                 geom = "point", 
                 shape = 95, 
                 size = 10, 
                 color = "black", 
                 show.legend = FALSE) +
    labs(y="log2 expression", x = "sample",
         title="Log2 Counts per Million (CPM)",
         subtitle=subtitle,
         caption=paste0("produced on ", Sys.time())) +
    theme_bw()
  return(plt)
}

filter_dge_list <- function(dge_list, min_cpm, min_samples_with_min_cpm) {
  # get a counts per million matrix where each record is a gene and each col
  #   is a sample
  cpm_mtx <- cpm(dge_list)
  # filter dge_list 
  msk <- rowSums(cpm_mtx > min_cpm) >= min_samples_with_min_cpm
  dge_list <- dge_list[msk, ]
  return(dge_list)
}

get_gene_level_stats_dfs <- function(study_design, tx_to_gene_df) {
  # params
  # abundance_paths: c vector of paths of abundance files
  # tx_to_gene_df: df mapping a transcript identifier to a gene name
  
  # get list of matrixes summarizing information across the data in abundance_paths
  gene_stats_list <- tximport(file = study_design$abundance_paths,
                              type = "kallisto",
                              tx2gene = tx_to_gene_df,
                              txOut = FALSE, # false means gene level, not tx level
                              countsFromAbundance = "lengthScaledTPM",
                              ignoreTxVersion = TRUE,
                              abundanceCol = 'tpm')
  sample_labels <- study_design$sample_label
  
  # add column headers and convert to tibble
  gene_counts <- convert_tx_gene_mtx_to_tibble(gene_stats_list$counts, 
                                               sample_labels)
  gene_lengths <- convert_tx_gene_mtx_to_tibble(gene_stats_list$length, 
                                               sample_labels)
  gene_abunds <- convert_tx_gene_mtx_to_tibble(gene_stats_list$abundance, 
                                               sample_labels)
  
  return(list(gene_counts, gene_lengths, gene_abunds))
}


assign_abundance_paths_to_study_design_old <-
  function(study_design,
           abundance_paths) {
    abundance_df <- data.frame(abundance_paths)
    abundance_df$sample_label <-
      str_replace(abundance_df$abundance_paths, abundance_root_dir, '')
    abundance_df$sample_label <-
      str_replace(abundance_df$sample_label, '/abundance.tsv', '')
    
    study_design <-
      merge(study_design, abundance_df, by.x = 'sra_accession', 
            by.y='sample_label', all = TRUE)
    assert_col_not_null(study_design$abundance_paths)
    
    return(study_design)
  }

assign_abundance_paths_to_study_design <- function(study_design,
                                                   abundance_paths) {
  abundance_df <- data.frame(abundance_paths)
  abundance_df$sample_label <-
    str_replace(abundance_df$abundance_paths, abundance_root_dir, '')
  abundance_df$sample_label <-
    str_replace(abundance_df$sample_label, '/abundance.tsv', '')
  
  study_design <-
    merge(study_design, abundance_df, by = 'sample_label', all = TRUE)
  assert_col_not_null(study_design$abundance_paths)
  
  return(study_design)
}

assert_col_unique <- function(col) {
  stopifnot(length(unique(col)) == length(col))
}

assert_col_not_null <- function(col) {
  stopif(any(is.na(col)))
  }

get_study_design_df_old <- function(study_design_path) {
  study_design <- read_csv(study_design_path, col_types = cols(.default = 'c'))
  study_design <- dplyr::rename(study_design, sample_label = sample)
  
  return(study_design)
}

get_study_design_df <- function(study_design_path) {
  study_design <- read_csv(study_design_path, col_types = cols(.default = 'c'))
  assert_col_unique(study_design$sample_label)
  
  return(study_design)
}

stopif <- function(logic) {
  stopifnot(!logic)
}


# abundance_root_dir <- '/media/amundy/Windows/bio/diyt/rna_txs/fastq_folders/'
# study_design_path <- '/media/amundy/Windows/bio/diyt/studydesign.csv'
# study_design <- get_study_design_df_old(study_design_path)
# sample_labels <- study_design$sample_label
# sra_accessions <- study_design$sra_accession
# sample_groups <- study_design$group
# abundance_paths <- get_abundance_paths_old(study_design, abundance_root_dir)
# study_design <- assign_abundance_paths_to_study_design_old(study_design, abundance_paths)
# tx_to_gene_df <- get_transcript_to_gene_df(EnsDb.Hsapiens.v86)
# c(gene_counts, gene_lengths, gene_abunds) %<-% get_gene_level_stats_dfs(study_design, tx_to_gene_df)
# dge_list <- build_digital_gene_expression_list(gene_counts, study_design$sample_label)
# log_cpm_long <- build_log_cpm_df_old(dge_list, long = TRUE)




abundance_root_dir <- '/media/amundy/Windows/bio/ac_thymus/rna_txs/fastq_folders/'
study_design_path <- '/media/amundy/Windows/bio/ac_thymus/study_design_removed_bad_one.csv'


study_design <- get_study_design_df(study_design_path)
abundance_paths <- get_abundance_paths(abundance_root_dir)

study_design <- assign_abundance_paths_to_study_design(study_design, abundance_paths)

tx_to_gene_df <- get_transcript_to_gene_df(EnsDb.Mmusculus.v79)

c(gene_counts, gene_lengths, gene_abunds) %<-% get_gene_level_stats_dfs(study_design, tx_to_gene_df)

#TODO build some sort of relevant plot for these abundances
# gene_abunds <- add_row_descriptive_stats(gene_abunds, study_design$sample_label)

dge_list <- build_digital_gene_expression_list(gene_counts, study_design$sample_label)
log_cpm_long <- build_log_cpm_df(dge_list, long = TRUE)

dge_list_filt <- filter_dge_list(dge_list,
                                 min_cpm = 1,
                                 min_samples_with_min_cpm = 5)
log_cpm_filt_long <- build_log_cpm_df(dge_list_filt, long = TRUE)

dge_list_filt_norm <- calcNormFactors(dge_list_filt, method = 'TMM')
log_cpm_filt_norm_long <- build_log_cpm_df(dge_list_filt_norm, long = TRUE)

plt_1 <- build_log_cpm_plot(log_cpm_long, "unfiltered, non-normalized")
plt_2 <- build_log_cpm_plot(log_cpm_filt_long, "filtered, non-normalized")
plt_3 <- build_log_cpm_plot(log_cpm_filt_norm_long, "filtered, normalized")
plot_grid(plt_1, plt_2, plt_3, labels = c('A', 'B', 'C'), label_size = 12)

#TODO confirm this needs to be a factor
# sample_groups <- factor(sample_groups)
log_cpm_filt_norm <- build_log_cpm_df(dge_list_filt_norm, long = FALSE)

#TODO probably just remove this, these objects are weird and the plot doesnt even work
# distance <- dist(t(log_cpm_filt_norm), method = "maximum")
# clusters <- hclust(distance, method = "average") 
# plot(clusters, labels=sample_labels)

pca <- prcomp(t(log_cpm_filt_norm), scale.=F, retx=T)
summary(pca)
