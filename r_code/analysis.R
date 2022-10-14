suppressPackageStartupMessages({
	library(tidyverse)
	library(tximport)
	library(ensembldb)
	library(EnsDb.Mmusculus.v79)
	library(EnsDb.Hsapiens.v86)
	library(edgeR)
	library(matrixStats)
	library(gridExtra)
	library(zeallot)
	library(cowplot)
	library(rhdf5)
	library(limma)
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

  if (long == TRUE) {
    sample_cols <- colnames(log_cpm)[-1]
    log_cpm <- pivot_longer(log_cpm,
                            cols = all_of(sample_cols),
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
  msk <- rowSums(cpm_mtx >= min_cpm) >= min_samples_with_min_cpm
  dge_list <- dge_list[msk, ]
  return(dge_list)
}

get_gene_level_stats_dfs <- function(abundance_paths, sample_labels, tx_to_gene_df) {
	# params:
	# abundance_paths: c vector of paths of abundance files
	# tx_to_gene_df: df mapping a transcript identifier to a gene name
	# returns:
	# gene_counts: count of genes aggregated from txs- not exactly equal to 
	#			   a million- possibly bc some txs don't have a gene in the reference?
	# gene_lengths
	# gene_abunds
  
  # get list of matrixes summarizing information across the data in abundance_paths
  gene_stats_list <- tximport(file = abundance_paths,
                              type = "kallisto",
                              tx2gene = tx_to_gene_df,
                              txOut = FALSE, # false means gene level, not tx level
                              countsFromAbundance = "lengthScaledTPM",
                              ignoreTxVersion = TRUE,
                              abundanceCol = 'tpm')
  
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
  abundance_df <-  dplyr::rename(abundance_df, abundance_path = abundance_paths)
  abundance_df$sample_label <-
    str_replace(abundance_df$abundance_path, abundance_root_dir, '')
  abundance_df$sample_label <-
    str_replace(abundance_df$sample_label, '/abundance.tsv', '')
  
  study_design <-
    merge(study_design, abundance_df, by = 'sample_label', all = TRUE)
  assert_col_not_null(study_design$abundance_path)
  
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

convert_tibble_to_mtx <- function(tbl, row_names_col) {
	mtx <- as.matrix(tbl)
	
	# set rownames
	row.names(mtx) = mtx[, row_names_col]
	
	# remove row name col
	mtx <- mtx[, ! colnames(mtx) %in% c(row_names_col)]
	return(mtx)
}

write_cluster_dendogram_plot <- function(tbl, sample_labels, cluster_out_path) {
	## writes a dendogram cluster  plot for a gene level dataset where 
	##	the headers are sample_labels
	
	# build matrix since distance function requires is
	mtx <- convert_tibble_to_mtx(tbl, "gene_id")
	stopifnot(setequal(colnames(mtx), sample_labels))
	
	# calc distance between samples
	distance <- dist(t(mtx), method = "maximum")
	# agglomeration
	clusters <- hclust(distance, method = "average")
	
	# build dendogram and write
	pdf(cluster_out_path)
	plot(clusters, labels=sample_labels)
	dev.off()
}

write_pca_scatter_plots <- function(pca_metrics, sample_dimensions,
									study_design, pca_out_path) {
	sample_labels <- study_design$sample_label
	# pca_var_pct <- pca_metrics[[1]]
	# pca_loadings <- pca_metrics[[2]]
	plot_list = list()
	for(sample_dimension in sample_dimensions) {
		# build factor for coloring
		sample_dimension_factor <- factor(study_design[, sample_dimension])
		# length 1 means factor failed, e.g. doesn't work on tibbles
		stopif(length(sample_dimension_factor) == 1)
		
		# use just odd elements of range so that PCs aren't repeated across plots
		for(i in seq(1, 6, 2)) {
			first_pc_label = paste('PC', i, sep='')
			second_pc_label = paste('PC', i + 1, sep='')
			
			plt <- ggplot(pca_metrics$scores) +
				aes(x=get(first_pc_label), y=get(second_pc_label),
					label=sample_labels, color=sample_dimension_factor) +
				geom_point(size=4) +
				# geom_label(nudge_y = 10) +
				xlab(paste0(first_pc_label, "(", pca_metrics$var_explained[i]*100,"%",")")) +
				ylab(paste0(second_pc_label, "(", pca_metrics$var_explained[i + 1]*100,"%",")")) +
				labs(title="PCA plot",
					 caption=paste0("produced on ", Sys.time())) +
				coord_fixed() +
				# overwrites legend title to be the actual dimension label
				guides(color=guide_legend(sample_dimension)) +
			theme_bw()
			
			# print plt so that recordPlot can capture it
			print(plt)
			plot_list[[paste(i, sample_dimension, sep='')]] <- recordPlot()
		}}
	
	# write out
	pdf(pca_out_path, onefile=TRUE, width=4, height=4)
	for (plt in plot_list) {
		replayPlot(plt)
	}
	graphics.off()
}

get_pca_metrics <- function(cpm_matrix) {
	## Performs Principal Component Analysis on genetic counts per million
	##	sample data
	## Returns: 
	##	variance explained: numeric series like object with pct variance described
	##				 by each PC
	##	scores: tibble, shows how much each sample influenced each PC
	##  loadings: show much each variable (gene) influenced each PC
	
	#  [-1] removes gene_id column
	pca <- prcomp(t(cpm_matrix[-1]), scale.=F, retx=T)
	pca_sum <- summary(pca)
	var_explained <- as.data.frame(t(pca_sum$importance))$'Proportion of Variance'
	scores <- as_tibble(pca$x)
	loadings <- pca$rotation
	
	pca_metrics <- list(var_explained = var_explained,
						loadings = loadings,
						scores = scores
						)
	
	return(pca_metrics)
}

write_pca_small_multiples_plots <- function(pca_metrics,
											sample_dimensions,
											study_design,
											pca_small_multiples_out_path) {
	
	sample_labels <- study_design$sample_label
	scores <- pca_metrics$scores[, 1:6]
	
	plot_list = list()
	for (sample_dimension in sample_dimensions) {
		sample_dim_vals <- study_design[, sample_dimension]
		# prep and pivot long
		pca_pivot <- add_column(scores, sample = sample_labels,
								group = sample_dim_vals)
		pca_pivot <- pivot_longer(
			pca_pivot,
			cols = PC1:PC6,
			names_to = "PC",
			values_to = "scores"
		)
		
		# build small multiples plot
		plt <- ggplot(pca_pivot) +
			aes(x = sample,
				y = scores,
				fill = group) +
			geom_bar(stat = "identity") +
			facet_wrap( ~ PC) +
			labs(title = "PCA small multiples plot",
				 caption = paste0("produced on ", Sys.time())) +
			theme_bw() +
			coord_flip() +
			guides(color = guide_legend(sample_dimension))
		# print plt so that recordPlot can capture it
		print(plt)
		plot_list[[sample_dimension]] <- recordPlot()
	}
	
	pdf(pca_small_multiples_out_path, onefile = TRUE)
	for (plt in plot_list) {
		replayPlot(plt)
	}
	graphics.off()
}

describe <- function(df_col) {
	min_val <- min(df_col)
	max_val <- max(df_col)
	mean_val <- mean(df_col)
	q5 <- quantile(df_col, 0.05)
	q10 <- quantile(df_col, 0.10)
	q25 <- quantile(df_col, 0.25)
	q50 <- quantile(df_col, 0.50)
	q75 <- quantile(df_col, 0.75)
	q90 <- quantile(df_col, 0.90)
	q95 <- quantile(df_col, 0.95)
	out <- data.frame(min=min_val, max=max_val, mean=mean_val, 
					  q5=q5, q10=q10, q25=q25, q50=q50, q75=q75, 
					  q90=q90, q95=q95, row.names = NULL)
	out <- t(out)
	return(out)
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
pca_scatter_out_path <- "/home/amundy/Desktop/pca_scatter.pdf"
pca_small_multiples_out_path <- "/home/amundy/Desktop/pca_small_multiples.pdf"
pca_scatter_ext_out_path <- "/home/amundy/Desktop/pca_scatter_ext.pdf"
pca_small_multiples_ext_out_path <- "/home/amundy/Desktop/pca_small_multiples_ext.pdf"
cluster_out_path <- "/home/amundy/Desktop/cluster.pdf"

study_design <- get_study_design_df(study_design_path)
abundance_paths <- get_abundance_paths(abundance_root_dir)
study_design <- assign_abundance_paths_to_study_design(study_design, abundance_paths)
sample_labels <- study_design$sample_label
population <- study_design$population
abundance_paths <- study_design$abundance_path
sample_dimensions <- c('population', 'age')


tx_to_gene_df <- get_transcript_to_gene_df(EnsDb.Mmusculus.v79)

c(gene_counts, gene_lengths, gene_abunds) %<-% 
	get_gene_level_stats_dfs(abundance_paths, sample_labels, tx_to_gene_df)

#TODO build some sort of relevant plot for these abundances
# gene_abunds <- add_row_descriptive_stats(gene_abunds, sample_labels)

dge_list <- build_digital_gene_expression_list(gene_counts, sample_labels)
# log_cpm_long <- build_log_cpm_df(dge_list, long = TRUE)

dge_list_filt <- filter_dge_list(dge_list,
                                 min_cpm = 2,
                                 min_samples_with_min_cpm = 5)
# log_cpm_filt_long <- build_log_cpm_df(dge_list_filt, long = TRUE)

# adjusts norm factors in dge list samples df, allows comparison across samples
dge_list_filt_norm <- calcNormFactors(dge_list_filt, method = 'TMM')
# log_cpm_filt_norm_long <- build_log_cpm_df(dge_list_filt_norm, long = TRUE)

# plt_1 <- build_log_cpm_plot(log_cpm_long, "unfiltered, non-normalized")
# plt_2 <- build_log_cpm_plot(log_cpm_filt_long, "filtered, non-normalized")
# plt_3 <- build_log_cpm_plot(log_cpm_filt_norm_long, "filtered, normalized")
# plot_grid(plt_1, plt_2, plt_3, labels = c('A', 'B', 'C'), label_size = 12)


log_cpm_filt_norm <- build_log_cpm_df(dge_list_filt_norm, long = FALSE)
write_cluster_dendogram_plot(log_cpm_filt_norm, sample_labels, cluster_out_path)

pca_metrics <- get_pca_metrics(log_cpm_filt_norm)
write_pca_scatter_plots(pca_metrics, sample_dimensions, study_design, pca_scatter_out_path)
write_pca_small_multiples_plots(pca_metrics, sample_dimensions, study_design,
								pca_small_multiples_out_path)

mouse_archs4_rnaseq_path = '/media/amundy/Windows/bio/archs4_rnaseq/mouse_matrix_v10.h5'
# h5ls(mouse_archs4_rnaseq_path)
all_arch_sample_geos <- h5read(mouse_archs4_rnaseq_path, name="meta/samples/geo_accession")
ext_sample_geos <- c("GSM2310941", "GSM2310942", "GSM2310943", "GSM2310944", "GSM2310945", "GSM2310946", "GSM2310947", "GSM2310948", "GSM2310949", "GSM2310950", "GSM2310951", "GSM2310952")
ext_idxs <- which(all_arch_sample_geos %in% ext_sample_geos)

ext_sample_source_name <- 
	h5read(mouse_archs4_rnaseq_path, "meta/samples/source_name_ch1")[ext_idxs]
ext_sample_labels <- 
	h5read(mouse_archs4_rnaseq_path, name="meta/samples/title")[ext_idxs]
#TODO consider trying to automate the production of these, but dont do it prematurely
# These are derived from the title in this case, but it varies where to find them
ext_sample_population <- c("WT", "WT", "Ripk3", "Ripk3", "Ripk3Casp8", "Ripk3Casp8", "WT", "WT", "Ripk3", "Ripk3", "Ripk3Casp8", "Ripk3Casp8")
ext_sample_treatment <- c("unstim", "unstim", "unstim", "unstim", "unstim", "unstim", "LPS", "LPS", "LPS", "LPS", "LPS", "LPS")

ext_study_design <- data.frame(
	population = ext_sample_population,
	treatment = ext_sample_treatment,
	sample_label = ext_sample_labels,
	source_name = ext_sample_source_name
)

gene_ids <- h5read(mouse_archs4_rnaseq_path, "meta/genes/gene_symbol")

# get gene expressions, wide by sample, then transcribe and add col/rownames
expression <- h5read(mouse_archs4_rnaseq_path, "data/expression", 
					 index=list(ext_idxs, 1:length(gene_ids)))
expression <- t(expression)
rownames(expression) <- gene_ids
colnames(expression) <- all_arch_sample_geos[ext_idxs]
ext_dge_list <- DGEList(expression)
ext_dge_list <- filter_dge_list(ext_dge_list,
								min_cpm = 2,
								min_samples_with_min_cpm = 2)

ext_dge_list <- calcNormFactors(ext_dge_list, method = "TMM")
ext_cpm <- build_log_cpm_df(ext_dge_list, long = FALSE)

ext_pca_metrics <- get_pca_metrics(ext_cpm)
write_pca_scatter_plots(ext_pca_metrics, c('population', 'treatment'), 
						ext_study_design, pca_scatter_ext_out_path)
write_pca_small_multiples_plots(ext_pca_metrics, c('population', 'treatment'), 
								ext_study_design, pca_small_multiples_ext_out_path)

age_factor <- factor(study_design$age)
design <- model.matrix(age_factor)
colnames(design) <- levels(age_factor)


