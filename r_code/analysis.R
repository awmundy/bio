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
	library(gt)
	library(plotly)
	library(DT)
	library(IsoformSwitchAnalyzeR)
	library(gplots)
	library(RColorBrewer)
}) 

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

write_log_cpm_filter_norm_impact_plots <- function(dge_list, 
												   dge_list_filt, 
												   dge_list_filt_norm,
												   log_cpm_filter_norm_out_path) {
	
	log_cpm_long <- build_log_cpm_df(dge_list, long = TRUE)
	log_cpm_filt_long <- build_log_cpm_df(dge_list_filt, long = TRUE)
	log_cpm_filt_norm_long <- build_log_cpm_df(dge_list_filt_norm, long = TRUE)
	
	plt_1 <- build_log_cpm_plot(log_cpm_long, "unfiltered, non-normalized")
	plt_2 <- build_log_cpm_plot(log_cpm_filt_long, "filtered, non-normalized")
	plt_3 <- build_log_cpm_plot(log_cpm_filt_norm_long, "filtered, normalized")
	
	all_plts <- plot_grid(plt_1, plt_2, plt_3, 
						  labels = c('A', 'B', 'C'), 
						  label_size = 12)
	ggsave(file=log_cpm_filter_norm_out_path, all_plts)
	
}

get_external_sample_data_and_study_design <- function(mouse_archs4_rnaseq_path,
													  ext_sample_metadata) {
	## retrieves sample data from an external source for comparison purposes
	## returns 
	ext_sample_geos <-
		ext_sample_metadata$ext_sample_geos
	ext_sample_population <-
		ext_sample_metadata$ext_sample_population
	ext_sample_treatment <-
		ext_sample_metadata$ext_sample_treatment
	
	# h5ls(mouse_archs4_rnaseq_path)
	
	all_arch_sample_geos <-
		h5read(mouse_archs4_rnaseq_path, name = "meta/samples/geo_accession")
	
	ext_idxs <- which(all_arch_sample_geos %in% ext_sample_geos)
	stopifnot(length(ext_sample_geos) == length(ext_idxs))
	
	ext_sample_source_name <- 
		h5read(mouse_archs4_rnaseq_path, "meta/samples/source_name_ch1")[ext_idxs]
	ext_sample_labels <- 
		h5read(mouse_archs4_rnaseq_path, name="meta/samples/title")[ext_idxs]
	ext_study_design <- data.frame(
		population = ext_sample_population,
		treatment = ext_sample_treatment,
		sample_label = ext_sample_labels,
		source_name = ext_sample_source_name
	)
	
	gene_ids <- h5read(mouse_archs4_rnaseq_path, "meta/genes/gene_symbol")
	
	# get gene expressions, wide by sample, then transcribe and add col/rownames
	# TODO determine if these are in cpm already
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
	
	ext_data <- list(ext_cpm = ext_cpm,
					 ext_study_design = ext_study_design)
	
	return(ext_data)
}


write_external_sample_pca <- function(ext_data, pca_scatter_ext_out_path,
									  pca_small_multiples_ext_out_path) {
	
	ext_pca_metrics <- get_pca_metrics(ext_data$ext_cpm)
	write_pca_scatter_plots(ext_pca_metrics, 
							c('population', 'treatment'), 
							ext_data$ext_study_design, 
							pca_scatter_ext_out_path)
	write_pca_small_multiples_plots(ext_pca_metrics, 
									c('population', 'treatment'), 
									ext_data$ext_study_design, 
									pca_small_multiples_ext_out_path)
}

get_design_matrix <- function(study_design, has_intercept, 
							  explanatory_variable) {
	variable_factor <- factor(study_design[, explanatory_variable])
	if (has_intercept == TRUE) {
		design_matrix <- model.matrix(variable_factor)
	} else {
		design_matrix <- model.matrix(~0 + variable_factor)
	}
	
	# interaction
	# design_matrix <- model.matrix(~age_factor*population_factor)
	# additive
	# design_matrix <- model.matrix(~age_factor + population_factor)
	
	# better column names
	prefix <- paste(explanatory_variable, '_', sep='')
	colnames(design_matrix) <- sub("variable_factor", prefix, colnames(design_matrix))
	
	return(design_matrix)
}

write_dge_volcano_plot <- function(dge_top_volcano, dge_volcano_out_path) {
	
	dge_top_volcano <- as_tibble(dge_top_volcano, rownames='gene_id')
	# gt(deg) # pretty table
	
	# now plot
	vplot <- ggplot(dge_top_volcano) +
		# aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
		aes(y=adj.P.Val, x=logFC, text = paste("Symbol:", gene_id)) +
		geom_point(size=2) +
		#geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
		#geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
		#geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
		#annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
		#annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
		labs(title="Volcano plot",
			 # subtitle = "Insert Subtitle",
			 caption=paste0("produced on ", Sys.time())) +
		theme_bw()
	
	# write interactive plot
	interactive_vplot <- plotly::ggplotly(vplot)
	htmlwidgets::saveWidget(interactive_vplot, dge_volcano_out_path)
}

write_dge_csv_and_datatable <- function(dge, elist, dge_csv_out_path, 
										dge_datatable_out_path) {

	
	write_csv(dge, file=dge_csv_out_path)
	
	# build and write datatable for a pretty output
	dtable <- datatable(dge, 
						extensions = c('KeyTable', "FixedHeader"), 
						caption = 'Differentially Expressed Genes',
						rownames = FALSE,
						options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, 
									   lengthMenu = c("10", "25", "50", "100")))
	round_cols <- names(dtable$x$data)[! names(dtable$x$data) %in% c('gene_id')]
	dtable <- formatRound(dtable, columns=round_cols, digits=2)
	htmlwidgets::saveWidget(dtable, dge_datatable_out_path)
}

temp_isoform_analysis <- function(study_design, explanatory_variable,
								  abundance_paths, isoform_annotation_path,
								  fasta_reference_path,
								  isoform_analysis_out_dir) {
	
	# build design matrix specific to isoform anaysis library
	isoform_design_matrix <- 
		dplyr::select(study_design, sample_label, explanatory_variable)
	isoform_design_matrix <- 
		dplyr::rename(isoform_design_matrix, 
					  sampleID=sample_label, 
					  condition=explanatory_variable)
	
	# build object containing abundances needed for the next step
	abundance_list <- importIsoformExpression(sampleVector = abundance_paths)
	colnames(abundance_list$abundance) <- c("isoform_id", sample_labels) 
	colnames(abundance_list$counts) <- c("isoform_id", sample_labels) 
	
	# TODO not working with the mouse reference data I've tried
	switch_list <- importRdata(
		isoformCountMatrix   = abundance_list$counts,
		isoformRepExpression = abundance_list$abundance,
		designMatrix         = isoform_design_matrix,
		removeNonConvensionalChr = TRUE,
		addAnnotatedORFs=TRUE,
		ignoreAfterPeriod=TRUE,
		isoformExonAnnoation = isoform_annotation_path,
		isoformNtFasta       = fasta_reference_path,
		showProgress = TRUE)
	
	# performs analysis of isoform switches and writes to output location
	switch_list <- isoformSwitchAnalysisCombined(
		switchAnalyzeRlist   = switch_list,
		pathToOutput = isoform_analysis_out_dir)
	
	switch_summary <- extractSwitchSummary(switch_list)
	print(switch_summary)
	
	# get the top switches by usage or q-values
	top_switches <- extractTopSwitches(
		switch_list, 
		filterForConsequences = TRUE, 
		n = 50, 
		sortByQvals = FALSE) 
	
	
	switchPlot(
		switch_list,
		gene='insert_gene_id_here',
		condition1 = 'condition_1_here',
		condition2 = 'condition_2_here',
		localTheme = theme_bw())
	
}

# input paths
abundance_root_dir <- '/media/amundy/Windows/bio/ac_thymus/rna_txs/fastq_folders/'
study_design_path <- '/media/amundy/Windows/bio/ac_thymus/study_design_removed_bad_one.csv'
isoform_annotation_path <- '/media/amundy/Windows/bio/reference_genomes/mouse/gencode.vM31.chr_patch_hapl_scaff.annotation.gtf.gz'
fasta_reference_path <- '/media/amundy/Windows/bio/reference_genomes/mouse/Mus_musculus.GRCm39.cdna.all.fa.gz'
# external validation inputs
mouse_archs4_rnaseq_path = '/media/amundy/Windows/bio/archs4_rnaseq/mouse_matrix_v10.h5'
# TODO replace these with actual genes of interest
ext_sample_metadata <- 
	list(ext_sample_geos = 
		 	c("GSM2310941", "GSM2310942", "GSM2310943", "GSM2310944", "GSM2310945",
		 	  "GSM2310946", "GSM2310947", "GSM2310948", "GSM2310949", "GSM2310950",
		 	  "GSM2310951", "GSM2310952"),
		 ext_sample_population = 
		 	c("WT", "WT", "Ripk3", "Ripk3", "Ripk3Casp8", "Ripk3Casp8", "WT", 
		 	  "WT", "Ripk3", "Ripk3", "Ripk3Casp8", "Ripk3Casp8"),
		 ext_sample_treatment = 
		 	c("unstim", "unstim", "unstim", "unstim", "unstim", "unstim",
		 	  "LPS", "LPS", "LPS", "LPS", "LPS", "LPS"))

# output_paths
log_cpm_filter_norm_out_path <- "/home/amundy/Documents/dge/ac_thymus/filter_norm_impact.pdf"
pca_scatter_out_path <- "/home/amundy/Documents/dge/ac_thymus/pca_scatter.pdf"
pca_small_multiples_out_path <- "/home/amundy/Documents/dge/ac_thymus/pca_small_multiples.pdf"
pca_scatter_ext_out_path <- "/home/amundy/Documents/dge/ac_thymus/pca_scatter_ext.pdf"
pca_small_multiples_ext_out_path <- "/home/amundy/Documents/dge/ac_thymus/pca_small_multiples_ext.pdf"
cluster_out_path <- "/home/amundy/Documents/dge/ac_thymus/cluster.pdf"
dge_volcano_out_path <- "/home/amundy/Documents/dge/ac_thymus/dge_volcano.html"
dge_csv_out_path <- "/home/amundy/Documents/dge/ac_thymus/dge_table.csv"
dge_datatable_out_path <- "/home/amundy/Documents/dge/ac_thymus/dge_table.html"
isoform_analysis_out_dir <- "/home/amundy/Documents/dge/ac_thymus/isoform_analysis/"


study_design <- get_study_design_df(study_design_path)
abundance_paths <- get_abundance_paths(abundance_root_dir)
study_design <- assign_abundance_paths_to_study_design(study_design, abundance_paths)
sample_labels <- study_design$sample_label
population <- study_design$population
abundance_paths <- study_design$abundance_path
sample_dimensions <- c('population', 'age')
explanatory_variable <- c('population')


tx_to_gene_df <- get_transcript_to_gene_df(EnsDb.Mmusculus.v79)

c(gene_counts, gene_lengths, gene_abunds) %<-% 
	get_gene_level_stats_dfs(abundance_paths, sample_labels, tx_to_gene_df)

#TODO build some sort of relevant plot for these abundances
# gene_abunds <- add_row_descriptive_stats(gene_abunds, sample_labels)

dge_list <- build_digital_gene_expression_list(gene_counts, sample_labels)

dge_list_filt <- filter_dge_list(dge_list,
                                 min_cpm = 2,
                                 min_samples_with_min_cpm = 5)

# adjusts norm factors in dge list samples df, allows comparison across samples
dge_list_filt_norm <- calcNormFactors(dge_list_filt, method = 'TMM')
write_log_cpm_filter_norm_impact_plots(dge_list, dge_list_filt, dge_list_filt_norm,
									   log_cpm_filter_norm_out_path)

log_cpm_filt_norm <- build_log_cpm_df(dge_list_filt_norm, long = FALSE)
write_cluster_dendogram_plot(log_cpm_filt_norm, sample_labels, cluster_out_path)

pca_metrics <- get_pca_metrics(log_cpm_filt_norm)
write_pca_scatter_plots(pca_metrics, sample_dimensions, study_design, pca_scatter_out_path)
write_pca_small_multiples_plots(pca_metrics, sample_dimensions, study_design,
								pca_small_multiples_out_path)

## comparing to external sample
ext_data <- get_external_sample_data_and_study_design(mouse_archs4_rnaseq_path,
										  ext_sample_metadata)
write_external_sample_pca(ext_data, pca_scatter_ext_out_path,
						  pca_small_multiples_ext_out_path)

design_matrix <- get_design_matrix(study_design, FALSE, explanatory_variable)

# voom requires the input to be counts, not CPM, TPM etc
# lot2cpm of the counts then variance stabilize them
# $E contains the resulting counts
# object is an Expression List (EList)
elist <- voom(dge_list_filt_norm, design_matrix, plot = TRUE)

# builds a linear model with coefficients for each column in design matrix
# each row is a gene, each value is the average of the transformed 
# counts from voom for that category/gene combo, if explanatory variables 
# are factors
linear_fit <- lmFit(elist, design_matrix)

# makes matrix representing the contrasts that will to be evaluated
# contrast_matrix <- makeContrasts(age_delta=age_old-age_young,
								 # levels=design_matrix)
contrast_matrix <- makeContrasts(cx3_effect=population_cx3_pos-population_cx3_neg,
								 levels=design_matrix)

# coefficients here are the differences between the coefficients of the categories
contrast_fit <- contrasts.fit(linear_fit, contrast_matrix)

# uses empirical bayes method to produce p values showing whether each gene
# has a log fold change greater than 0
bayes_fit <- eBayes(contrast_fit)

# adjust p values and get dataframe of genes sorted by abs log fold change
dge_top_volcano <- topTable(bayes_fit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
write_dge_volcano_plot(dge_top_volcano, dge_volcano_out_path)


# TODO the transformations to the counts here should be enforced to be 
# the same as for the volcano plot
# using the fitted model object (that includes the contrast matrix), 
# get a table-like gene level object that records whether the gene was 
# significantly negative, sig positive, or not sig
all_dge <- decideTests(bayes_fit, method="global", adjust.method="BH",
					   p.value=0.05, lfc=2)
# subset to just the genes that were significantly differently expressed
dge_mtx <- elist$E[all_dge[,1] !=0,]
dge <- as_tibble(dge_mtx, rownames = "gene_id")

write_dge_csv_and_datatable(dge, elist, dge_csv_out_path, 
							dge_datatable_out_path)

# TODO not working yet
# temp_isoform_analysis(study_design, explanatory_variable,
# 					  abundance_paths, isoform_annotation_path,
# 					  fasta_reference_path,
# 					  isoform_analysis_out_dir)

#TODO remove, or at least be more deliberate about how many and 
#which genes to evaluate
dge_mtx = dge_mtx[1:10,]

# Clustering
# TODO consider using the clust command line tool instead
#	- different clustering algorithms produce very different clusters
#	- clust uses an ensembl method to get the consensus clustering across
#	  multiple clustering algorithms
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1536-8

# get correlations at the gene and sample level
gene_cor <- cor(t(dge_mtx), method='pearson')
sample_cor <- cor(dge_mtx, method='spearman')

# compute distance matrixes for each (1-cor is to make the vals be 0 to 2)
# - euclidian distance is the default method
gene_cor_dist <- as.dist(as.dist(1-gene_cor))
sample_cor_dist <- as.dist(1-sample_cor)

# compute heirarchical clusters
gene_clust <- hclust(gene_cor_dist, method='complete')
sample_clust <- hclust(sample_cor_dist, method='complete')

# group the clusters, k is the number of sample categories
gene_clust_groups <- cutree(gene_clust, k=2)

# convert the cluster groups to colors
cluster_colors <- rainbow(length(unique(gene_clust_groups)), start=0.1, end=0.9) 
cluster_colors <- cluster_colors[as.vector(gene_clust_groups)] 

heat_colors <- rev(brewer.pal(name="RdBu", n=11))


# TODO can't easily save this
# dge static heatmap , scale='row' computes z score that 
# scales the expression of the rows (genes) to better highlight 
# between gene differences
heatmap.2(dge_mtx,
		  Rowv = as.dendrogram(gene_clust),
		  Colv = as.dendrogram(sample_clust),
		  RowSideColors = cluster_colors,
		  col = heat_colors, scale = 'row', labRow = NA,
		  density.info = "none", trace = "none",
		  cexRow = 1, cexCol = 1, margins = c(5, 5))
dev.off()

# same as above except no row-wise scaling, making the 
# between sample differences more obvious
heatmap.2(dge_mtx,
		  Rowv = as.dendrogram(gene_clust),
		  Colv = as.dendrogram(sample_clust),
		  RowSideColors = cluster_colors,
		  col = heat_colors, labRow = NA,
		  density.info = "none", trace = "none",
		  cexRow = 1, cexCol = 1, margins = c(5, 5))
dev.off()
