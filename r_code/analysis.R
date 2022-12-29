#' ---
#' title: Thymus Analysis 
#' date: "<center> _`r Sys.Date()`_ <center>"
#' output:
#'   html_document:
#'     code_folding: hide
#'     df_print: paged
#'     theme: yeti
#'     highlight: tango
#'     toc: yes
#'     number_sections: true
#' ---
#' ```{r setup, include=FALSE}
#' knitr::opts_chunk$set(warning=FALSE, message=FALSE)
#' ```

#' # Libraries
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
  library(rmarkdown)
  library(ggplot2)
  library(renv)
}) 

#' # Functions
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
  row.names(dge_list$counts) <- gene_counts$gene
  
  return(dge_list)
}


build_log_cpm_df <- function(dge_list, long) {
  # Construct a logged counts per million df
  
  log_cpm <- cpm(dge_list, log=TRUE)
  log_cpm <- as_tibble(log_cpm, rownames = "gene_id")
  log_cpm$gene_id <- gsub('"', "", log_cpm$gene_id)

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

plot_sample_cluster_dendogram <- function(tbl, sample_labels, sample_cluster_out_path,
                                          write_output) {
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
	if (write_output) {
	  pdf(sample_cluster_out_path)
	  plot(clusters, labels=sample_labels)
	  dev.off()
	} else {
	  plot(clusters, labels=sample_labels)
	}	

}

plot_pca_scatter <- function(pca_metrics, sample_dimensions,
									study_design, pca_out_path, write_output) {
	sample_labels <- study_design$sample_label

		plot_list = list()
	for(sample_dimension in sample_dimensions) {
		# build factor for coloring
		sample_dimension_factor <- factor(study_design[, sample_dimension])
		# length 1 means factor failed, e.g. doesn't work on tibbles
		stopif(length(sample_dimension_factor) == 1)
		
		# use just odd elements of range so that PCs aren't repeated across plots
		for(i in seq(1, 4, 2)) {
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
	
	if (write_output==TRUE) {
	  pdf(pca_out_path, onefile=TRUE, width=4, height=4)
	  for (plt in plot_list) {
	    replayPlot(plt)
	  }
	  graphics.off()
	}
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

plot_pca_small_multiples <- function(pca_metrics,
                                     sample_dimensions,
                                     study_design,
                                     pca_small_multiples_out_path,
                                     write_output) {
	
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
	
	if (write_output==TRUE) {
	  pdf(pca_small_multiples_out_path, onefile = TRUE)
	  for (plt in plot_list) {
	    replayPlot(plt)
	  }
	  graphics.off()
	}
}

describe <- function(df_col) {
	min_val <- min(df_col, na.rm = T)
	max_val <- max(df_col, na.rm = T)
	mean_val <- mean(df_col, na.rm = T)
	q5 <-  quantile(df_col, 0.05, na.rm = T)
	q10 <- quantile(df_col, 0.10, na.rm = T)
	q25 <- quantile(df_col, 0.25, na.rm = T)
	q50 <- quantile(df_col, 0.50, na.rm = T)
	q75 <- quantile(df_col, 0.75, na.rm = T)
	q90 <- quantile(df_col, 0.90, na.rm = T)
	q95 <- quantile(df_col, 0.95, na.rm = T)
	out <- data.frame(min=min_val, max=max_val, mean=mean_val, 
					  q5=q5, q10=q10, q25=q25, q50=q50, q75=q75, 
					  q90=q90, q95=q95, row.names = NULL)
	out <- t(out)
	return(out)
}

plot_impact_of_filtering_and_normalizing <- 
  function(
    dge_list,
    dge_list_filt,
    dge_list_filt_norm,
    filtering_and_normalizing_impact_out_path,
    write_output
  ) {
	
	log_cpm_long <- build_log_cpm_df(dge_list, long = TRUE)
	log_cpm_filt_long <- build_log_cpm_df(dge_list_filt, long = TRUE)
	log_cpm_filt_norm_long <- build_log_cpm_df(dge_list_filt_norm, long = TRUE)
	
	plt_1 <- build_log_cpm_plot(log_cpm_long, "unfiltered, non-normalized")
	plt_2 <- build_log_cpm_plot(log_cpm_filt_long, "filtered, non-normalized")
	plt_3 <- build_log_cpm_plot(log_cpm_filt_norm_long, "filtered, normalized")
	
	all_plts <- plot_grid(plt_1, plt_2, plt_3, 
						  labels = c('A', 'B', 'C'), 
						  label_size = 12)
	if (write_output){
	  ggsave(file=log_cpm_filter_norm_out_path, all_plts)  
	} else {
	  print(all_plts)
	}
}

get_external_sample_data_and_study_design <- function(archs4_rnaseq_path) {
	
  # TODO if function in use regularly, expose this as an input
  # three vectors of geo ids, cell population labels, and treatment labels,
  # all of same length
  ext_sample_metadata <- 
    list(ext_sample_geos = 
           c(),
         ext_sample_population = 
           c(),
         ext_sample_treatment = 
           c())
  
  ## retrieves sample data from an external source for comparison purposes
	## returns 
	ext_sample_geos <-
		ext_sample_metadata$ext_sample_geos
	ext_sample_population <-
		ext_sample_metadata$ext_sample_population
	ext_sample_treatment <-
		ext_sample_metadata$ext_sample_treatment
	
	# h5ls(archs4_rnaseq_path)
	
	all_arch_sample_geos <-
		h5read(archs4_rnaseq_path, name = "meta/samples/geo_accession")
	
	ext_idxs <- which(all_arch_sample_geos %in% ext_sample_geos)
	stopifnot(length(ext_sample_geos) == length(ext_idxs))
	
	ext_sample_source_name <- 
		h5read(archs4_rnaseq_path, "meta/samples/source_name_ch1")[ext_idxs]
	ext_sample_labels <- 
		h5read(archs4_rnaseq_path, name="meta/samples/title")[ext_idxs]
	ext_study_design <- data.frame(
		population = ext_sample_population,
		treatment = ext_sample_treatment,
		sample_label = ext_sample_labels,
		source_name = ext_sample_source_name
	)
	
	gene_ids <- h5read(archs4_rnaseq_path, "meta/genes/gene_symbol")
	
	# get gene expressions, wide by sample, then transcribe and add col/rownames
	# TODO determine if these are in cpm already
	expression <- h5read(archs4_rnaseq_path, "data/expression", 
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
	
	external_data <- list(ext_cpm = ext_cpm,
					 ext_study_design = ext_study_design)
	
	return(external_data)
}


plot_external_sample_pca <- function(external_data, pca_scatter_ext_out_path,
                                     pca_small_multiples_ext_out_path,
                                     write_output) {
	
	ext_pca_metrics <- get_pca_metrics(external_data$ext_cpm)
	plot_pca_scatter(ext_pca_metrics,
	                 c('population', 'treatment'),
	                 external_data$ext_study_design,
	                 pca_scatter_ext_out_path,
	                 write_output)
	plot_pca_small_multiples(ext_pca_metrics,
	                        c('population', 'treatment'),
	                        external_data$ext_study_design,
	                        pca_small_multiples_ext_out_path,
	                        write_output)
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

get_topTable_dge <- function(bayes_stats,
                             multiple_testing_correction_method,
                             min_lfc,
                             max_p_val) {
  # adjust p values and get dataframe of genes sorted by abs log fold change,
  sig_dge <- topTable(
    bayes_stats,
    adjust = multiple_testing_correction_method,
    coef = 1,
    number = 999999,
    lfc = min_lfc
  )
  
  sig_dge <- as_tibble(sig_dge, rownames='gene_id')
  sig_dge <- dplyr::filter(sig_dge, adj.P.Val <= max_p_val)
  
  return(sig_dge)
}

plot_dge_volcano <- function(dge,
                             title,
                             dge_volcano_out_path,
                             write_output) {
  
  vplot <- ggplot(dge) +
    aes(y = -log10(adj.P.Val),
        x = logFC,
        text = paste("Symbol:", gene_id)) +
    geom_point(size = 2) +
    labs(title = title,
         # subtitle = "Insert Subtitle",
         caption = paste0("produced on ", Sys.time())) +
    theme_bw()
  
  # write interactive plot
  vplot <- plotly::ggplotly(vplot)
  
	if (write_output) {
	  htmlwidgets::saveWidget(vplot, dge_volcano_out_path)
	  } else {
	    vplot
	  }
}

plot_dge_datatable_and_write_csv <- function(sig_dge,
                                             dge_csv_out_path,
                                             dge_datatable_out_path, 
                                             write_output) {
  
  sig_dge <- select(sig_dge, -any_of(c('t', 'P.Value', 'B')))
  
	write_csv(sig_dge, file=dge_csv_out_path)
	
	# build and write datatable for a pretty output
	dtable <- datatable(sig_dge, 
						extensions = c('KeyTable', "FixedHeader"), 
						caption = 'Differentially Expressed Genes',
						rownames = FALSE,
						selection = 'multiple',
						filter = 'top',
						options = list(keys = TRUE, searchHighlight = TRUE,
						               pageLength = 10, orderMulti = TRUE,
						               scrollX='400px',
						               lengthMenu = c("10", "25", "50", "100")))
	round_cols <- names(dtable$x$data)[! names(dtable$x$data) %in% c('gene_id')]
	dtable <- formatRound(dtable, columns=round_cols, digits=2)
	
	if (write_output) {
	  htmlwidgets::saveWidget(dtable, dge_datatable_out_path)
	} else {
	  dtable
	}
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

get_mean_variance_weights <- function(dge_list_filt_norm,
                                      design_matrix) {
  # voom requires the input to be counts, not CPM, TPM etc
  # performs log2cpm of the counts then variance stabilizes them,
  # then estimates the mean-variance relationship and produces observation
  # level weights for linear modelling
  # object is an Expression List (EList)
  
  mean_variance_weights <- voom(dge_list_filt_norm, design_matrix,
                                # plot=TRUE,
                                save.plot = TRUE)
  return(mean_variance_weights)
}

plot_mean_variance_distribution <- function(mean_variance_weights,
                                            mean_variance_plot_out_path,
                                            write_output){
  # points
  df = data.frame(x = mean_variance_weights$voom.xy$x,
                  y = mean_variance_weights$voom.xy$y)
  
  # trend line
  df.line = data.frame(x = mean_variance_weights$voom.line$x,
                       y = mean_variance_weights$voom.line$y)
  
  plt <- 
    ggplot(df, aes(x,y)) + 
    geom_point(size=.1) + 
    theme_bw(15) + 
    geom_line(data=df.line, aes(x,y), color="red") +
    ylim(0, max(df$y))
  interactive_plot <- plotly::ggplotly(plt)
  # add titles/labels here due to bug in ggplotly
  interactive_plot <- plotly::layout(interactive_plot, 
                                     title='Gene Mean vs Variance Across Samples',
                                     xaxis=list(title='log2 count + 0.5'),
                                     yaxis=list(title="Sqrt(standard deviation)"),
                                     margin=list(l=50, r=50, b=100, t=100))
  
  if (write_output) {
    htmlwidgets::saveWidget(interactive_plot, mean_variance_plot_out_path)
  } else {
    interactive_plot
  }
}

get_empirical_bayes_differential_expression_stats <- 
  function(mean_variance_weights, explanatory_variable,
           experimental_label, control_label) {
    
    # builds a linear model with coefficients for each column in design matrix
    # each row is a gene, each value is the average of the transformed 
    # counts from voom for that category/gene combo, if explanatory variables 
    # are factors
    linear_stats <- lmFit(mean_variance_weights, design_matrix)
    
    # makes matrix representing the contrasts that will to be evaluated
    experimental_column <- paste(explanatory_variable, '_', experimental_label, 
                                 sep='')
    control_column <- paste(explanatory_variable, '_', control_label, sep='')
    contrast_string <- paste(experimental_column, '-', control_column, sep='')
    contrast_matrix <- makeContrasts(contrasts=contrast_string,
                                     levels=design_matrix)
    
    # coefficients here are the differences between the coefficients of each 
    # category in the contrast
    contrast_stats <- contrasts.fit(linear_stats, contrast_matrix)
    
    # uses empirical bayes method to produce p values and other statistics 
    # showing whether each gene has a log fold change greater than 0
    bayes_stats <- eBayes(contrast_stats)
    
    return(bayes_stats)
  }

get_sig_dif_expressed_genes <- function(bayes_stats, 
                                        multiple_testing_correction_method,
                                        significance_function,
                                        min_lfc,
                                        max_p_val
) {
  
  
  if (significance_function == 'decideTests') {
    stop('decideTests deg creation not fully implemented yet')
    #TODO dge output needs logFC attached to it, which the decideTests output
    # does not have
    
    # using bayes_stats- the fitted model object (that includes the
    # contrast matrix)- get a gene-level table showing whether the gene was
    # significantly negative, sig positive, or not sig with respect to each
    # coefficient in the contrast matrix
    dge <- decideTests(
      bayes_stats,
      method = "global",
      adjust.method = multiple_testing_correction_method,
      p.value = max_p_val,
      lfc = min_lfc
    )
    dge <- data.frame(dge)
    dge <- as_tibble(dge, rownames = 'gene_id')
    dge$gene_id <- gsub('"', "", dge$gene_id)
    
  } else if (significance_function == 'topTable') {
    dge <- get_topTable_dge(bayes_stats,
                            multiple_testing_correction_method,
                            min_lfc,
                            max_p_val)
  }
  else {
    stop('Invalid significance function')
  }
  
  return(dge)
}

get_clusters <- function(cluster_type, sig_dge_mtx) {
  
  # TODO consider using the clust command line tool instead
  #	- different clustering algorithms produce very different clusters
  #	- clust uses an ensembl method to get the consensus clustering across
  #	  multiple clustering algorithms
  # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1536-8
  
  stopifnot(cluster_type %in% c('gene', 'sample'))
  
  if (cluster_type == 'gene') {
    # use pearson bc gene expression is continuous
    cor_mtx <- cor(t(sig_dge_mtx), method='pearson')
    #TODO comment on why as.dist needs to be done twice
    # (1-cor is to make the vals be 0 to 2)
    cor_dist <- as.dist(as.dist(1-cor_mtx))
  } else if (cluster_type == 'sample') {
    # spearman bc samples discrete and should be ranked
    cor_mtx <- cor(sig_dge_mtx, method='spearman')
    cor_dist <- (as.dist(1-cor_mtx))
  }
  
  clust <-hclust(cor_dist, method='complete')
}

plot_gene_cluster_heatmaps <- 
  function(sig_dge,
           log_cpm_filt_norm,
           gene_cluster_heatmap_gene_scaling_out_path,
           gene_cluster_heatmap_sample_scaling_out_path,
           write_output) {
    
    # subset to genes with highest/lowest lfc if necessary
    if (nrow(sig_dge) > 20) {
      sig_dge_high_lfc <- head(sig_dge[order(-sig_dge$logFC),], 10)
      sig_dge_low_lfc <- head(sig_dge[order(sig_dge$logFC),], 10)
      sig_dge_subset <- rbind(sig_dge_high_lfc, sig_dge_low_lfc)
    } else {
      sig_dge_subset <- sig_dge
    }
    # subset to just sample logfc cols to do cluster correlations
    sig_dge_subset <- sig_dge_subset[, colnames(log_cpm_filt_norm)]
    row.names(sig_dge_subset) <- sig_dge_subset$gene_id
    sig_dge_subset <- subset(sig_dge_subset, select=-c(gene_id))
    sig_dge_subset <- as.matrix(sig_dge_subset)
    
    gene_clust <- get_clusters('gene', sig_dge_subset)
    sample_clust <- get_clusters('sample', sig_dge_subset)
    
    # group the clusters, k is the number of sample categories
    gene_clust_groups <- cutree(gene_clust, k=2)
    
    # convert the cluster groups to colors
    gene_clust_colors <- rainbow(length(unique(gene_clust_groups)), 
                                 start=0.1, end=0.9)
    gene_clust_colors <- gene_clust_colors[as.vector(gene_clust_groups)]
    
    heat_colors <- rev(brewer.pal(name="RdBu", n=11))
    
    # dge static heatmap , scale='row' computes z score that
    # scales the expression of the rows (genes) to better highlight
    # between gene differences
    gene_scaling_heatmap <-
      heatmap.2(sig_dge_subset,
                Rowv = as.dendrogram(gene_clust),
                Colv = as.dendrogram(sample_clust),
                RowSideColors = gene_clust_colors,
                col = heat_colors, scale = 'row',
                srtCol = 45,
                density.info = 'none', trace = 'none',
                cexRow = 1, cexCol = 1, margins = c(8, 8),
                main=paste('Gene Cluster Heatmap (gene scaling)', 
                           'for up to 10 highest/lowest LFC genes', 
                           sep= "\n"))
    if (write_output) {
      png(gene_cluster_heatmap_gene_scaling_out_path)
      dev.off()}
    
    # same as above except no row-wise scaling, making the
    # between sample differences more obvious
    sample_scaling_heatmap <-
      heatmap.2(sig_dge_subset,
                Rowv = as.dendrogram(gene_clust),
                Colv = as.dendrogram(sample_clust),
                RowSideColors = gene_clust_colors,
                col = heat_colors, scale = 'none',
                srtCol = 45,
                density.info = 'none', trace = 'none',
                cexRow = 1, cexCol = 1, margins = c(8, 8),
                main=paste('Gene Cluster Heatmap (sample scaling)', 
                           'for up to 10 highest/lowest LFC genes', 
                           sep= "\n"))
    if (write_output) {
      png(gene_cluster_heatmap_sample_scaling_out_path)
      dev.off()}
  }

get_dge_list_filt_norm <- function(gene_counts, sample_labels) {
  
  dge_list <-
    build_digital_gene_expression_list(gene_counts, sample_labels)
  
  dge_list_filt <- filter_dge_list(dge_list,
                                   min_cpm = 2,
                                   min_samples_with_min_cpm = 5)
  
  # adjusts norm factors in dge list samples df, allows comparison across samples
  dge_list_filt_norm <-
    calcNormFactors(dge_list_filt, method = 'TMM')
  
  return(list(dge_list, dge_list_filt, dge_list_filt_norm))
}

# import configuration information
source('~/code/bio/r_code/config.R')

#' # Output Paths
filtering_and_normalizing_impact_out_path <- 
  paste0(output_dir, "filter_norm_impact.pdf")
pca_scatter_out_path <- 
  paste0(output_dir, "pca_scatter.pdf")
pca_small_multiples_out_path <- 
  paste0(output_dir, "pca_small_multiples.pdf")
pca_scatter_ext_out_path <- 
  paste0(output_dir, "pca_scatter_ext.pdf")
pca_small_multiples_ext_out_path <- 
  paste0(output_dir, "pca_small_multiples_ext.pdf")
sample_cluster_out_path <- 
  paste0(output_dir, "sample_cluster.pdf")
dge_volcano_out_path <- 
  paste0(output_dir, "dge_volcano.html")
dge_volcano_sig_out_path <- 
  paste0(output_dir, "dge_volcano_sig.html")
dge_csv_out_path <- 
  paste0(output_dir, "dge_table.csv")
dge_datatable_out_path <- 
  paste0(output_dir, "dge_table.html")
isoform_analysis_out_dir <- 
  paste0(output_dir, "isoform_analysis/")
mean_variance_plot_out_path <- 
  paste0(output_dir, "mean_variance_trend.png")
gene_cluster_heatmap_gene_scaling_out_path <- 
  paste0(output_dir, "gene_cluster_heatmap_gene_scaling.png")
gene_cluster_heatmap_sample_scaling_out_path <- 
  paste0(output_dir, "gene_cluster_heatmap_sample_scaling.png")

study_design <- get_study_design_df(study_design_path)
abundance_paths <- get_abundance_paths(abundance_root_dir)
study_design <- assign_abundance_paths_to_study_design(study_design, abundance_paths)
design_matrix <- get_design_matrix(study_design, FALSE, explanatory_variable)
sample_labels <- study_design$sample_label
abundance_paths <- study_design$abundance_path


#' # Read abundances and build digital gene expression lists
tx_to_gene_df <- get_transcript_to_gene_df(EnsDb.Mmusculus.v79)

c(gene_counts, gene_lengths, gene_abunds) %<-%
	get_gene_level_stats_dfs(abundance_paths, sample_labels, tx_to_gene_df)


#' # Normalization, Filtering, Logging, Converting to Counts per Million
c(dge_list, dge_list_filt, dge_list_filt_norm) %<-% 
  get_dge_list_filt_norm(gene_counts, sample_labels)

log_cpm_filt_norm <- build_log_cpm_df(dge_list_filt_norm, long = FALSE)

#' # Build Mean/Variance Weights Across Samples for Each Gene 
mean_variance_weights <- get_mean_variance_weights(dge_list_filt_norm, 
                                                   design_matrix)

#' # Principal Component Analysis
pca_metrics <- get_pca_metrics(log_cpm_filt_norm)
plot_pca_scatter(pca_metrics, sample_dimensions, study_design,
                 pca_scatter_out_path, write_output)
plot_pca_small_multiples(pca_metrics, sample_dimensions, study_design,
                         pca_small_multiples_out_path, write_output)

#' # Build Differential Gene Expression Dataframe
bayes_stats <-
  get_empirical_bayes_differential_expression_stats(mean_variance_weights,
                                                    explanatory_variable,
                                                    experimental_label,
                                                    control_label)

sig_dge <- get_sig_dif_expressed_genes(bayes_stats,
                                       multiple_testing_correction_method,
                                       'topTable',
                                       min_lfc = 0,
                                       max_p_val = 0.05)
all_dge <- get_sig_dif_expressed_genes(bayes_stats,
                                         multiple_testing_correction_method,
                                         'topTable',
                                         min_lfc = 0,
                                         max_p_val = 1.0)

#' # Differential Gene Expression Volcano Plots
plot_dge_volcano(sig_dge, 
                 'Significantly Differentially Expressed Genes',
                 dge_volcano_sig_out_path, 
                 write_output)
plot_dge_volcano(all_dge, 
                 'All Genes',
                 dge_volcano_out_path, 
                 write_output)

#' # Differential Gene Expression Table
# merge gene level df with sample cols with gene level df with significance info
sig_dge <- merge(sig_dge, log_cpm_filt_norm, by='gene_id', all.x = TRUE)
# throw error if any row of sig_dge failed to find a merge
stopif(sum(is.na(sig_dge[colnames(log_cpm_filt_norm[2])])) > 0)

plot_dge_datatable_and_write_csv(sig_dge, dge_csv_out_path,
                                 dge_datatable_out_path, write_output)

#' # Gene/Sample Cluster Heatmaps
plot_gene_cluster_heatmaps(sig_dge, 
                           log_cpm_filt_norm,
                           gene_cluster_heatmap_gene_scaling_out_path,
                           gene_cluster_heatmap_sample_scaling_out_path,
                           write_output)

#' # QC
plot_impact_of_filtering_and_normalizing(dge_list, dge_list_filt, 
                                         dge_list_filt_norm,
                                         filtering_and_normalizing_impact_out_path, 
                                         write_output)

plot_sample_cluster_dendogram(log_cpm_filt_norm, sample_labels, 
                              sample_cluster_out_path, write_output)

plot_mean_variance_distribution(mean_variance_weights, 
                                mean_variance_plot_out_path,
                                write_output)

if (compare_to_external_data) {
  external_data <-
    get_external_sample_data_and_study_design(external_validation_archs4_rnaseq_path)
  
  plot_external_sample_pca(
    external_data,
    pca_scatter_ext_out_path,
    pca_small_multiples_ext_out_path,
    write_output
  )
}

# # TODO not working yet
# # temp_isoform_analysis(study_design, explanatory_variable,
# # 					  abundance_paths, isoform_annotation_path,
# # 					  external_validation_fasta_reference_path,
# # 					  isoform_analysis_out_dir)
# 
