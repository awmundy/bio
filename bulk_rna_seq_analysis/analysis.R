#' ---
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
  library(gprofiler2)
  library(msigdbr)
  library(clusterProfiler)
  library(enrichplot)
  library(grid)
  library(fgsea)
  library(ggpubr)
  library(cowplot)
  library(knitr)
}) 

#' # Functions and Output Paths
get_abundance_paths <- function(abundance_root_dir, sample_labels) {
  
  abundance_dirs <- list.dirs(abundance_root_dir)[-1]
  
  abundance_paths <- c()
  for (abundance_dir in abundance_dirs)
    # retrieve only the abundance paths corresponding to the sample labels
    if (basename(abundance_dir) %in% sample_labels) {
      abundance_path <- paste(normalizePath(abundance_dir),
                              '/abundance.tsv',
                              sep = '')
      abundance_paths <- c(abundance_paths, abundance_path)
    }
  
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
  tx <- tx[tx$gene_name != '',]
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
  #         placeholder normalization factor of 1
  
  # subset to sample labels cols only, and convert to matrix
  dge_list <- DGEList(as.matrix(gene_counts[, sample_labels]))
  row.names(dge_list$counts) <- gene_counts$gene
  
  return(dge_list)
}


build_log_cpm_df <- function(dge_list, control_label, long) {
  # Construct a logged counts per million df
  
  # uses normalized library sizes by default, as calculated earlier
  log_cpm <- cpm(dge_list, log=TRUE)
  log_cpm <- as_tibble(log_cpm, rownames = "gene_id")
  log_cpm$gene_id <- gsub('"', "", log_cpm$gene_id)
  # move control cols to the front so that later plots are 
  # more straightforward
  log_cpm <- select(log_cpm, starts_with(control_label), everything())
  log_cpm <- select(log_cpm, 'gene_id', everything())

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

get_gene_counts <- function(abundance_paths, sample_labels, tx_to_gene_df) {
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
  # gene_lengths <- convert_tx_gene_mtx_to_tibble(gene_stats_list$length, 
  #                                              sample_labels)
  # gene_abunds <- convert_tx_gene_mtx_to_tibble(gene_stats_list$abundance, 
  #                                              sample_labels)
  
  return(gene_counts)
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
				labs(title="Sample PCA") +
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
    control_label,
    filtering_and_normalizing_impact_out_path,
    write_output
  ) {
	
	log_cpm_long <- 
	  build_log_cpm_df(dge_list, control_label, long = TRUE)
	log_cpm_filt_long <- 
	  build_log_cpm_df(dge_list_filt, control_label, long = TRUE)
	log_cpm_filt_norm_long <- 
	  build_log_cpm_df(dge_list_filt_norm, control_label, long = TRUE)
	
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

build_datatable <- function(data, caption_text) {
  dtbl <- datatable(
    data,
    extensions = c('KeyTable', "FixedHeader"),
    caption = caption_text,
    rownames = FALSE,
    selection = 'multiple',
    filter = 'top',
    options = list(
      keys = TRUE,
      searchHighlight = TRUE,
      pageLength = 10,
      orderMulti = TRUE,
      scrollX = '400px',
      lengthMenu = c("10", "25", "50", "100")
    ))
  
  return(dtbl)
}

plot_dge_datatable <- function(all_dge,
                               sig_dge_datatable_out_path,
                               write_output) {
  all_dge <- dplyr::select(all_dge, -any_of(c('t', 'P.Value', 'B', 'AveExpr')))
  # add flag for if logfc is statistically significant
  all_dge$significant <- with(all_dge, ifelse(adj.P.Val <= .05, 'True', 'False'))
  all_dge <- relocate(all_dge, significant, .after = adj.P.Val)
  
  dtable <- build_datatable(all_dge, 
                           'Gene Expression (counts are log2cpm)')
  round_cols <-
    names(dtable$x$data)[!names(dtable$x$data) %in% c('gene_id', 'significant')]
  dtable <- formatRound(dtable, columns = round_cols, digits = 2)
  
  if (write_output) {
    htmlwidgets::saveWidget(dtable, sig_dge_datatable_out_path)
  } else {
    dtable
  }
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
    
    # makes matrix representing the contrasts that will be evaluated
    experimental_column <- paste0(explanatory_variable, '_', experimental_label)
    control_column <- paste0(explanatory_variable, '_', control_label)
    contrast_string <- paste0(experimental_column, '-', control_column)
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

prep_dge_subset_for_heatmap <- function(dge_subset, log_cpm_filt_norm) {
  # subset to just sample logfc cols to do cluster correlations
  dge_subset <- dge_subset[, colnames(log_cpm_filt_norm)]
  row.names(dge_subset) <- dge_subset$gene_id
  dge_subset <- subset(dge_subset, select=-c(gene_id))
  dge_subset <- as.matrix(dge_subset)
  
  return(dge_subset)
}

get_dge_subset_by_gene_set_for_heatmap <- function(dge, gene_set, 
                                                   log_cpm_filt_norm) {
  
  dge_subset <- dplyr::filter(dge, gene_id %in% gene_set$gene_symbol)
  dge_subset <- prep_dge_subset_for_heatmap(dge_subset, log_cpm_filt_norm)
  
  return(dge_subset)
}

get_most_expressed_dge_subset_for_heatmap <- function(dge, log_cpm_filt_norm) {
  # subset to genes with highest/lowest lfc if necessary
  if (nrow(dge) > 20) {
    dge_high_lfc <- head(dge[order(-dge$logFC),], 10)
    dge_low_lfc <- head(dge[order(dge$logFC),], 10)
    dge_subset <- rbind(dge_high_lfc, dge_low_lfc)
  } else {
    dge_subset <- dge
  }
  
  dge_subset <- prep_dge_subset_for_heatmap(dge_subset, log_cpm_filt_norm)
  
  return(dge_subset)
}

get_gene_cluster_colors <- function(gene_clust) {
  # group the clusters, k is the number of sample categories
  gene_clust_groups <- cutree(gene_clust, k=2)
  
  # convert the cluster groups to colors
  gene_clust_colors <- rainbow(length(unique(gene_clust_groups)), 
                               start=0.1, end=0.9)
  gene_clust_colors <- gene_clust_colors[as.vector(gene_clust_groups)]
  
  return(gene_clust_colors)
}

plot_gene_cluster_heatmap <- function(dge, log_cpm_filt_norm, gene_set=NULL) {
    
    if (is.null(gene_set)) {
      dge_subset <- 
        get_most_expressed_dge_subset_for_heatmap(dge, log_cpm_filt_norm)
      plot_title <- paste('Gene Cluster Heatmap for up to', 
                          '10 highest/lowest LFC genes', 
                          sep= "\n")
    } else {
      gene_set_name <-unique(gene_set$gs_name)
      stopifnot(length(gene_set_name) == 1)
      dge_subset <- 
        get_dge_subset_by_gene_set_for_heatmap(dge, gene_set, log_cpm_filt_norm)
      plot_title <- paste('Gene Cluster Heatmap for', 
                          gene_set_name, sep= "\n")
    }
    
    gene_clust <- get_clusters('gene', dge_subset)
    sample_clust <- get_clusters('sample', dge_subset)
    gene_clust_colors <- get_gene_cluster_colors(gene_clust)
    heat_colors <- rev(brewer.pal(name="RdBu", n=11))
    
    # dge static heatmap , scale='row' computes z score that
    # scales the expression of the rows (genes) to better highlight
    # between gene differences
    gene_scaling_heatmap <-
      heatmap.2(dge_subset,
                Rowv = as.dendrogram(gene_clust),
                dendrogram='row',
                # Colv = as.dendrogram(sample_clust),
                Colv = FALSE,
                RowSideColors = gene_clust_colors,
                col = heat_colors, scale = 'row',
                srtCol = 45,
                density.info = 'none', trace = 'none',
                cexRow = 1, cexCol = 1, 
                # margins = c(4, 4),
                main=plot_title)
    
    # same as above except no row-wise scaling, making the
    # between sample differences more obvious
    # sample_scaling_heatmap <-
    #   heatmap.2(sig_dge_subset,
    #             Rowv = as.dendrogram(gene_clust),
    #             Colv = as.dendrogram(sample_clust),
    #             RowSideColors = gene_clust_colors,
    #             col = heat_colors, scale = 'none',
    #             srtCol = 45,
    #             density.info = 'none', trace = 'none',
    #             cexRow = 1, cexCol = 1, margins = c(8, 8),
    #             main=paste('Gene Cluster Heatmap (sample scaling)', 
    #                        'for up to 10 highest/lowest LFC genes', 
    #                        sep= "\n"))
    # if (write_output) {
    #   png(gene_cluster_heatmap_sample_scaling_out_path)
    #   dev.off()}
  }

get_dge_list_filt_norm <- function(gene_counts, sample_labels, min_cpm, 
                                   min_samples_with_min_cpm) {
  
  dge_list <-
    build_digital_gene_expression_list(gene_counts, sample_labels)
  
  dge_list_filt <- filter_dge_list(dge_list,
                                   min_cpm = min_cpm,
                                   min_samples_with_min_cpm = min_samples_with_min_cpm)
  
  # adjusts norm factors in dge list samples df, allows comparison across samples
  dge_list_filt_norm <-
    calcNormFactors(dge_list_filt, method = 'TMM')
  
  return(list(dge_list, dge_list_filt, dge_list_filt_norm))
}

get_msig_gene_sets <- function(species) {
  msig_hallmarks <- get_msig_hallmark_labels_of_interest()
  
  # load the gene ontology sets, and filter to the ones that are a part of 
  # the hallmarks we care about
  msig_gene_sets <- msigdbr(species=species, category='H') # H for hallmarks
  msig_gene_sets <- dplyr::filter(msig_gene_sets, gs_name %in% msig_hallmarks)
  msig_gene_sets <- dplyr::select(msig_gene_sets, gs_name, gene_symbol)
  
  return(msig_gene_sets)
}

plot_gost_gene_set_enrichment <- function(sig_dge, organism, out_path, 
                                          write_output, title) {
  # Builds and writes an interactive plot showing functional 
  # enrichment by various categories. A gene set is considered functionally 
  # enriched if it statstically (p-value) significantly enriched compared to 
  # all other genes. Uses Gene Ontology mappings (http://geneontology.org/)
  # sig_dge: dataframe of differentially expressed genes
  
  unique_genes <- unique(sig_dge$gene_id)
  
  # TODO determine what correction method is typically used and use that
  # default args return significantly enriched gene sets at a 0.05 threshold
  gost_res <- gost(unique_genes,
                   organism = organism,
                   correction_method = "fdr", 
                   sources= c('GO:BP', 'GO:MF', 'GO:CC'))
  
  out_plot <- gostplot(gost_res, interactive = T, capped = T)
  out_plot <- layout(out_plot, title=title)
  
  if (write_output) {
    htmlwidgets::saveWidget(out_plot, out_path)  
  } else {
    out_plot
  }
}

get_msig_hallmark_labels_of_interest <- function() {

  msig_hallmarks <- c(
  	'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
  	'HALLMARK_ALLOGRAFT_REJECTION',
  	'HALLMARK_COAGULATION',
  	'HALLMARK_INFLAMMATORY_RESPONSE',
  	'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
  	'HALLMARK_P53_PATHWAY',
  	'HALLMARK_IL6_JAK_STAT3_SIGNALING',
  	'HALLMARK_KRAS_SIGNALING_UP',
  	'HALLMARK_INTERFERON_GAMMA_RESPONSE',
  	'HALLMARK_APOPTOSIS',
  	'HALLMARK_COMPLEMENT',
  	'HALLMARK_KRAS_SIGNALING_DN',
  	'HALLMARK_ANGIOGENESIS',
  	'HALLMARK_XENOBIOTIC_METABOLISM',
  	'HALLMARK_APICAL_SURFACE',
  	'HALLMARK_APICAL_JUNCTION',
  	'HALLMARK_HYPOXIA',
  	'HALLMARK_TGF_BETA_SIGNALING',
  	'HALLMARK_UV_RESPONSE_UP',
  	'HALLMARK_MYOGENESIS',
  	'HALLMARK_E2F_TARGETS',
  	'HALLMARK_G2M_CHECKPOINT',
  	'HALLMARK_MITOTIC_SPINDLE',
  	'HALLMARK_MYC_TARGETS_V1',
  	'HALLMARK_SPERMATOGENESIS',
  	'HALLMARK_ANDROGEN_RESPONSE',
  	'HALLMARK_MTORC1_SIGNALING',
  	'HALLMARK_INFLAMMATORY_RESPONSE',
  	'HALLMARK_P53_PATHWAY'
  )
  return(msig_hallmarks)
}

get_go_gene_sets <- function(ont) {
  # ont is GO category (ontology)
  gene_set_list <- getGO(org='mouse', ont=ont)
  
  # merge human readable gene set name onto gene set df
  gene_sets <- merge(gene_set_list$geneset, gene_set_list$geneset_name,
                     by=1, all.x = TRUE)
  
  # cleanup
  gene_sets <- dplyr::rename(gene_sets, gs_name=Term, gene_symbol=gene)
  gene_sets <- dplyr::select(gene_sets, 'gs_name', 'gene_symbol')
  gene_sets <- as_tibble(gene_sets)
  gene_sets <- drop_na(gene_sets)
  
  return(gene_sets)
}

get_custom_gene_sets <- function(custom_gene_sets_path) {
  gene_sets_custom <- read_csv(custom_gene_sets_path, show_col_types = FALSE)
  
  return(gene_sets_custom)
}

get_gene_sets <- function(custom_gene_sets_path) {
  # get msig gene sets
  gene_sets <- get_msig_gene_sets('Mus musculus')

  # get go gene sets, concat to msig gene sets
  go_bp_gene_sets <- msigdbr(species = "Mus musculus", category = "C5",
                             subcategory = "BP")
  go_bp_gene_sets <- dplyr::select(go_bp_gene_sets, "gs_name", "gene_symbol")
  gene_sets <- rbind(gene_sets, go_bp_gene_sets)

  # if there are custom gene sets, attach them too
  if (!is.null(custom_gene_sets_path)) {
    gene_sets_custom <- get_custom_gene_sets(custom_gene_sets_path)
    gene_sets <- rbind(gene_sets, gene_sets_custom)
  }
  
  return(gene_sets)  
}

filter_gsea_df_to_most_sig_pos_and_neg_enriched_pathways <- function(gsea_df) {
  gsea_df_down <- dplyr::filter(gsea_df, NES<0)
  gsea_df_down <- head(gsea_df_down[order(gsea_df_down$p.adjust),], 10)
  gsea_df_up <- dplyr::filter(gsea_df, NES>=0)
  gsea_df_up <- head(gsea_df_up[order(gsea_df_up$p.adjust),], 10)
  gsea_df_filtered <- rbind(gsea_df_down, gsea_df_up)
  
  return(gsea_df_filtered)
}

get_fgsea_gene_set_input <- function(gene_sets) {
  # convert to named list
  gene_sets_fgsea <- split(gene_sets$gene_symbol, gene_sets$gs_name)
  return(gene_sets_fgsea)
}

get_fgsea_input <- function(all_dge) {
  
  fgsea_input <- dplyr::select(all_dge, 'gene_id', 'logFC')
  # Todo convert this to an assertion that it's unique by gene and
  # add .0001 to logFC that are the same to distinguish them
  fgsea_input <- distinct(fgsea_input, gene_id, .keep_all=TRUE) %>%
    distinct(logFC, .keep_all=TRUE)

  # order descending by the ranking metric although it doesn't seem 
  # to make a big difference for fgsea
  fgsea_input <- fgsea_input[order(fgsea_input$logFC),]
  
  # convert to named list
  fgsea_input <- deframe(fgsea_input)
  
  return(fgsea_input)
}

get_fgsea_response_df <- function(gene_sets_fgsea, fgsea_input) {

  fgsea_res <- fgsea(pathways=gene_sets_fgsea, stats=fgsea_input)
  fgsea_df <- as_tibble(fgsea_res)
  
  return(fgsea_df)  
}

plot_fgsea <- function(gene_sets_fgsea_filtered, fgsea_df, fgsea_input, 
                       grid_title_text) {
  plot_list = list()
  for (gs_name in names(gene_sets_fgsea_filtered)) {
    adj_p_value = fgsea_df[fgsea_df$pathway == gs_name,]$padj
    adj_p_value = format(round(adj_p_value, digits=4), nsmall=4)

    nes <- fgsea_df[fgsea_df$pathway == gs_name,]$NES
    nes <- format(round(nes, digits=3), nsmall=3)
    subtitle_str <- paste0('(FDR PVal: ', adj_p_value, ', NES: ', nes, ')')
    title_str <- paste(gs_name, subtitle_str, sep='\n')
    
    plt <- plotEnrichment(pathway=gene_sets_fgsea_filtered[[gs_name]], 
                          stats=fgsea_input) +
      labs(title=title_str) + 
      theme(plot.title = element_text(hjust = 0.5, size=12),
            # plot.margin = margin(0, 0, 0, 0)
            ) +
      xlab('Rank by LogFC')
    plot_list[[gs_name]] <- plt
  }
  
  grid_title <- text_grob(grid_title_text, size = 15, face = "bold")
  ncol <- min(length(plot_list), 2)
  nrow <- floor(length(plot_list)/2) + 1
  # heights <- rep(unit(7, 'cm'), nrow)
  out <- plot_grid(plotlist=plot_list, nrow=nrow, ncol=ncol)
  title <- ggdraw() + draw_label(grid_title_text, fontface='bold')
  print(plot_grid(title, out, ncol=1, rel_heights=c(0.1, 1)))

  # grid.arrange(grobs=plot_list, 
  #              ncol=ncol, 
  #              nrow=nrow,
  #              top=grid_title,
  #              # heights=list(unit(25, 'cm'), unit(10, 'cm'), unit(10, 'cm'))
  #              heights=heights,
  #              )
}

build_and_plot_fgsea <- function(gene_sets, all_dge, grid_title_text) {
  gene_sets_fgsea <- get_fgsea_gene_set_input(gene_sets)
  fgsea_input <- get_fgsea_input(all_dge)
  fgsea_df <- get_fgsea_response_df(gene_sets_fgsea, fgsea_input)
  fgsea_df <- dplyr::select(fgsea_df, -any_of(c('pval', 'log2err')))
  dtbl <- build_datatable(fgsea_df, paste0(grid_title_text, ' Table'))
  dtbl <- formatRound(dtbl, columns = c('padj', 'ES', 'NES'), digits = 4)

  if (nrow(fgsea_df) > 0) {
    fgsea_df_down <- dplyr::filter(fgsea_df, NES<0, padj<.05)
    fgsea_df_down <- head(fgsea_df_down[order(fgsea_df_down$padj),], 5)
    fgsea_df_up <- dplyr::filter(fgsea_df, NES>=0, padj<.05)
    fgsea_df_up <- head(fgsea_df_up[order(fgsea_df_up$padj),], 5)
    fgsea_df_filtered <- rbind(fgsea_df_down, fgsea_df_up)
    gene_sets_fgsea_filtered <- gene_sets_fgsea[fgsea_df_filtered$pathway]
    
    plot_fgsea(gene_sets_fgsea_filtered, fgsea_df, fgsea_input, grid_title_text)
  }
  dtbl
}

# import configuration information and write to disk for recordkeeping
source('./bulk_rna_seq_analysis/config.R')
# source('/home/awmundy/code/bio/bulk_rna_seq_analysis/config.R')
capture.output(cfgs[[run]], file=paste0(output_dir, "config.csv"))

# Output Paths (only used if not knitting to rmarkdown)
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
all_dge_csv_out_path <- 
  paste0(output_dir, "all_dge_table.csv")
sig_dge_datatable_out_path <- 
  paste0(output_dir, "dge_table.html")
isoform_analysis_out_dir <- 
  paste0(output_dir, "isoform_analysis/")
mean_variance_plot_out_path <- 
  paste0(output_dir, "mean_variance_trend.png")
gene_cluster_heatmap_gene_scaling_out_path <- 
  paste0(output_dir, "gene_cluster_heatmap_gene_scaling.png")
gene_cluster_heatmap_sample_scaling_out_path <- 
  paste0(output_dir, "gene_cluster_heatmap_sample_scaling.png")
gost_plot_up_path <- 
  paste0(output_dir, 'gene_ontology_plot_up.html')
gost_plot_down_path <- 
  paste0(output_dir, 'gene_ontology_plot_down.html')
gsea_table_path <- 
  paste0(output_dir, 'gsea_table.xlsx')
gsea_line_plot_path <- 
  paste0(output_dir, 'gsea_line_plot.pdf')
gsea_bubble_plot_path <- 
  paste0(output_dir, 'gsea_bubble_plot.pdf')

#' # Data Preparation
#' ### First we read in our study design file and construct a design matrix.
#' ### These objects contain information about our samples, 
#' ### our variable of interest, and our model specification.
study_design <- get_study_design_df(study_design_path)
knitr::kable(study_design)
sample_labels <- study_design$sample_label
abundance_paths <- get_abundance_paths(abundance_root_dir, sample_labels)
study_design <- assign_abundance_paths_to_study_design(study_design, abundance_paths)
design_matrix <- get_design_matrix(study_design, FALSE, explanatory_variable)
knitr::kable(design_matrix)
abundance_paths <- study_design$abundance_path


# Read abundances and build digital gene expression lists
tx_to_gene_df <- get_transcript_to_gene_df(EnsDb.Mmusculus.v79)

gene_counts <- get_gene_counts(abundance_paths, sample_labels, tx_to_gene_df)


# Normalization, Filtering, Logging, Converting to Counts per Million
c(dge_list, dge_list_filt, dge_list_filt_norm) %<-% 
  get_dge_list_filt_norm(gene_counts, sample_labels, min_cpm, 
                         min_samples_with_min_cpm)

log_cpm_filt_norm <- build_log_cpm_df(dge_list_filt_norm, control_label, 
                                      long = FALSE)

#' # Flowchart
# error=FALSE is required be there is not actually an error
include_graphics(paste0(getwd(), "/documentation/pipeline_flowchart.png"), 
                 dpi = 110, error = FALSE)


#' # Differential Gene Expression
# Build Mean/Variance Weights Across Samples for Each Gene 
mean_variance_weights <- get_mean_variance_weights(dge_list_filt_norm, 
                                                   design_matrix)

# Build Differential Gene Expression Dataframe
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

# Differential Gene Expression Volcano Plots
plot_dge_volcano(sig_dge,
                 'Significantly Differentially Expressed Genes',
                 dge_volcano_sig_out_path,
                 write_output)

# Differential Gene Expression Table
# merge gene level df with sample cols with gene level df with significance info
sig_dge <- merge(sig_dge, log_cpm_filt_norm, by='gene_id', all.x = TRUE)
all_dge <- merge(all_dge, log_cpm_filt_norm, by='gene_id', all.x = TRUE)
# throw error if any row of sig_dge failed to find a merge
stopif(sum(is.na(sig_dge[colnames(log_cpm_filt_norm[2])])) > 0)

write_csv(sig_dge, file=dge_csv_out_path)
write_csv(all_dge, file=all_dge_csv_out_path)
plot_dge_datatable(all_dge, sig_dge_datatable_out_path, write_output)

#' # Gene Set Enrichment Analysis (GSEA) for Custom Gene Sets
gene_sets_custom <- get_custom_gene_sets(custom_gene_sets_path)
build_and_plot_fgsea(gene_sets_custom, all_dge, 'GSEA for Custom Gene Sets')

#' # GSEA for MSIG Gene Sets
gene_sets_msig <- get_msig_gene_sets('Mus musculus')
build_and_plot_fgsea(gene_sets_msig, all_dge, 
                     'GSEA for Most Significant Highest/Lowest NES MSIGDB Gene Sets')

#' # Gene Cluster Differential Expression Heatmap for Custom Gene Sets
gene_sets_custom <- get_custom_gene_sets(custom_gene_sets_path)
for (gene_set_label in unique(gene_sets_custom$gs_name)) {
  gene_set <- 
    dplyr::filter(gene_sets_custom, gene_sets_custom$gs_name == gene_set_label)
  plot_gene_cluster_heatmap(all_dge, log_cpm_filt_norm, gene_set)
}

#' # Gene Cluster Differential Expression Heatmap for Highest/Lowest LFC Genes
plot_gene_cluster_heatmap(sig_dge, log_cpm_filt_norm)

#' # Gene Ontology GOST plots
# split out into upregulated and downregulated sets
sig_dge_up <- dplyr::filter(sig_dge, logFC >= 0)
sig_dge_down <- dplyr::filter(sig_dge, logFC < 0)
plot_gost_gene_set_enrichment(sig_dge_up, 'mmusculus',
                              gost_plot_up_path, write_output,
                              'Significantly Upregulated Pathways')
Sys.sleep(5)
plot_gost_gene_set_enrichment(sig_dge_down, 'mmusculus',
                              gost_plot_down_path, write_output,
                              'Significantly Downregulated Pathways')

#' # QC
# Principal Component Analysis
pca_metrics <- get_pca_metrics(log_cpm_filt_norm)
plot_pca_scatter(pca_metrics, sample_dimensions, study_design,
                 pca_scatter_out_path, write_output)
plot_pca_small_multiples(pca_metrics, sample_dimensions, study_design,
                         pca_small_multiples_out_path, write_output)
# sleep to stop knitr bug where plots repeat
Sys.sleep(5)
plot_impact_of_filtering_and_normalizing(dge_list, dge_list_filt, 
                                         dge_list_filt_norm,
                                         control_label,
                                         filtering_and_normalizing_impact_out_path, 
                                         write_output)

plot_sample_cluster_dendogram(log_cpm_filt_norm, sample_labels, 
                              sample_cluster_out_path, write_output)
plot_mean_variance_distribution(mean_variance_weights, 
                                mean_variance_plot_out_path,
                                write_output)

