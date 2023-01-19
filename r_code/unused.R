get_external_sample_data_and_study_design <- function(archs4_rnaseq_path,
                                                      control_label) {
  
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
  ext_cpm <- build_log_cpm_df(ext_dge_list, control_label, long = FALSE)
  
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

isoform_analysis <- function(study_design, explanatory_variable,
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

get_gsea_input <- function(sig_dge) {
  # builds a sorted differentially expressed gene input for gsea analysis
  
  # drop dupes on both gene id and logfc bc gsea function doesn't handle 
  # them consistently
  no_dupes_sig_dge <- distinct(sig_dge, gene_id, .keep_all=TRUE) %>% 
    distinct(logFC, .keep_all=TRUE) 
  
  # build sorted list since the GSEA function needs it that way
  gsea_input <- no_dupes_sig_dge$logFC
  names(gsea_input) <- as.character(no_dupes_sig_dge$gene_id)
  gsea_input <- sort(gsea_input, decreasing = TRUE)
  
  return(gsea_input)
}

get_gsea_res <- function(gsea_input, gene_sets) {
  # competitive GSEA with gene level permutations
  gsea_res <- GSEA(gsea_input, TERM2GENE=gene_sets)
  
  return(gsea_res)  
}

plot_gsea_bubble <- function(gsea_df_filtered, write_output, gsea_bubble_plot_path) {
  # TODO figure out if bubble plot axis labels are correct/phenotype in df is correct
  
  # bubble plot
  # - bubble size: number of genes in the gene set
  # - color: enrichment score
  # - transparency: -log10 adjusted p value
  bubble_plot_df <- mutate(gsea_df_filtered,
                           phenotype = case_when(NES > 0 ~ "intervention",
                                                 NES < 0 ~ "control"))
  plt <- ggplot(bubble_plot_df, aes(x = phenotype, y = ID)) +
    geom_point(aes(
      size = setSize,
      color = NES,
      alpha = -log10(p.adjust)
    )) +
    scale_color_gradient(low = "blue", high = "red") +
    theme_bw()
  
  if (write_output) {
    pdf(gsea_bubble_plot_path)
    plt
    dev.off()
  } else {
    plt
  }
}

plot_gsea_line <- function(gsea_res, gsea_df_filtered, write_output, 
                           gsea_line_plot_path) {
  
  for (idx in 1:nrow(gsea_df_filtered)) {
    print(gseaplot2(gsea_res, geneSetID = gsea_df_filtered$ID[idx],
                    title = gsea_df_filtered$Description[idx]))
  } 
  
  # if (write_output) {
  #   pdf(gsea_line_plot_path)
  #   plt
  #   dev.off()
  # } else {
  #   plt
  # }
}

plot_gsea_datatable <- function(gsea_df, write_output, 
                                gsea_datatable_out_path) {
  dtable <- datatable(gsea_df, 
                      extensions = c('KeyTable', "FixedHeader"), 
                      caption = 'Differentially Expressed Genes',
                      rownames = FALSE,
                      selection = 'multiple',
                      filter = 'top',
                      options = list(keys = TRUE, searchHighlight = TRUE,
                                     pageLength = 10, orderMulti = TRUE,
                                     scrollX='400px',
                                     lengthMenu = c("10", "25", "50", "100")))
  # round_cols <- names(dtable$x$data)[! names(dtable$x$data) %in% c('gene_id')]
  # dtable <- formatRound(dtable, columns=round_cols, digits=2)
  
  if (write_output) {
    htmlwidgets::saveWidget(dtable, gsea_datatable_out_path)
  } else {
    dtable
  }
}
# Alternate GSEA construction using the GSEA library
# gene_sets <- get_gene_sets(custom_gene_sets_path)
# gsea_input <- get_gsea_input(sig_dge)
# gsea_res <- get_gsea_res(gsea_input, gene_sets)
# gsea_df <- as_tibble(gsea_res@result)
# plot_gsea_datatable(gsea_df, write_output, gsea_datatable_out_path)
# if (nrow(gsea_df) > 0) {
#   # subset to a reasonable number of gene sets for plotting
#   gsea_df_filtered <-
#     filter_gsea_df_to_most_sig_pos_and_neg_enriched_pathways(gsea_df)
#   plot_gsea_line(gsea_res, gsea_df_filtered, write_output, gsea_line_plot_path)
#   plot_gsea_bubble(gsea_df_filtered, write_output, gsea_bubble_plot_path)
# }

