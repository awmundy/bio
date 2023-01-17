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