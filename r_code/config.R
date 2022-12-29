
get_inputs_ac_thymus <- function() {
  inputs =   list(
    abundance_root_dir =
      '/media/awmundy/Windows/bio/ac_thymus/rna_txs/fastq_folders/',
    study_design_path =
      '/media/awmundy/Windows/bio/ac_thymus/study_design/study_design_removed_bad_one.csv'
    # isoform_annotation_path =
    # '/media/awmundy/Windows/bio/reference_genomes/mouse/gencode.vM31.chr_patch_hapl_scaff.annotation.gtf.gz',
    # external_validation_fasta_reference_path =
    # '/media/awmundy/Windows/bio/reference_genomes/mouse/Mus_musculus.GRCm39.cdna.all.fa.gz',
    # external_validation_archs4_rnaseq_path =
    #   '/media/awmundy/Windows/bio/archs4_rnaseq/mouse_matrix_v10.h5'
  )
  return(inputs)
}

sample_dimensions <- c('population', 'age')
# explanatory_variable <- c('age')
explanatory_variable <- c('age')
if (explanatory_variable == 'age') {
  output_name <- 'ac_thymus_young_vs_adult_comparison.html'
  control_label <- 'old'
  experimental_label <- 'young'
} else if (explanatory_variable == 'population') {
  control_label <- 'cx3_neg'
  experimental_label <- 'cx3_pos'
  output_name <- 'ac_thymus_population_comparison.html'  
} else{
  stop('Invalid explanatory variable')
}

multiple_testing_correction_method <- "BH"
# Must be FALSE if knitting to Rmarkdown
write_output <- FALSE
compare_to_external_data <- FALSE
c(abundance_root_dir, study_design_path) %<-% get_inputs_ac_thymus()
output_dir <- '/media/awmundy/Windows/bio/ac_thymus/outputs/'



