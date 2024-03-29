library(zeallot)

cfgs <- 
  list('ac_thymus' =
         list(
           abundance_root_dir =
             '/media/awmundy/Windows/bio/ac_thymus/rna_txs/fastq_folders/',
           output_dir = '/media/awmundy/Windows/bio/ac_thymus/outputs/',
           sample_dimensions = c('population', 'age'),
           custom_gene_sets_path = NULL,
           explanatory_variables =
             list(
               'population' =
                 list(
                   study_design_path = 
                     '/media/awmundy/Windows/bio/ac_thymus/study_design/study_design_removed_bad_one.csv',
                   control_label = 'cx3_neg',
                   experimental_label = 'cx3_pos',
                   output_name = 'ac_thymus_population_comparison.html',
                   output_title = 'Thymus Population Comparison'
                 ),
               'age' =
                 list(
                   study_design_path = 
                     '/media/awmundy/Windows/bio/ac_thymus/study_design/study_design_removed_bad_one.csv',
                   control_label = 'old',
                   experimental_label = 'young',
                   output_name = 'ac_thymus_young_vs_adult_comparison.html',
                   output_title = 'Thymus Young vs Adult Comparison'
                 )
             )
         ),
       'senescence_skeletal_muscle_myofiber' =
         list(
           abundance_root_dir =
             '/media/awmundy/Windows/bio/senescence_skeletal_muscle_myofiber/rna_txs/fastq_folders/',
           sample_dimensions = c('cohort'),
           run_subtype =
             list(
               'old_veh' =
                 list(
                   explanatory_variable = 'cohort',
                   study_design_path =
                     '/media/awmundy/Windows/bio/senescence_skeletal_muscle_myofiber/study_design/study_design_old_veh.csv',
                   control_label = 'young',
                   experimental_label = 'old_veh',
                   output_name = 'senescence_skeletal_muscle_myofiber_old_veh_comparison.html',
                   output_title = 'Young vs Old-veh Comparison Senescent Skeletal Muscle Myofiber, Zhang et al. Replication',
                   output_dir = '/media/awmundy/Windows/bio/senescence_skeletal_muscle_myofiber/outputs/old_veh/'
                 ),
               'old_dq' =
                 list(
                   explanatory_variable = 'cohort',
                   study_design_path =
                     '/media/awmundy/Windows/bio/senescence_skeletal_muscle_myofiber/study_design/study_design_old_dq.csv',
                   control_label = 'young',
                   experimental_label = 'old_dq',
                   output_name = 'senescence_skeletal_muscle_myofiber_old_dq_comparison.html',
                   output_title = 'Young vs Old-DQ Comparison Senescent Skeletal Muscle Myofiber, Zhang et al. Replication',
                   output_dir = '/media/awmundy/Windows/bio/senescence_skeletal_muscle_myofiber/outputs/old_dq/'
                 )
             )
         ),
       'age_related_steatohepatitis' =
         list(
           abundance_root_dir =
             '/media/awmundy/TOSHIBA EXT/age_related_steatohepatitis/rna_txs/fastq_folders/',
           sample_dimensions = c('cohort'),
           custom_gene_sets_path = 'https://raw.githubusercontent.com/awmundy/bio/master/custom_gene_sets/custom_gene_sets_duan_et_al.csv',
           # cpm chosen to be about 10-15 reads per million, based on guidance here: 
           # https://f1000research.com/articles/5-1438
           min_cpm = .5,
           min_samples_with_min_cpm = 2,
           run_subtype =
             list('old' =
                    list(
                      explanatory_variable = 'cohort',
                      study_design_path =
                        '/media/awmundy/Windows/bio/age_related_steatohepatitis/study_design/study_design.csv',
                      control_label = 'young',
                      experimental_label = 'old',
                      output_title = 'Differential Gene Expression Replication of Duan et al Steatohepatitis Aging',
                      output_dir = '/media/awmundy/Windows/bio/age_related_steatohepatitis/outputs/young_vs_old/',
                      #TODO generalize this
                      output_html_file = '~/code/bio/docs/index_temp.html'
                    )
       )))


# run <- 'ac_thymus'
# run_subtype <- c('age')
# run_subtype <- c('population')

run <- 'age_related_steatohepatitis'
run_subtype <- 'old'
multiple_testing_correction_method <- "BH"
# Must be FALSE if knitting to Rmarkdown
write_output <- FALSE
compare_to_external_data <- FALSE
  
abundance_root_dir <- cfgs[[run]]$abundance_root_dir
study_design_path <- cfgs[[run]]$run_subtype[[run_subtype]]$study_design_path
sample_dimensions <- cfgs[[run]]$sample_dimensions
explanatory_variable <- cfgs[[run]]$run_subtype[[run_subtype]]$explanatory_variable
control_label <- cfgs[[run]]$run_subtype[[run_subtype]]$control_label
experimental_label <- cfgs[[run]]$run_subtype[[run_subtype]]$experimental_label
custom_gene_sets_path <- cfgs[[run]]$custom_gene_sets_path
min_cpm <- cfgs[[run]]$min_cpm
min_samples_with_min_cpm <- cfgs[[run]]$min_samples_with_min_cpm
output_dir <- cfgs[[run]]$run_subtype[[run_subtype]]$output_dir
output_title <- cfgs[[run]]$run_subtype[[run_subtype]]$output_title
output_html_file <- cfgs[[run]]$run_subtype[[run_subtype]]$output_html_file

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}






