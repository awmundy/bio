library(zeallot)

cfgs <- 
  list('ac_thymus' =
         list(
           abundance_root_dir =
             '/media/awmundy/Windows/bio/ac_thymus/rna_txs/fastq_folders/',
           output_dir = '/media/awmundy/Windows/bio/ac_thymus/outputs/',
           sample_dimensions = c('population', 'age'),
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
           output_dir = '/media/awmundy/Windows/bio/senescence_skeletal_muscle_myofiber/outputs/',
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
                   output_title = 'Young vs Old-veh Comparison Senescent Skeletal Muscle Myofiber, Zhang et al. Replication'
                 ),
               'old_dq' =
                 list(
                   explanatory_variable = 'cohort',
                   study_design_path =
                     '/media/awmundy/Windows/bio/senescence_skeletal_muscle_myofiber/study_design/study_design_old_dq.csv',
                   control_label = 'young',
                   experimental_label = 'old_dq',
                   output_name = 'senescence_skeletal_muscle_myofiber_old_dq_comparison.html',
                   output_title = 'Young vs Old-DQ Comparison Senescent Skeletal Muscle Myofiber, Zhang et al. Replication'
                 )
             )
         ))

# run <- 'ac_thymus'
# run_subtype <- c('age')
# run_subtype <- c('population')

run <- 'senescence_skeletal_muscle_myofiber'
run_subtype <- 'old_dq'

abundance_root_dir <- cfgs[[run]]$abundance_root_dir
study_design_path <- cfgs[[run]]$run_subtype[[run_subtype]]$study_design_path
sample_dimensions <- cfgs[[run]]$sample_dimensions
explanatory_variable <- cfgs[[run]]$run_subtype[[run_subtype]]$explanatory_variable
control_label <- cfgs[[run]]$run_subtype[[run_subtype]]$control_label
experimental_label <- cfgs[[run]]$run_subtype[[run_subtype]]$experimental_label
output_dir <- cfgs[[run]]$output_dir
output_name <- cfgs[[run]]$run_subtype[[run_subtype]]$output_name
output_title <- cfgs[[run]]$run_subtype[[run_subtype]]$output_title


multiple_testing_correction_method <- "BH"
# Must be FALSE if knitting to Rmarkdown
write_output <- FALSE
compare_to_external_data <- FALSE



