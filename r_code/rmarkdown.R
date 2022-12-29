suppressPackageStartupMessages(
  library(rmarkdown)
  )

# import configuration details tp get output_dir and output_name
source('~/code/bio/r_code/config.R')

# render the r script as an rmarkdown file (or html file)
output_file = paste0(output_dir, output_name)
rmarkdown::render('/home/awmundy/code/bio/r_code/analysis.R',
                  output_file=output_file)

browseURL(output_file)

# for storing/recording the libraries used in this project, 
# for docker/reproducibility purposes
# renv::init(project='/home/awmundy/code/bio/r_code/')
