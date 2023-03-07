#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(rmarkdown)
  })

# TODO needed to add this to get it runnable from command line- figure out a
# way around this (maybe with command line args)
setwd("~/code/bio/")

# import configuration details to get output_dir and output_name/title
source('./bulk_rna_seq_analysis/config.R')
title_arg = paste0("--metadata=title:", output_title)

# render the r script as an rmarkdown file (or html file)
output_file = paste0(output_dir, output_name)
rmarkdown::render('./bulk_rna_seq_analysis/analysis.R',
                  output_file=output_file,
                  output_options = list(pandoc_args = c(title_arg)),
                  knit_root_dir = getwd())

browseURL(output_file)

# for storing/recording the libraries used in this project, 
# for docker/reproducibility purposes
# renv::init(project='/home/awmundy/code/bio/r_code/')
