# R
### Installing R
https://cran.rstudio.com/

### Installing rstudio
https://posit.co/download/rstudio-desktop/

### Installing r-base-dev
- In order to install some R libraries, various build tools need to be installed
- `r-base-dev` contains these tools

### Installing bioconductor packages
- Bioconductor is a suite of R packages useful for biological work
- Change version number below as needed, i.e. for compatability with different R versions, and 
  then paste the following into the rstudio console

`if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")`
- may need to change file permissions in package directory on linux
- requires `make`, and `gcc`, `liblapack-dev`, `libopenblas-dev`, `libcurl4-openssl-dev`, `libxml2-dev`, 
  `libfontconfig1-dev`, `r-cran-rjava`, to be installed (some just on linux)
- `cluserProfiler` issue with `patchworkGrob`
  - `patchwork`, a dependency for `clusterprofiler`, is also an unrelated package in another repository
  - need to uninstall the other `patchwork`, `setRepositories()` to just CRAN and Bioc software, and then install
    the desired patchwork to get the correct one installed
  - then `install.packages("clusterProfiler")` should work

# Conda environment
### Installing conda
https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

### Setting conda channels
`conda config --add channels bioconda`
`conda config --add channels conda-forge`

### Creating conda environment
`conda env create -f <repo_dir>/bio/env/bio.yml`

