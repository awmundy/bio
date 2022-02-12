## Conda Guide
https://protocols.hostmicrobe.org/conda

### Conda envs ubuntu
- `rnaseq`
  - kallisto, fastqc, multiqc
- `bio` (python 3.8 due to dependency on kb-python)
  - kb-python 
- `sourmash`
  - sourmash
- `centrifuge`
  - centrifuge

#R
Installed Rstudio
- https://www.rstudio.com/products/rstudio/download/#download

Upgraded R 
- after installing R studio I was only at R 3.6 or something
- copied the commands from here, adding sudo to the front of them
  - https://cran.r-project.org/ (Ubuntu)

Installed Bioconductor for R, following lecture
- pasted the following in to Rstudio console
  - `if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
      BiocManager::install(version = "3.13"`

Style
- Rstudio can style code (default shortcut Ctrl + Alt + A)
- Default styling was weird, installed the `styler` package and it automatically did indentation better