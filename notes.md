###Conda envs WSL
- `rnaseq`
  - Kallisto, fastqc, multiqc
- `centrifuge`
  - Centrifuge
- TODO: install kb-python somewhere on WSL if working in wsl for that

### Conda envs Windows
- `app`
  - kb-python
- `rnaseq_win`
  - multiqc
- `sourmash`
  - sourmash

###Tools
- `fastqc`
  - WSL
    - installed in to conda env
  - Windows 
    - can run the GUI version by using the `run_fastqc.bat` file in the fastqc installation dir
    - needed java (https://www.java.com/en/download/)
    - installation guide: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt
    - could not figure out how to get it to run from the command line, 
      the installation guide is inaccurate
    - could not get it in to conda env
- `kallisto`
  - WSL
    - installed in to conda env
  - Windows
    - Had to install outside of conda, by following these instructions: https://chmi-sops.github.io/mydoc_kallisto.html
        - Download files
        - Extract
        - Add C:/Program Files/kallisto to Path variable
        - Did something with command prompt or something to get it to recognize kallisto? i dont remember
- `multiqc` 
  - WSL
    - installed in to conda env
  - Windows
    - installed in to conda env
- `centrifuge`
  - WSL
    - installed in to conda env
  - Windows
    - could not figure out how to get it installed
    - conda doesn't work
    - tried to follow guide here: https://ccb.jhu.edu/software/centrifuge/manual.shtml#obtaining-centrifuge
      - MinGW link was broken, downloaded it from here:https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win32/Personal%20Builds/mingw-builds/installer/mingw-w64-install.exe/download
      - also installed cigwin
      - uncertain how to use those tools to build the binary
- `sourmash`
  - WSL
    - not installed (haven't tried yet)
  - Windows
    - installed in to conda env
- `fastp`
  - WSL
    - not installed (already hav fastqc there)
  - Windows
    - cant install through conda



