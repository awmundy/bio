### Kallisto Index
- Description:
  - Constructs an kallisto index file based on a reference cDNA genome
- Inputs:
  - Referece genome fasta file
Purpose:
  - Speeds up downstream kallisto abundance calculations

### Kallisto Quant
- Description:
  - Calculates transcript abundance, which is the estimated frequency of transcripts in the sample
- Inputs:
  - Kallisto Index file
  - Bulk RNAseq fastq files
- Purpose:
  - To summarize the transcript distribution
- Record Level:
  - Transcript
- Unit:
  - Transcripts per Kilobase Million (TPM)
    - TPM causes the within sample values to sum to 1,000,000, which is one step towards
      being able to make cross sample comparisons
    - Normalizes for transcript length and read depth
- Documentation
  - https://www.nature.com/articles/nbt.3519
  - https://sci-hub.se/10.1038/nbt.3519


### tximport
- Description:
  - Reads in upstream data and converts to counts if that data is abundances
  - Aggregates abundance transcript level data to gene level
- Purpose:
- Record Level:
  - Gene
- Unit:
  - Scaled TPM
- Documentation
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4712774/
  - https://f1000research.com/articles/7-952/v3
    - Both scaledTPM and lengthScaledTPM are appropriate for differential gene expression
  - https://f1000research.com/articles/4-1521
    - Motivation for us Using gene level data rather than transcript level data
- Notes
  - Multiple Unit options
    - scaledTPM: scaled up to library size by multiplying TPM counts by the library size in
      millions. These counts sum to the mapped read count (library size)
    - lengthScaledTPM: scaled using the average transcript length over samples and then scaled
      as in scaledTPM. These counts sum to the mapped read count (library size). This is what is 
      used in this analysis.
    - Can optionally derive counts from the abundances (this is what is done for this analysis)


### DGEList
- Description:
  - Creates a list for downstream processing that contains
    - The counts from tximport
    - A sample level dataframe with a library size column showing the library size for each sample
- Purpose:
  - Prepare the data for the next step
- Record Level:
  - Gene
- Unit:
  - Scaled TPM


### Filter DGEList
- Description:
  - DGEList gene counts that don't pass a volume threshold within sample and between
    samples are filtered out
- Purpose:
  - To filter out genes that are so lowly expressed across samples that analysis
    on them is not useful. Filtering out lowly expressed genes can also result in a higher
    quantity of differentially expressed genes, because false discovery rate correction then removes 
    fewer genes. The thresholds are set so that at least a few samples have ~10-15 reads per million 
    links: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2906865/ https://f1000research.com/articles/5-1438
- Record Level:
  - Gene
- Unit:
  - Scaled TPM
- Notes:
  - May need to adjust filter cutoff depending on read depth, to maximize
    number of genes with a reasonable false discovery rate downstream


### Normalize DGEList
- Description:
  - Uses Trimmed Mean of M values (TMM) to build sample level normalization factors and
    store them in the sample level df
- Purpose:
  - TMM controls for differences in library sizes (read depth) and trims highly expressed genes' counts to mitigate crowding out caused by very highly expressed genes
  - One step towards allowing comparisons across samples
- Record Level (Expression):
  - Gene
- Unit (Expression):
  - Scaled TPM (unchanged from previous step)
- Record Level (Normalization):
  - Sample
- Unit (Normalization):
  - Normalization Factors and Library Sizes


### CPM/Log 2 transform
- Description: Converts counts to normalzed Counts Per Million (CPM) and Log2 transforms it
  - CPM calc: Normalized counts from DGEList$counts / normalized library size from DGEList Samples * 1,000,000
- Purpose:
  - Logging helps mitigate heteroskedasticity (variance increases as the genes count increases)
    - Heteroskedasticity makes comparing highly expressed genes across samples challenging
    - Log transforming the counts mitigates this, although the resulting dataset is
      slightly heteroskedastic towards lowly expressed genes
  - In conjuction with the TMM values from the previous normalization, CPM allows for comparisons between samples by controlling for different amounts of reads 
    in each sample


### voom
Description:
  - Using a supplied design matrix and DGEList object, variance stabilizes the CPM.
  - Currently requires counts not CPM and performs the log2CPM itself, duplicatively
    with the previous CPM step
  - Voom fits a model to how the average logcpm across samples for each gene predicts 
    the square root of the standard deviation across samples for each gene
  - This mean-variance trend line is then used to produce weights for each gene/sample
    observation
Purpose:
  - Variance of counts may be different across samples, so for cross-sample comparisons 
    this is useful. The weights from voom are used downstream
- Record Level:
  - Gene
- Unit :
  - Log2CPM, with weights in the object to be used for variance stabilization


### lmFit
- Description:
  - Produces coefficients at the gene level using the model specified in the design matrix
    - i.e. if the design matrix has a positive and negative group, there will be a
      column of gene level coefficients for positive and a column for negative
  - Each coefficient is the average expression of that gene for that sample category
    (e.g. positive)
  - Uses the Elist object from voom containing log2 transformed, variance stabilized counts
- Purpose:
  - The resulting average expression for each gene for each category can have statistical
    tests run on them to see if they are significantly different
- Record Level:
  - Gene
- Unit :
  - Coefficient/average expression value for each sample category


### contrast fit
- Description:
  - Constructs gene level coefficient differences between the categories passed in
  - e.g. age_old - age_young = the difference between the coefficients from linear fit
    for these categories
- Purpose:
  - The differences here are produced so that they can be tested to see if they're
    significantly different from 0
- Record Level:
  - Gene
- Unit:
  - Difference in coefficients/average expression values for the given contrast(s)
- Documentation:
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/


### eBayes
Description:
  - Does an empirical Bayes calculation
  - Produces t statistics for if each of the gene level contrasts produced in the
    previous step are equal to 0
    - T statistics are modified somehow (fill in how)
  - Produces F statistics for each gene testing whether all contrasts for that gene
    are equal to 0 (need to figure out how to interpret the F object values, there are
    6 for each gene for a 1 contrast run)
  - Produces P values for whether the log2fold change is greater than 0
Purpose:
  - To produce test statistics and p values for whether genes were differentially expressed
    for the given contrasts
- Record Level:
  - Gene
- Unit (for various gene level outputs):
  - T statistics, F statistics, or P values, among others
- Documentation:
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/

### topTable
- Description:
  - Adjusts P values for false discovery and returns a sorted gene level dataframe for genes
    that meet a maximum P value and minimum log2 fold change
    - P value adjustment method, max P value, min logfc, statistic to sort by, and max number
      of genes to return are all params
    - In this analysis the Benjamini Hoschberg method is used to adjust the P values
- Purpose:
  - To get a sorted table of the genes that are the most highly or most significantly expressed,
    subject to some cutoffs
- Record Level:
  - Gene
- Unit:
  - logfc, average expression, t statistic, adjusted and unadjusted P values

### decideTests
- Description:
  - Similar to topTable, but used to produce datatables (pretty tables) for some reason

### Clustering
- Distances between genes are calculated
- Spearman Correlation
  - Evaluates the correlation between the ranks of the variables
  - Assumes there is a monotonic relationship between the variables (as one variable goes up (down)
    the other only goes up (down) or stays constant)
  - Works with ordinal data (data that has directionality and is discrete, like a pain score
    with values of low, medium, and high), in addition to continuous data
  - Does not assume normal distribution
- Pearson Correlation
  - Evaluates the correlation between the values of the variables
  - The most frequently used correlation measure
  - Does not work with ordinal data (it involves averaging the values, which may not be accurate for
    ordinal data where the distance between some adjacent values (low vs medium, medium vs high) may not
    be the same
  - Assumes normal distribution
  

### Functional Enrichment analysis
- Used to determine if certain sets of genes corresponding to some function are differentially 
  expressed (enriched)
- Documentation:
  - https://sci-hub.se/https://www.nature.com/articles/s41596-018-0103-9
- Externally cosntructed gene sets can be used, these can be curated by experts
- g:Profiler (gost) functional enrichment
  - Uses Gene Ontology mappings from http://geneontology.org/
  - It takes just a list of genes and determines whether those genes represent relevant clusters of genes,
    as opposed to just being a random set of genes
  - In this analysis an R library is used to interface with g:Profiler
- clusterProfiler (GSEA) functional enrichment
  - Uses gene sets passed in by user
- GSEABase is an R libray that helps perform functional enrichment analysis as well
  - Requires a sorted list of genes (by e.g. log fold change)
  - Gene expression is compared across two classes (e.g. experimental and control) and the genes are 
    ranked by this comparison's metric (e.g. logfold change)
  - Outcome measure is whether each gene set is significantly differentially expressed


