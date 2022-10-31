###Kallisto Quant
- Description:
  - Calculates transcript abundance, which is the estimated frequency of transcripts in the sample
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


###tximport
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
  - https://f1000research.com/articles/4-1521
    - Motivation for us Using gene level data rather than transcript level data
- Notes
  - Multiple Unit options
    - scaledTPM: scaled up to library size by multiplying TPM counts by the library size in
      millions. These counts sum to the mapped read count (library size)
    - lengthScaledTPM: scaled using the average transcript length over samples and then scaled
      as in scaledTPM. These counts sum to the mapped read count (library size)
    - Can optionally derive counts from the abundances (this is what is done for this analysis)


###DGEList
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


###Filter DGEList
- Description:
  - DGEList gene counts that don't pass a volume threshold within sample and between
    samples are filtered out
- Purpose:
  - To filter out genes that are so lowly expressed across samples that analysis
    on them is not useful
- Record Level:
  - Gene
- Unit:
  - Scaled TPM
- Notes:
  - May need to adjust filter cutoff depending on read depth, to maximize
    number of genes with a reasonable false discovery rate downstream


###Normalize DGEList
- Description:
  - Uses Trimmed Mean of M values (TMM) to build sample level normalization factors and
    store them in the sample level df
- Purpose:
  - TMM controls for differences in library sizes (read depth) and trims highly expressed genes'
    counts to mitigate crowding out on the sequencing machine
  - Allows comparisons across samples
- Record Level (Expression):
  - Gene
- Unit (Expression):
  - Scaled TPM (unchanged from previous step)
- Record Level (Normalization):
  - Sample
- Unit (Normalization):
  - Normalization Factors and Library Sizes


###CPM/Log 2 transform
- Description: Converts counts to Normalzed Counts Per Million (CPM) and Log2 transforms it
  - CPM calc: Counts from DGEList$counts / library size from DGEList Samples * 1,000,000
  - Applies normalization factors too
- Purpose:
  - Logging helps mitigate heteroskedasticity (variance increases as the genes count increases)
    - Heteroskedasticity makes comparing highly expressed genes across samples challenging
    - Log transforming the counts mitigates this, although the resulting dataset is
      slightly heteroskedastic towards lowly expressed genes
  - Log2CPM TODO fill in


###voom
Description:
  - Using a supplied design matrix and DGEList object, variance stabilizes the CPM.
  - Currently requires counts not CPM and performs the log2CPM itself, duplicatively
    with the previous CPM step
Purpose:
- Record Level:
  - Gene
- Unit :
  - Log2CPM, variance stabilized


###lmFit
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


###contrast fit
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


###eBayes
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

###topTable
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

###decideTests
- Description:
  - Similar to topTable, but used to produce datatables (pretty tables) for some reason

###Clustering
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
  

### Gene Set Enrichment analysis
- Documentation:
  - https://sci-hub.se/https://www.nature.com/articles/s41596-018-0103-9


