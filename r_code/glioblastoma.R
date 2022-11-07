suppressPackageStartupMessages({
	library(tidyverse)
	library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
	library(Biobase) #base functions for bioconductor; required by GSEABase
	library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
	library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
	library(clusterProfiler) # GSEA tool/plot builder
	library(msigdbr) # download msigdb gene lists
	library(enrichplot) # enrichment plots
	# library(openxlsx)
}) 

# Gene Ontology Analysis (GO Analysis)
# - For a given list of genes (i.e. highly differentially expressed ones), 
#	determine whether certain categorizations (ontologies) of those genes are
#	the ones that are differentially expressed
# - The ontologies are hand curated by experts
# - Ontologies are specific to aspects, listed here:
#	cellular component (where in the cell the expressed products are), 
#	biological process (what function do the expressed product perform), and
#	molecular function (how does the expressed product carry out the function)

# Gene Set Enrichment Analysis (GSEA)
# - Gene expression is compared across two classes (e.g. experimental 
#	and control)
# - The difference in expression is calculated, and genes are then 
#	ranked and ordered by this difference
# - Large differences in the positive and negative direction are at the end,
#	genes that aren't differentially expressed are in the middle
# - Sets of genes (signatures) are evaluated to see if they are positively, 
#	or negatively expressed
# - Running sums are calculated such that if a gene in a gene set is positively
#	expressed, the sum gains +1 and if another gene is negatively expressed the
#	sum gains -1
# - Generally the analysis is for whether the selected gene set as a whole is 
#	more differentially expressed than genes outside the set (competitive), but 
#	sometimes the analysis is for whether any gene in the selected gene set is 
#	differentially expressed (self-contained)
# - A technique called permutation can be used with >= 7 samples per group, 
#	where samples are assigned randomly to either their correct group or the 
#	incorrect one, and the other parts of the GSEA are performed. If the gene 
#	set is truly differentially expressed, when the samples are scrambled into
#	random groups, they shouldn't be differentially expressed anymore. Permuting
#	can also be done at the gene level, although this can lead to more false 
#	false positives

### Inputs ###
previous_gsea_results_path <- 
	'/media/amundy/Windows/bio/glioblastoma_p_selectin/GSEA_Results_Arianna.xlsx'
deg_df_path <- 
	'/media/amundy/Windows/bio/glioblastoma_p_selectin/Gy10_vs_Gy0_DBTRG_log2FC_1_padj_0.05.csv'
non_msig_gene_lists_path <- 
	'/media/amundy/Windows/bio/glioblastoma_p_selectin/senescence_genelists.csv'

gost_plot_up_path <- 
	'/media/amundy/Windows/bio/glioblastoma_p_selectin/outputs/gene_ontology_plot_up.html'
gost_plot_up_path <- 
	'/media/amundy/Windows/bio/glioblastoma_p_selectin/outputs/gene_ontology_plot_down.html'
gsea_table_path <- 
	'/media/amundy/Windows/bio/glioblastoma_p_selectin/outputs/gsea_table.xlsx'
gsea_plot_path <- 
	'/media/amundy/Windows/bio/glioblastoma_p_selectin/outputs/gsea_plot.pdf'
gsea_bubble_plot_path <- 
	'/media/amundy/Windows/bio/glioblastoma_p_selectin/outputs/gsea_bubble_plot.pdf'


# build dataframe with gene sets of interest
msig_hallmarks <- c('HALLMARK_INFLAMMATORY_RESPONSE',
					'HALLMARK_P53_PATHWAY')
# msig_hallmarks <- c(
# 	'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
# 	'HALLMARK_ALLOGRAFT_REJECTION',
# 	'HALLMARK_COAGULATION',
# 	'HALLMARK_INFLAMMATORY_RESPONSE',
# 	'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
# 	'HALLMARK_P53_PATHWAY',
# 	'HALLMARK_IL6_JAK_STAT3_SIGNALING',
# 	'HALLMARK_KRAS_SIGNALING_UP',
# 	'HALLMARK_INTERFERON_GAMMA_RESPONSE',
# 	'HALLMARK_APOPTOSIS',
# 	'HALLMARK_COMPLEMENT',
# 	'HALLMARK_KRAS_SIGNALING_DN',
# 	'HALLMARK_ANGIOGENESIS',
# 	'HALLMARK_XENOBIOTIC_METABOLISM',
# 	'HALLMARK_APICAL_SURFACE',
# 	'HALLMARK_APICAL_JUNCTION',
# 	'HALLMARK_HYPOXIA',
# 	'HALLMARK_TGF_BETA_SIGNALING',
# 	'HALLMARK_UV_RESPONSE_UP',
# 	'HALLMARK_MYOGENESIS',
# 	'HALLMARK_E2F_TARGETS',
# 	'HALLMARK_G2M_CHECKPOINT',
# 	'HALLMARK_MITOTIC_SPINDLE',
# 	'HALLMARK_MYC_TARGETS_V1',
# 	'HALLMARK_SPERMATOGENESIS',
# 	'HALLMARK_ANDROGEN_RESPONSE',
# 	'HALLMARK_MTORC1_SIGNALING'
# )
msig_gene_sets <- msigdbr(species='Homo sapiens', category='H')
msig_gene_sets <- dplyr::filter(msig_gene_sets, gs_name %in% msig_hallmarks)
msig_gene_sets <- dplyr::select(msig_gene_sets, gs_name, gene_symbol)

other_gene_sets <- read_csv(non_msig_gene_lists_path, col_types=cols(.default = 'c'))
other_gene_sets <- dplyr::select(other_gene_sets, Senescence_UP, SASP)
other_gene_sets <- pivot_longer(other_gene_sets,
								cols = c(Senescence_UP, SASP),
								names_to = "gs_name",
								values_to = "gene_symbol")
other_gene_sets <- drop_na(other_gene_sets)
other_gene_sets <- other_gene_sets[order(other_gene_sets$gs_name),]
gene_sets <- rbind(msig_gene_sets, other_gene_sets)

# read differentially expressed genes df and subset to genes of interest
deg_df <- read_csv(deg_df_path, col_types=cols("Row.names"="c"))
# subset to genes in gene sets
deg_df <- as_tibble(merge(deg_df, gene_sets, by.x = 'Row.names',
					   by.y='gene_symbol'))
# TODO figure out why this was wrong
# test2 <- dplyr::filter(deg_df, Row.names %in% gene_sets$gene_symbol)
deg_df_up <- dplyr::filter(deg_df, log2FoldChange >= 0)
deg_df_down <- dplyr::filter(deg_df, log2FoldChange < 0)

# Gene Ontology analysis
gost_res_up <- gost(deg_df_up$Row.names, organism = "hsapiens", correction_method = "fdr")
gost_res_down <- gost(deg_df_down$Row.names, organism = "hsapiens", correction_method = "fdr")
# TODO figure out how to save these to html programatically
# need to export these to html using r studio viewer
gostplot(gost_res_up, interactive = T, capped = T) 
gostplot(gost_res_down, interactive = T, capped = T) 

# build a named vector to pass to the GSEA function
gsea_input <- deg_df$log2FoldChange
names(gsea_input) <- as.character(deg_df$Row.names)
gsea_input <- sort(gsea_input, decreasing = TRUE)

# competitive GSEA with gene level permutations
gsea_res <- GSEA(gsea_input, TERM2GENE=gene_sets)
gsea_df <- as_tibble(gsea_res@result)
print(gsea_df)
print(gsea_df$Description)
write.xlsx(gsea_df, gsea_table_path)

# TODO one plot for each hallmark
pdf(gsea_plot_path)
gseaplot2(gsea_res, geneSetID = gsea_df$ID,
		  title = gsea_df$Description)
dev.off()

# bubble plot
# - bubble size: number of genes in the gene set
# - color: enrichment score
# - transparency: -log10 adjusted p value
bubble_plot_df <- mutate(gsea_df,
				  phenotype = case_when(NES > 0 ~ "intervention",
				  					  NES < 0 ~ "control"))
pdf(gsea_bubble_plot_path)
ggplot(bubble_plot_df, aes(x=phenotype, y=ID)) +
	geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
	scale_color_gradient(low="blue", high="red") +
	theme_bw()
dev.off()
