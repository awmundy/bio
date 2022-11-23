suppressPackageStartupMessages({
	library(tidyverse)
  library(dplyr)
  #functions and methods for Gene Set Enrichment Analysis (i.e. functional enrichment analysis)
  library(GSEABase)
  # required by GSEABase
	library(Biobase) 
  # Gene Set Variation Analysis using an unsupervised method
	library(GSVA) 
  # tools for functional enrichment using g:Profiler
	library(gprofiler2) 
  # GSEA tool/plot builder
	library(clusterProfiler)
  # download msigdb gene lists
	library(msigdbr) 
  # enrichment plots
	library(enrichplot) 
  library(EnsDb.Hsapiens.v86)
  library(htmlwidgets)
}) 

write_gost_plot <- function(deg_df, out_path) {
  # builds and writes an interactive plot showing functional 
  # enrichment by various categories
  # deg_df: dataframe of differentially expressed genes
  unique_genes <- unique(deg_df$Row.names)
  
  # TODO determine what correction method is typically used and use that
  # default args return significantly enriched gene sets at a 0.05 threshold
  gost_res <- gost(unique_genes,
                   organism = "hsapiens",
                   correction_method = "fdr")
  
  out_plot <- gostplot(gost_res, interactive = T, capped = T) 
  htmlwidgets::saveWidget(out_plot, out_path)
}

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
	'/media/awmundy/Windows/bio/glioblastoma_p_selectin/GSEA_Results_Arianna.xlsx'
deg_df_path <- 
	'/media/awmundy/Windows/bio/glioblastoma_p_selectin/Gy10_vs_Gy0_DBTRG_log2FC_1_padj_0.05.csv'
non_msig_gene_lists_path <- 
  '/media/awmundy/Windows/bio/glioblastoma_p_selectin/senescence_genelists.csv'

gost_plot_up_path <- 
	'/media/awmundy/Windows/bio/glioblastoma_p_selectin/outputs/gene_ontology_plot_up.html'
gost_plot_down_path <- 
	'/media/awmundy/Windows/bio/glioblastoma_p_selectin/outputs/gene_ontology_plot_down.html'
gsea_table_path <- 
	'/media/awmundy/Windows/bio/glioblastoma_p_selectin/outputs/gsea_table.xlsx'
gsea_plot_path <- 
	'/media/awmundy/Windows/bio/glioblastoma_p_selectin/outputs/gsea_plot.pdf'
gsea_bubble_plot_path <- 
	'/media/awmundy/Windows/bio/glioblastoma_p_selectin/outputs/gsea_bubble_plot.pdf'


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

# load the gene ontology sets, and filter to the ones that are a part of 
# the hallmarks we care about
msig_gene_sets <- msigdbr(species='Homo sapiens', category='H')
msig_gene_sets <- dplyr::filter(msig_gene_sets, gs_name %in% msig_hallmarks)
msig_gene_sets <- dplyr::select(msig_gene_sets, gs_name, gene_symbol)

#TODO explain provenance of these gene sets and why we care about them
# load senescence related gene sets
other_gene_sets <- read_csv(non_msig_gene_lists_path, 
                            col_types=cols(.default = 'c'))
other_gene_sets <- dplyr::select(other_gene_sets, Senescence_UP, SASP)
other_gene_sets <- pivot_longer(other_gene_sets,
                                cols = c(Senescence_UP, SASP),
                                names_to = "gs_name",
                                values_to = "gene_symbol")
other_gene_sets <- drop_na(other_gene_sets)
other_gene_sets <- other_gene_sets[order(other_gene_sets$gs_name),]

# append gene sets together
gene_sets <- rbind(msig_gene_sets, other_gene_sets)

# read significantly differentially expressed genes df and subset 
# to genes of interest
deg_df <- read_csv(deg_df_path, col_types=cols("Row.names"="c"))
# subset to genes in gene sets
deg_df <- as_tibble(merge(deg_df, gene_sets, by.x = 'Row.names',
					   by.y='gene_symbol'))
# add label for whether the gene is in ensembl
all_genes <- as_tibble(genes(EnsDb.Hsapiens.v86))
deg_df$in_ensembl <- as.numeric(deg_df$Row.names %in% all_genes$gene_name)
genes_not_in_ensembl <- 
  list(unique(dplyr::filter(deg_df, deg_df$in_ensembl == 0)$Row.names))
if (length(genes_not_in_ensembl) > 0) {
  print(paste('Genes from given gene sets that are not in esembldb:', 
              genes_not_in_ensembl))}

# split out into upregulated and downregulated sets
deg_df_up <- dplyr::filter(deg_df, log2FoldChange >= 0)
deg_df_down <- dplyr::filter(deg_df, log2FoldChange < 0)

write_gost_plot(deg_df_up, gost_plot_up_path)
write_gost_plot(deg_df_down, gost_plot_down_path)

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
