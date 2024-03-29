suppressPackageStartupMessages({
  library(tidyverse)
  library(zeallot)
  library(DropletUtils)
  library(Seurat)
  library(Matrix)
  library(scales)
  library(rjson)
  library(R2HTML)
  library(DT)
  library(arrow)
  library(tidyr)
  library(purrr)
  library(patchwork)
  library(scater)
  library(scran)
  # tensorflow installation instructions:
  # - install.packages tensorflow 
  # - then tensorflow::install_tensorflow(extra_packages='tensorflow-probability')
  #   - it may prompt about miniconda, if python is installed outside of 
  #     miniconda hit no and it should find that python and then install
  library(tensorflow) 
  # Note: cellassign not currently working (so no need for tensorflow either)
  # install instructions: https://github.com/irrationone/cellassign
  library(cellassign) 
  library(SingleR) 
  library(celldex) 
  library(pheatmap)
})

plain_format <- function(x,...) {
  # Format for labels
  # Lightly edited version of function in: 
  #   https://github.com/Sarah145/scRNA_pre_process
  format(x, ..., scientific = FALSE, drop0trailing = TRUE, big.mark = ",")
}

write_barcode_rank_plot <- function(barcode_ranks, barcode_rank_plot_path){
  # Lightly edited version of function in: 
  #   https://github.com/Sarah145/scRNA_pre_process
  
  plot_df <- data.frame(Rank=barcode_ranks$rank, Counts=barcode_ranks$total, 
                        cell=barcode_ranks$has_cell)
  
  # TODO add note for why we drop cells with duped counts
  keep <- !duplicated(barcode_ranks$total)
  plot_df <- plot_df[keep, ]
  plot_df <- subset(plot_df, plot_df$Counts >0)
  
  png(barcode_rank_plot_path, width = 700, height = 500)
  plt <-ggplot(plot_df, 
           aes(x=Rank, y=Counts, col=cell, alpha=cell)) + 
      geom_point(size=3) + 
      geom_hline(yintercept = barcode_ranks@metadata$knee, lty = 2, 
                 col = '#0972D5', size=1.5) +
      annotate("text", x=max(plot_df$Counts), y=barcode_ranks@metadata$knee+10000,
               label="Knee", color = "#0972D5", size=5.5) +
      geom_hline(yintercept = barcode_ranks@metadata$inflection, lty = 2,
                 col = '#09CFD5', size = 1.5) +
      annotate("text", x=max(plot_df$Counts),
               y=barcode_ranks@metadata$inflection+500, label="Inflection",
               color = "#09CFD5", size=5.5) +
      scale_x_log10(labels = plain_format,
                    breaks = trans_breaks("log10", function(x) round(10^x, 0))) +
      scale_y_log10(breaks = trans_breaks('log10', function(x) floor(10^x)),
                    labels = plain_format) +
      scale_color_manual(values = c('#8595A8', '#6406B6'), name = NULL,
                         labels = c("Background", "Cells")) +
      scale_alpha_manual(values = c(0.5,1)) +
      labs(x = 'Barcodes', y = 'UMI counts', title = 'Barcode Rank Plot') +
      guides(alpha = "none",
             colour = guide_legend(reverse = TRUE, override.aes=list(size = 5))) +
      theme_linedraw() +
      theme(plot.title = element_text(size=20, hjust=0.5, face = 'bold'),
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 15),
            legend.text = element_text(size=19),
            legend.background = element_rect(fill = 'transparent'),
            legend.position = c(0.15,0.15))
  print(plt)
  dev.off()
}

write_barcode_rank_plot_html <- function(seq_meta_df, cell_stats_df, 
                                         analysis_output_dir, report_label){
  # Lightly edited version of function in: 
  #   https://github.com/Sarah145/scRNA_pre_process
  system(paste0('base64 ', analysis_output_dir, '/barcode_rank.png > ', 
                analysis_output_dir, '/barcode_rank.txt'))
  b64_bc <- readChar(paste0(analysis_output_dir, '/barcode_rank.txt'), 
                     file.info(paste0(analysis_output_dir, 
                                      '/barcode_rank.txt'))$size)
  target <- HTMLInitFile(analysis_output_dir, 
                         filename=paste0(report_label, '_summary'))
  HTML('<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Roboto">', 
       file=target)
  HTML("<div class='title'>", file=target)
  HTML.title(' Pre-Processing Summary', HR=1, file = target)
  HTML("</div>", file = target)
  HTML.title(report_label, HR=2, file = target)
  HTML("<div id='wrapper'>", file=target)
  HTML("<div class='boxed' id='left' align='center'>", file=target)
  HTML.title('Sequencing/Alignment Stats', HR=3, file=target)
  HTML('<table style="width:100%">', file=target)
  HTML(paste('<tr> <td>', seq_meta_df$metric, '</td> 
             <td align="left">', seq_meta_df$value, 
             '</td> </tr>'), file=target)
  HTML('</table> <hr>', file=target)
  HTML.title('Cell Stats', HR=3, file=target)
  HTML('<table style="width:100%">', file=target)
  HTML(paste('<tr> <td>', cell_stats_df$metric, '</td> <td align="right">', 
             cell_stats_df$value, '</td> </tr>'), file=target)
  HTML('</table>', file=target)
  HTML("</div>", file = target)
  HTML("<div class='boxed' id='right' align='center'>", file=target)
  HTML(paste0("<img src='data:image/png;base64,", b64_bc, "' width=90%>"), file=target)
  HTML("</div>", file = target)
  HTML("</div>", file = target)
  HTML('<style type="text/css">
		.title {
    			background-color: #0972D5;
    			padding: 8px;
    			color: white;
    			position: fixed;
    			top: 0;
    			left: 0;
    			z-index: 999;
    			width: 100%;
		}
		.boxed {
  			border: 1px solid #868D96;
  			padding: 10px;
  			margin: 20px;
		}
		h1 {
			font-family: "Roboto";
			font-size: 33px;
		}
		h2 {
			font-family: "Roboto";
			font-size: 26px;
		}
		h3 {
			font-family: "Roboto";
			font-size: 18px;
		}
		#wrapper {
  			display: flex;
		}
		#left {
  			width: 50%;
		}
		#right {
  			width: 50%;
		}
		table tr:nth-child(even) {
  			background-color: #eee;
		}
		table tr:nth-child(odd) {
  			background-color: #fff;
		}
		table {
  			font-family: "Roboto";
			font-size: 20px;
			border: 1px solid #868D96;
		}
		#mathplayer{
  			height: 80px;
		}
		</style> </head>', file=target)
  HTMLEndFile()
}

plot_seurat_violin <- function(srt, seurat_data_labels) {
  # TODO consider renaming data objects to Total Molecules Detected (ncount) and 
  # Genes Detected (nfeature)
  # TODO for custom x and y labels would need to build violin plots manually
  # X axis is meaningless, the data is jittered to better see overlapping points
  VlnPlot(srt, seurat_data_labels, pt.size = 0.1) + 
  theme(axis.title.x = element_blank()) # removing Identity label at bottom
}

plot_gene_vs_molecule_count <- function(srt) {
  # Gene vs Molecule Count QA
  # A lot of points in lower left may be poor quality cells, a lot in upper 
  # right may mean multi cell drops or truly different cell populations with a 
  # lot going on in them
  ggplot(srt@meta.data, aes(nCount_RNA, nFeature_RNA)) +
    geom_point(alpha = 0.7, size = 0.5) +
    labs(x = "Molecule counts per cell", y = "Number of genes detected")
}

filter_seurat_metrics <- function(srt) {
  srt <- subset(srt,
                subset = nCount_RNA < 20000 &
                  nCount_RNA > 1000 &
                  nFeature_RNA > 1000 &
                  Mitochondria_pct < 20)
  
  return(srt)
}

build_differential_gene_expression_tibble <- function(srt) {
  # get differentially epressed genes dataframe for each cluster compared 
  # to all cells in the other clusters
  # 
  dge <- FindAllMarkers(srt, only.pos = TRUE, min.pct = 0.25,
                        logfc.threshold = 0.25)
  
  dge$pct_point_dif <- dge$pct.1 - dge$pct.2
  dge <- as_tibble(dge)
  dge <- dge %>% arrange(desc(avg_log2FC))
  
  return(dge)
}

plot_seurat_dge_datatable <- function(dge) {
  dt <- datatable(dge,
                  extensions = c('KeyTable', "FixedHeader"),
                  caption = 'Cluster Comparisons by Gene',
                  filter = 'top',
                  options = list(
                    keys = TRUE,
                    searchHighlight = TRUE,
                    pageLength = 10,
                    lengthMenu = c("10", "25", "50", "100")
                  )) 
  numeric_cols <- purrr::map_lgl(dt$x$data, is.numeric)
  pct_cols <- grepl('pct', colnames(dt$x$data))
  dt <- formatRound(dt, numeric_cols, digits = 3)
  dt <- formatPercentage(dt, pct_cols, digits = 2)
  print(dt)
}

plot_seurat_genes_of_interest <- function(srt, genes_of_interest) {
  # plot genes of interest on umap, color scale is gene expression
  FeaturePlot(srt, 
              reduction = "umap", 
              features = genes_of_interest,
              pt.size = 0.4, 
              min.cutoff = 'q10',
              label = TRUE)  
}

plot_seurat_top_genes_heatmap <- function(dge, srt){
  top_dge_by_cluster <- group_by(dge, cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
  DoHeatmap(srt, features=top_dge_by_cluster$gene)
}

plot_genes_of_interest_overlap_across_categories <- 
  function(genes_of_interest_sce) {
    pheatmap(genes_of_interest_sce)
  }

scale_and_cluster_srt <- function(srt) {
  # Create scaling factor object that will center counts around 0, with variance 1
  srt <- ScaleData(srt, verbose = FALSE)
  
  # add PCA and UMAP metrics to the seurat object
  srt <- RunPCA(srt, npcs = 40, verbose = FALSE)
  srt <- RunUMAP(srt, reduction = "pca", dims = 1:40, verbose=FALSE)
  
  srt <- FindNeighbors(srt, reduction = "pca", dims = 1:40, verbose=FALSE)
  srt <- FindClusters(srt, resolution = 0.5)
  
  return(srt)
}

create_and_prep_seurat_object <- function(exp_mtx, project_name) {
  
  srt <- CreateSeuratObject(counts = exp_mtx, min.cells = 3)
  # project.name gets added to plots and can't be removed
  srt@project.name <- project_name
  
  # normalize the counts stored in srt@assays$RNA@data@x
  srt <- NormalizeData(srt, verbose = FALSE)
  
  # identify outliers
  srt <- FindVariableFeatures(srt, verbose = FALSE)
  
  # get pct mitochondrial reads
  # stores in srt@meta.data
  # TODO what are other useful prefixes/patterns in the gene labels
  # TODO consider just using Gene Ontology ones for consistency
  srt[["Mitochondria_pct"]] <- PercentageFeatureSet(srt, pattern = "^MT-")
  
  srt <- scale_and_cluster_srt(srt)
  
  return(srt)
}

plot_cell_assign_umap <- function(sce, sce_subset) {
  # NOTE: Can't get this to work, some bug the cellassign devs attempted to fix,
  #   may be some sort of interaction with the tensorflow version
  
  # determine cell groupings
  fit <- cellassign(
    exprs_obj = sce_subset,
    marker_gene_info = genes_of_interest_dummies,
    s = sizeFactors(sce),
    shrinkage = TRUE,
    max_iter_adam = 50,
    min_delta = 2,
    verbose = TRUE)
  
  # add cellassign mappings to sce object
  sce$cell_type <- fit$cell_type
  
  plotUMAP(sce, colour_by = "cell_type")
}

add_cluster_predictions_based_on_public_datasets <- 
  function(sce, public_dataset_label) {
    # Downloads normalized expression data from public single 
    # cell rna seq datasets. Then applies cluster annotations 
    # to a singleCellExperiment object based on the downloaded data
    # NOTE: library also has functions to run across multiple datasets 
    # and take best score
    
    if (public_dataset_label == 'encode') {
      data <- BlueprintEncodeData(ensembl = FALSE)
    } else if (public_dataset_label == 'hpca') {
      data <- HumanPrimaryCellAtlasData(ensembl = FALSE)
    } else if (public_dataset_label == 'dice') {
      data <- DatabaseImmuneCellExpressionData(ensembl = FALSE)
    } else if (public_dataset_label == 'immgen') {
      data <- ImmGenData(ensembl = FALSE)
    } else if (public_dataset_label == 'monaco') {
      data <- MonacoImmuneData(ensembl = FALSE)
    } else if (public_dataset_label == 'mouserna') {
      data <- MouseRNAseqData(ensembl = FALSE)
    } else if (public_dataset_label == 'nover') {
      data <- NovershternHematopoieticData(ensembl = FALSE)
    } else {
      stop('invalid public_dataset_label')
    }
    
    predictions <- SingleR(
      test = sce,
      assay.type.test = 1,
      ref = data,
      labels = data$label.main
    )
    sce[["singler_labels"]] <- predictions$labels
    
    return(list(sce, predictions))
  }

get_cellranger_counts_barcodes_and_genes <-
  function(cellranger_counts_dir) {
    barcodes <-
      read.table(
        paste0(cellranger_counts_dir, 'barcodes.tsv'),
        sep = '\t',
        header = F
      )
    barcodes <- dplyr::rename(barcodes, barcode = V1)
    
    genes <-
      read.table(paste0(cellranger_counts_dir, 'genes.tsv'),
                 sep = '\t',
                 header = F)
    
    # must keep as matrix for memory reasons
    # columns are drops, rows are genes
    counts_raw <- readMM(paste0(cellranger_counts_dir, 'matrix.mtx'))
    
    # assign gene ids to each record, assign cell barcodes to each column
    rownames(counts_raw) <- genes[, 1]
    colnames(counts_raw) <- barcodes[, 1]
    
    return(list(counts_raw, barcodes, genes))
  }

build_drop_has_cell_probabilities <- function(counts, 
                                            probabilities_drop_has_cell_path) {
  
  if (file.exists(probabilities_drop_has_cell_path)) {
    cell_prob <- read_parquet(probabilities_drop_has_cell_path)
  } else {
    # get FDR adjusted (by BH) p-values for whether the drop contains a cell
    cell_prob <- as_tibble(emptyDrops(counts))
    # write out bc it takes some time to run
    write_parquet(prob, probabilities_drop_has_cell_path)
    
  }
  return(cell_prob)
}

build_has_cell_mask <- function(cell_prob) {
  # build mask for drops that have a cell (using arbitrary .05 adj p 
  # value cutoff)
  has_cell_msk <- cell_prob$FDR <= 0.05 
  
  # drops with low counts of UMI (indicators of RNA present) are 
  # considered empty and assigned null values. Mark FALSE here to 
  # remove them in counts
  has_cell_msk[is.na(has_cell_msk)] <- FALSE 
  
  return(has_cell_msk)
  
}

filter_counts_to_drops_with_cell <- function(cell_prob, counts_raw) {
  has_cell_msk <- build_has_cell_mask(cell_prob)
  
  # subset to just the columns (drops) that have cells
  counts <- counts_raw[, has_cell_msk] 
  
  return(counts)
}

build_barcode_ranks <- function(counts_raw, cell_prob) {
  # get data for a plot showing droplet ranks vs counts of 
  # umi quantity (logged for each)
  # - needs raw counts because non-cell drops have background 
  #   rna we want to include in plot
  barcode_ranks <- barcodeRanks(counts_raw)
  
  # add this so that plot can color datapoints based on presence of a cell
  has_cell_msk <- build_has_cell_mask(cell_prob)
  barcode_ranks$has_cell <- has_cell_msk
  
  return(barcode_ranks)
}

get_sequencing_metadata_dataframe <- function(kallisto_bus_output_dir) {
  # build sequencing metadata dataframe
  # load run info from JSON files produced by kallisto bus
  kb_meta <- c(fromJSON(file = paste0(kallisto_bus_output_dir, 'inspect.json')), 
               fromJSON(file = paste0(kallisto_bus_output_dir, 'run_info.json')))
  # get technology used from the kallisto bus call string, e.g. 10XV3
  single_cell_tech <- strsplit(strsplit(kb_meta$call, '-x ')[[1]][2], ' ')[[1]][1]
  kb_meta <- append(kb_meta, list('single_cell_tech'=single_cell_tech))
  kb_meta[['call']] = NULL
  kb_meta <- purrr::map_chr(kb_meta, as.character)
  kb_meta <- unlist(kb_meta)
  seq_meta_df <- data.frame(metric = names(kb_meta), value = kb_meta)
  seq_meta_df <- map_df(seq_meta_df, prettyNum, big.interval = 3,  big.mark = ",")
  
  return(seq_meta_df)  
}

get_cell_counts_stats_metadata_dataframe <- function(counts, counts_raw) {
  n_cells <- ncol(counts)
  # percentage of gene counts that are in drops with cells
  pct_counts_in_cells <- round((sum(counts)/sum(counts_raw))*100, 2) 
  median_total_counts_per_drop_w_cell <- median(colSums(counts))
  med_num_genes_expressed_per_cell <- 
    median(apply(counts, 2, function(x) sum(x >= 1)))
  tot_genes_detected_in_any_cell <- sum(rowSums(counts)>=1)
  
  cell_stats_df <- 
    data.frame(metric = c('Number of cells', 
                          'Pct counts in drops with cells', 
                          'Median counts per cell', 
                          'Median genes per cell', 
                          'Total genes detected in any cell'), 
               value = prettyNum(
                 c(n_cells, 
                   pct_counts_in_cells, 
                   median_total_counts_per_drop_w_cell,
                   med_num_genes_expressed_per_cell, 
                   tot_genes_detected_in_any_cell), 
                 big.mark = ','))
  return(cell_stats_df)
}

plot_seurat_pca_and_umap_clustering <- function(srt) {
  DimPlot(srt, reduction = "pca", split.by = "orig.ident", label = TRUE) +
    ggtitle('PCA')
  DimPlot(srt, reduction = "umap", split.by = "orig.ident", label = TRUE) +
    ggtitle('UMAP')
}

cellranger_counts_dir <- 
  '/media/awmundy/Windows/bio/diyt/single_cell_data/kallisto_bus_outputs/counts_unfiltered/cellranger/'
kallisto_bus_output_dir <- 
  '/media/awmundy/Windows/bio/diyt/single_cell_data/kallisto_bus_outputs/'
analysis_output_dir <- '/media/awmundy/Windows/bio/diyt/single_cell_data/analysis/'
probabilities_drop_has_cell_path <- 
  paste0(analysis_output_dir, 'probabilities_that_drop_has_cell.parquet')
cellranger_filtered_10x_count_dir <- paste0(analysis_output_dir, 'counts_filtered/')
barcode_rank_plot_path <- paste0(analysis_output_dir, 'barcode_rank.png')

c(counts_raw, barcodes, genes) %<-% 
  get_cellranger_counts_barcodes_and_genes(cellranger_counts_dir)

cell_prob <- 
  build_drop_has_cell_probabilities(counts_raw, probabilities_drop_has_cell_path)

counts <- filter_counts_to_drops_with_cell(cell_prob, counts_raw)

# write out cellranger style filtered counts output
write10xCounts(cellranger_filtered_10x_count_dir, gene.symbol = genes[,2],
               counts, overwrite=T) 

barcode_ranks <- build_barcode_ranks(counts_raw, cell_prob)

seq_meta_df <- get_sequencing_metadata_dataframe(kallisto_bus_output_dir)
cell_stats_df <- get_cell_counts_stats_metadata_dataframe(counts, counts_raw)

write_barcode_rank_plot(barcode_ranks, barcode_rank_plot_path)
write_barcode_rank_plot_html(seq_meta_df, cell_stats_df,
                             analysis_output_dir, 'barcode_plot_html_example')

# build an expression matrix from the cellranger files on disk
#   - rows are genes, columns are drops (barcodes)
exp_mtx <- Read10X(
  cellranger_filtered_10x_count_dir,
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

srt <- create_and_prep_seurat_object(exp_mtx,'insert_project_name')

plot_seurat_violin(srt, c("nCount_RNA", "nFeature_RNA", "Mitochondria_pct"))

#TODO determine correct filtering strategy (i.e. when to exclude outliers)
# srt <- filter_seurat_metrics(srt)

# TODO this is QC
plot_gene_vs_molecule_count(srt)

# plot clusters using umap and pca
plot_seurat_pca_and_umap_clustering(srt)

dge <- build_differential_gene_expression_tibble(srt)
plot_seurat_dge_datatable(dge)

genes_of_interest_seurat <- c("IGHM", 'CD79A')
plot_seurat_genes_of_interest(srt, genes_of_interest_seurat)
plot_seurat_top_genes_heatmap(dge, srt)

# TODO not currently working, look in to
# read in counts as singleCellExperiment object (using DropletUtils)
# sce <- read10xCounts(cellranger_filtered_10x_count_dir)

# convert seurat object to singleCellExperiment object
sce <- as.SingleCellExperiment(srt)
# for QA checks later
rowData(sce)$Symbol <- rownames(sce)

# create a list of genes of interest
# cell type specific markers here: http://biocc.hrbmu.edu.cn/CellMarker/
genes_of_interest_dummies <- list(
  Monocytes = c("CD14", "CD68"),
  `T cells` = c("CD2", "CD3D", "TRAC", "IL32", "CD3E", "PTPRC"),
  `NK cells` = c("GZMK", "KLRF1", "CCL3", "CMC1", "NKG7", "PTPRC"),
  `Plasma cells` = c("CD27", "IGHG1", "CD79A", "IGHG2", "PTPRC", "IGKC"),
  `Mature B cells` = c("MS4A1", "LTB", "CD52", "IGHD", "CD79A", "PTPRC", "IGKC"))

# make dummies matrix (rows=genes, cols=gene categories)
genes_of_interest_dummies <- marker_list_to_mat(genes_of_interest_dummies,
                                                include_other = FALSE)
plot_genes_of_interest_overlap_across_categories(genes_of_interest_dummies)

gene_in_sce_idx <- match(rownames(genes_of_interest_dummies), rowData(sce)$Symbol)
stopifnot(all(!is.na(gene_in_sce_idx)))

#subset to genes of interest
sce_subset <- sce[gene_in_sce_idx, ]
stopifnot(all.equal(rownames(genes_of_interest_dummies), rowData(sce_subset)$Symbol))

# add size factors to sce that are scaling factors used to normalize the data
#   NOTE: non-positive size factor warning can trigger even when all size 
#   factors are positive, for some reason
sce <- computeSumFactors(sce)

# cellassign not working
# plot_cell_assign_umap(sce, sce_subset)

c(sce, predictions) %<-% 
  add_cluster_predictions_based_on_public_datasets(sce, 'nover')
plotScoreHeatmap(predictions)
plotUMAP(sce, colour_by = "singler_labels")

