# Import data into R and filter out empty drops ----
# Begin by setting up a new RProject in the folder where you just processed your scRNA-seq data with Kb
suppressPackageStartupMessages({
  library(tidyverse)
  library(DropletUtils)
  library(Seurat) # a huge, powerful, and popular library for analyzing single cell genomic data
  library(Matrix)
  library(scales)
  library(rjson)
  library(R2HTML)
  library(DT)
  library(arrow)
  library(tidyr)
  library(purrr)
  library(patchwork)
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


cellranger_counts_dir <- 
  '/media/awmundy/Windows/bio/diyt/single_cell_data/kallisto_bus_outputs/counts_unfiltered/cellranger/'
kallisto_bus_output_dir <- 
  '/media/awmundy/Windows/bio/diyt/single_cell_data/kallisto_bus_outputs/'
analysis_output_dir <- '/media/awmundy/Windows/bio/diyt/single_cell_data/analysis/'
probabilities_drop_has_cell_path <- 
  paste0(analysis_output_dir, 'probabilities_that_drop_has_cell.parquet')
cellranger_10x_count_dir <- paste0(analysis_output_dir, 'counts_filtered/')
barcode_rank_plot_path <- paste0(analysis_output_dir, 'barcode_rank.png')

barcodes <- read.table(paste0(cellranger_counts_dir, 'barcodes.tsv'), sep = '\t', header = F)
barcodes <- dplyr::rename(barcodes, barcode=V1)
genes <- read.table(paste0(cellranger_counts_dir, 'genes.tsv'), sep='\t', header=F)
# must keep as matrix for memory reasons
# columns are drops, rows are genes
counts_raw <- readMM(paste0(cellranger_counts_dir, 'matrix.mtx'))

# assign gene ids to each record, assign cell barcodes to each column
rownames(counts_raw) <- genes[,1]
colnames(counts_raw) <- barcodes[,1]

# # get FDR adjusted (by BH) p-values for whether the drop contains a cell
# prob <- as_tibble(emptyDrops(counts))
# # write out bc it takes some time to run
# write_parquet(prob, probabilities_drop_has_cell_path)

# TODO how to identify and remove drops with multiple cells?
prob <- read_parquet(probabilities_drop_has_cell_path)
# subset to drops that have a cell (using arbitrary .05 adj p value cutoff)
has_cell_msk <- prob$FDR <= 0.05 
# drops with low counts of UMI (indicators of RNA present) are 
# considered empty and assigned null values. Mark FALSE here to 
# remove them in counts
has_cell_msk[is.na(has_cell_msk)] <- FALSE 
# subset to just the columns (drops) that have cells
counts <- counts_raw[, has_cell_msk] 
# add flag for if the drop has a cell
# TODO remove this if we're not using the barcodes object anymore
barcodes$has_cell <- has_cell_msk
# write out cellranger style filtered counts output
write10xCounts(cellranger_filtered_10x_count_dir, gene.symbol = genes[,2],
               counts, overwrite=T) 
# data for a plot showing droplet ranks vs counts (logged for each)
# - needs raw counts because non-cell drops have background rna we want to 
#   include in plot
barcode_ranks <- barcodeRanks(counts_raw)
barcode_ranks$has_cell <- has_cell_msk

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

n_cells <- ncol(counts)
# percentage of gene counts that are in drops with cells
pct_counts_in_cells <- round((sum(counts)/sum(counts_raw))*100, 2) 
median_total_counts_per_drop_w_cell <- median(colSums(counts))
med_num_genes_expressed_per_cell <- median(apply(counts, 
                                                 2, 
                                                 function(x) sum(x >= 1)))
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



# create barcode rank plot png
write_barcode_rank_plot(barcode_ranks, barcode_rank_plot_path)

# # output a HTML summary of the run
print_HTML(seq_meta_df, cell_stats_df, 
           analysis_output_dir, sample_id = NULL)

 


