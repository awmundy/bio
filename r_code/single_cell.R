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
})

plain_format <- function(x,...) {
  # Format for labels
  # Lightly edited version of function in: 
  #   https://github.com/Sarah145/scRNA_pre_process
  format(x, ..., scientific = FALSE, drop0trailing = TRUE, big.mark = ",")
}

bc_rank_plot <- function(stats, raw_cells, filt_cells, save){
  # Lightly edited version of function in: 
  #   https://github.com/Sarah145/scRNA_pre_process
  cells <- raw_cells %in% filt_cells
  keep <- !duplicated(stats$total)
  plot_df <- data.frame(Rx=stats$rank, Tx=stats$total, cell=cells)
  plot_df <- plot_df[keep, ]
  png(save, width = 700, height = 500)
  print({
    ggplot(subset(plot_df, plot_df$Tx >0), aes(x=Rx, y=Tx, col=cell, alpha=cell))+ 
      geom_point(size=3) + 
      geom_hline(yintercept = stats@metadata$knee, lty = 2, col = '#0972D5', size=1.5) +
      annotate("text", x=max(plot_df$Tx), y=stats@metadata$knee+10000, label="Knee", color = "#0972D5", size=5.5) +
      geom_hline(yintercept = stats@metadata$inflection, lty = 2, col = '#09CFD5', size = 1.5) +
      annotate("text", x=max(plot_df$Tx), y=stats@metadata$inflection+500, label="Inflection", color = "#09CFD5", size=5.5) +
      scale_x_log10(labels = plain_format, breaks = trans_breaks("log10", function(x) round(10^x, 0))) + 
      scale_y_log10(breaks = trans_breaks('log10', function(x) floor(10^x)), labels = plain_format)+
      scale_color_manual(values = c('#8595A8', '#6406B6'), name = NULL, labels = c("Background", "Cells")) +
      scale_alpha_manual(values = c(0.5,1)) +
      labs(x = 'Barcodes', y = 'UMI counts', title = 'Barcode Rank Plot') +
      guides(alpha = "none", colour = guide_legend(reverse = TRUE, override.aes=list(size = 5))) +
      theme_linedraw() + 
      theme(plot.title = element_text(size=20, hjust=0.5, face = 'bold'),
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 15),
            legend.text = element_text(size=19),
            legend.background = element_rect(fill = 'transparent'),
            legend.position = c(0.15,0.15))
  })
  dev.off()
}

print_HTML <- function(seq_stats, cell_stats, dir, sample_id){
  # Lightly edited version of function in: 
  #   https://github.com/Sarah145/scRNA_pre_process
  system(paste0('base64 ', dir, '/barcode_rank.png > ', dir, '/barcode_rank.txt'))
  b64_bc <- readChar(paste0(dir, '/barcode_rank.txt'), file.info(paste0(dir, '/barcode_rank.txt'))$size)
  target <- HTMLInitFile(dir, filename=paste0(sample_id, '_summary'))
  HTML('<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Roboto">', file=target)
  HTML("<div class='title'>", file=target)
  HTML.title(' Pre-Processing Summary', HR=1, file = target)
  HTML("</div>", file = target)
  HTML.title(sample_id, HR=2, file = target)
  HTML("<div id='wrapper'>", file=target)
  HTML("<div class='boxed' id='left' align='center'>", file=target)
  HTML.title('Sequencing/Alignment Stats', HR=3, file=target)
  HTML('<table style="width:100%">', file=target)
  HTML(paste('<tr> <td>', seq_stats$stat, '</td> <td align="right">', seq_stats$value, '</td> </tr>'), file=target)
  HTML('</table> <hr>', file=target)
  HTML.title('Cell Stats', HR=3, file=target)
  HTML('<table style="width:100%">', file=target)
  HTML(paste('<tr> <td>', cell_stats$stat, '</td> <td align="right">', cell_stats$value, '</td> </tr>'), file=target)
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
drop_has_cell_probabilities_path <- 
  '/media/awmundy/Windows/bio/diyt/single_cell_data/analysis/drop_has_cell_probabilities.parquet'

barcodes <- read.table(paste0(cellranger_counts_dir, 'barcodes.tsv'), sep = '\t', header = F)
genes <- read.table(paste0(cellranger_counts_dir, 'genes.tsv'), sep='\t', header=F)
# must keep as matrix for memory reasons
counts <- readMM(paste0(cellranger_counts_dir, 'matrix.mtx'))

# assign gene ids to each record, assign cell barcodes to each column
rownames(counts) <- genes[,1]
colnames(counts) <- barcodes[,1]

# # get FDR adjusted (by BH) p-values for whether the drop contains a cell
# prob <- as_tibble(emptyDrops(counts))
# # write out bc it takes some time to run
# write_parquet(prob, drop_has_cell_probabilities_path)

prob <- read_parquet(drop_has_cell_probabilities_path)
# subset to drops that have a cell (using arbitrary .05 adj p value cutoff)
has_cell_msk <- prob$FDR <= 0.05 
# drops with low counts of UMI (indicators of RNA present) are 
# considered empty and assigned null values. Mark FALSE here to 
# remove them in counts
has_cell_msk[is.na(has_cell_msk)] <- FALSE 
# subset to just the columns (drops) that have cells
counts <- counts[, has_cell_msk] 

# load run info from JSON files produced by Kb
kb_meta <- c(fromJSON(file = paste0(kallisto_bus_output_dir, 'inspect.json')), 
              fromJSON(file = paste0(kallisto_bus_output_dir, 'run_info.json')))
# get technology used from the kallisto bus call string, e.g. 10XV3
single_cell_tech <- strsplit(strsplit(kb$call, '-x ')[[1]][2], ' ')[[1]][1]

seq_meta_df <- as_tibble(kb_meta)
seq_meta_df <- map_df(seq_meta_df, prettyNum, big.interval = 3,  big.mark = ",")
 


