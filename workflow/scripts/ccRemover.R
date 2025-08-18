library("optparse")

# input srt file
# output srt file
# use_ensembl_ids

option_list <- list(
  make_option(c("-i", "--input_file"), type = "character", default = NULL,
			  help = "Input file path", metavar = "character"),
  make_option(c("-o", "--output_file"), type = "character", default = NULL,
			  help = "Output file path", metavar = "character"),
make_option(c("-u", "--use_ensembl_ids"), action = "store_true", default = FALSE,
			help = "Use Ensembl IDs instead of gene symbols for cell cycle gene indexing")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

suppressPackageStartupMessages({
	library(Seurat)
	library(ccRemover)	
})

srt = readRDS(opt$input_file)

cell_cycle_gene_indices <- gene_indexer(rownames(srt), species = "human", name_type = ifelse(opt$use_ensembl_ids, "ensembl", "symbols"))
length(cell_cycle_gene_indices)

if_cc <- rep(FALSE,nrow(srt)) 
if_cc[cell_cycle_gene_indices] <- TRUE
summary(if_cc)

srt = ScaleData(srt)

dat <- list(x = GetAssayData(srt, 'RNA', "scale.data"), if_cc=if_cc)

xhat <- ccRemover(dat, bar=T)
xhat <- xhat + mean_gene_exp

saveRDS(xhat, opt$output_file)