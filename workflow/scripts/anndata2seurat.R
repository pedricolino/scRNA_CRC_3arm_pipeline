library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "Input h5ad file"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "Output Seurat RDS file"),
  make_option(c("-c", "--count_layer"), type = "character", default = "counts", help = "Layer name for counts"),
  make_option(c("-n", "--norm_layer"), type = "character", default = "log1p_norm_of_counts", help = "Layer name for normalized data"),
  make_option(c("-x", "--layer_in_X"), type = "character", default = "counts", help = "Which layer is in adata.X, 'counts' or 'data' (=normalized data)")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

suppressPackageStartupMessages({library(Seurat)})
if (!requireNamespace("schard", quietly = TRUE)) devtools::install_github("cellgeni/schard")

srt = schard::h5ad2seurat(opt$input)

if (opt$layer_in_X == "counts") {
	srt[['RNA']]$data = schard::h5ad2Matrix(opt$input, opt$norm_layer)
} else if (opt$layer_in_X == "data") {
	srt[['RNA']]$counts = schard::h5ad2Matrix(opt$input, opt$count_layer)
} else {
	stop("Invalid value for --layer_in_X. Use 'counts' or 'data'.")
}

saveRDS(srt, file = opt$output)