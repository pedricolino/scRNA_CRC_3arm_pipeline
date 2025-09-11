import sys as sys

print("Step 1: Imported sys")

# read in arguments
import argparse as parser

print("Step 2: Imported argparse")

parser = parser.ArgumentParser(description="Run CellUntangler on a dataset")
parser.add_argument("--use_ensembl_ids", action="store_true", help="Use Ensembl IDs instead of gene symbols")
parser.add_argument("--input_adata_path", type=str, required=True, help="Path to the input AnnData file")
parser.add_argument("--output_path", type=str, default="./", help="Path to save the output files")
args = parser.parse_args()

use_ensembl_ids = args.use_ensembl_ids
input_adata_path = args.input_adata_path
output_path = args.output_path

print("Step 3: Parsed arguments")

# print(sys.path)
# If running notebook outside of CellUntangler directory, add it to the path
sys.path.append('resources/CellUntangler/')
sys.path.append('resources/CellUntangler/')

print("Step 4: Updated sys.path")

import os

import torch

from src.data.umi_data import UMIVaeDataset
from src.celluntangler import utils
from src.celluntangler.models import Trainer
from src.celluntangler.models.nb_vae import NBVAE
import numpy as np

import pandas as pd
import scanpy as sc

print("Step 5: Imported required modules")

input_path = "/data/cephfs-1/home/users/cemo10_c/project_symlinks/crc/scratch/scRNA_CRC_3arm_pipeline/results/merged/adata_counts_scAutoQC.h5ad"
adata = sc.read_h5ad(input_path)

print("Step 6: Loaded AnnData")

# if layer log1p is not present, create it
if 'log1p_norm_of_counts' not in adata.layers:
	sc.pp.normalize_total(adata, target_sum=1e4)
	sc.pp.log1p(adata)
	adata.layers['log1p_norm_of_counts'] = adata.X.copy()
	print("Step 7: Created log1p_norm_of_counts layer")
else:
	adata.X = adata.layers['log1p_norm_of_counts'].copy()
	print("Step 7: Used existing log1p_norm_of_counts layer")

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = "sample")
print("Step 8: Computed highly variable genes")

adata.X = adata.layers["counts"].copy()
print("Step 9: Set adata.X to counts layer")

cell_cycle_genes_path = "resources/CellUntangler/cell_cycle_genes/celluntangler_human_cell_cycle_genes.tsv"
cell_cycle_genes = pd.read_csv(cell_cycle_genes_path, header=None, sep="\t")
print("Step 10: Loaded cell cycle genes")

if use_ensembl_ids:
	# add matching gene name from ensembl2symbol to adata.var
	ensembl2_symbol = pd.read_csv("resources/ensembl2genes_cellranger.tsv", sep="\t", header=None, index_col=0)
	# column one is ensembl id, column two is gene symbol
	ensembl2_symbol = ensembl2_symbol[1].to_dict()
	adata.var['gene_name'] = adata.var_names.map(ensembl2_symbol)
	print("Step 11: Mapped Ensembl IDs to gene names")

# Genes involved with the cell cycle present in adata
contained_cell_cycle_genes = adata.var["gene_name"].isin(cell_cycle_genes[0])
kept_genes = adata.var["highly_variable"] | contained_cell_cycle_genes
print("Number of kept genes after filtering for highly variable and cell cycle genes:", np.sum(kept_genes))
adata = adata[:, kept_genes]
print("Step 12: Filtered genes for highly variable and cell cycle genes")

# Genes involved with the cell cycle present in adata
contained_cell_cycle_genes = adata.var["gene_name"].isin(cell_cycle_genes[0])
print(f"Number of cell cycle genes present in adata: {np.sum(contained_cell_cycle_genes)}")
cycle_gene_indices = np.where(contained_cell_cycle_genes)[0]
non_cycle_gene_indices = np.where(~contained_cell_cycle_genes)[0]
rearranged_indices = np.hstack((cycle_gene_indices, non_cycle_gene_indices))
adata = adata[:, rearranged_indices]
adata.uns["new_gene_ordering"] = rearranged_indices
print("Step 13: Rearranged gene order")

from src.celluntangler.models.get_config import get_config

config = get_config()
print("Step 14: Loaded config")

config.model_name = "r2, e10"
config.seed = 68715
config.init = "custom"

x = adata.X.todense().astype(np.double)
chemistry_batch = adata.obs['sample'].rank(method='dense').astype(int).values.reshape(-1, 1)-1
# donor_batch = adata.obs['Donor'].rank(method='dense').astype(int).values.reshape(-1, 1)-1
# batch = np.hstack((chemistry_batch, donor_batch))
batch = chemistry_batch
# y holds the batch vector for the dataset
y = batch

# config.n_batch = [3, 12]
config.n_batch = [15]

print("Step 15: Prepared input data and batch info")

# Create the dataset and separate it into training and test sets
in_dim = x.shape[1]
print(f"in_dim = {in_dim}")
batch_size = config.batch_size
dataset = UMIVaeDataset(batch_size=batch_size, in_dim=in_dim)
print(f"dataset.in_dim = {dataset.in_dim}")
# Create the dataset loaders
train_loader = dataset.create_loaders(x, y, seed=config.seed)
print("Step 16: Created dataset and loaders")

# A mask that is 1 where the gene is a cell cycling gene and 0 otherwise
mask_cyc = np.zeros(adata.n_vars)
mask_cyc[adata.var["gene_name"].isin(cell_cycle_genes[0])] = 1

# A mask that is 0 where the gene is a cell cycling gene and 1 otherwise
mask_all = np.ones(adata.n_vars)
mask_all[adata.var["gene_name"].isin(cell_cycle_genes[0])] = 0

mask = torch.tensor([mask_cyc, mask_all])
print("Step 17: Created masks for cell cycle genes")

for row in mask:
	print(row)
	print(torch.sum(row).item())

torch.set_default_dtype(torch.float64)
print("Step 18: Set torch default dtype")

use_gpu = True
if use_gpu:
	print("Using GPU")
	config.device = torch.device("cuda")
else:
	print("Using CPU")
	config.device = torch.device("cpu")
print("Step 19: Set device")

# The path to save the intermediate embeddings to
epoch_embeddings_save_path = "~/scratch/celluntangler"
# Create the directory if it does not exist
if not os.path.exists(epoch_embeddings_save_path):
	os.makedirs(epoch_embeddings_save_path)
print("Step 20: Created embeddings save path")

visualize_information={}
# The epochs to save the intermediate embeddings for
visualize_information["epochs"]=[i for i in range(0, 500, 50)]
visualize_information["x"]=x
visualize_information["y"]=y
visualize_information["embeddings_save_path"]=epoch_embeddings_save_path
visualize_information["model_name"]=config.model_name
visualize_information["device"]=config.device
print("Step 21: Set visualize_information")

config

if config.seed:
	print(config.seed)
	torch.manual_seed(config.seed)
	np.random.seed(config.seed)
	np.random.default_rng(config.seed)
print("Step 22: Set random seeds")

model_name = config.model_name
components = utils.parse_components(model_name, config.fixed_curvature)
print("Step 23: Parsed model components")

model = NBVAE(h_dim=config.h_dim,
			  components=components,
			  mask=mask,
			  dataset=dataset,
			  config=config).to(config.device)
print("Step 24: Created NBVAE model")

trainer = Trainer(model)
print("Step 25: Created Trainer")

optimizer = trainer.build_optimizer(learning_rate=config.learning_rate,
										fixed_curvature=config.fixed_curvature,
										use_adamw=config.use_adamw,
										weight_decay=config.weight_decay)
print("Step 26: Built optimizer")

betas = utils.linear_betas(config.start,
                           config.end,
                           end_epoch=config.end_epoch,
                           epochs=config.epochs)

print("Step 27: Training now")

trainer.train_epochs(optimizer=optimizer,
                       train_data=train_loader,
                       betas=betas,
                       likelihood_n=0,
                       max_epochs=config.max_epochs,
                           visualize_information=visualize_information)


embeddings_save_path = output_path
# create the directory if it does not exist
if not os.path.exists(embeddings_save_path):
	os.makedirs(embeddings_save_path)
	
a = trainer.model(torch.log1p(torch.tensor(x,device=config.device)), torch.tensor(y,device=config.device))
np.savetxt(os.path.join(embeddings_save_path, f'{model_name}_all_encode_v63_z_params.txt'), a[4].detach().to(torch.device("cuda")).numpy())

print("Step 27: Saved embeddings")
