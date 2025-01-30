import glob
import re
import scanpy as sc
import anndata as ad

# get input files from snakemake
input_files = snakemake.input

# extract sample name
sample_ids = [match.group() for path in sample_paths if (match := re.search("CE[a-zA-Z0-9_]*", path))]

# remove CE_SC_ prefix
sample_ids = [re.sub("CE_SC_", "", name) for name in sample_ids]
sample_ids = [re.sub("5FU_", "", name) for name in sample_ids]

# from sample names split by _ and take the first element as treatment and the second as week number
sample_treatments = [name.split("_")[0] for name in sample_ids]
sample_weeks = [name.split("_")[1] for name in sample_ids]
sample_treatments = ["Control" if treatment == "C" else treatment for treatment in sample_treatments]

# create a samples dictionary
samples = dict(zip(sample_ids, zip(sample_paths, sample_treatments, sample_weeks)))

adatas = {}

for sample_id, (sample_path, treatment, week) in samples.items():
    sample_adata = sc.read_h5ad(sample_path)
    sample_adata.obs['treatment'] = treatment
    sample_adata.obs['week'] = week
    adatas[sample_id] = sample_adata
    
adata = ad.concat(adatas, label="sample")
adata.obs_names_make_unique()

# write the merged anndata object to file
adata.write(snakemake.output)