#----- Import all libraries
print('Import all libraries')

import sys as sys
import re
import os
import pandas as pd
import numpy as np
from scipy.io import mmread
import scipy.sparse as sp
import matplotlib.pyplot as plt
from IPython.display import Image
import scanpy as sc
from cnmf import cNMF
import shutil
import argparse
import multiprocessing


#----- Parse arguments
print('Parse arguments')

parser = argparse.ArgumentParser(description='Run cNMF analysis.')
parser.add_argument('--input_file', required=True, help='Path to input h5ad file')
parser.add_argument('--count_layer', required=True, help='Layer name in AnnData to use for counts')
parser.add_argument('--deviant_genes_file', required=True, help='Path to file with deviant genes')
args = parser.parse_args()

input_file = args.input_file
count_layer = args.count_layer
deviant_genes_file = args.deviant_genes_file

output_directory = re.sub('adata.*.h5ad', 'cNMF_' + count_layer, input_file)
if not os.path.exists(output_directory):
    os.mkdir(output_directory)
else:
    shutil.rmtree(output_directory)
    os.mkdir(output_directory)
    

#----- Set parameters
print('Set parameters')

run_name = 'cNMF'
seedx = 14 ## Specify a seed pseudorandom number generation for reproducibility
numiter=200 # Number of NMF replicates. Set this to a larger value ~200 for real data. We set this to a relatively low value here for illustration at a faster speed
numhvgenes=2000 ## Number of over-dispersed genes to use for running the actual factorizations
comp=np.arange(2,19) # Range of number of components to try
num_workers = 32  # Number of workers


#----- Run cNMF
print('Run cNMF')

## Initialize the cnmf object that will be used to run analyses
cnmf_obj = cNMF(output_dir=output_directory, name=run_name)

## Prepare the data, I.e. subset to 2000 high-variance genes, and variance normalize
cnmf_obj.prepare(counts_fn=input_file, components=comp, n_iter=numiter, seed=seedx, genes_file=deviant_genes_file)

## Specify workers for parallelization
def run_worker(worker_index, cnmf_obj, total_workers):
    cnmf_obj.factorize(worker_i=worker_index, total_workers=total_workers)

with multiprocessing.Pool(num_workers) as pool:
    pool.starmap(run_worker, [(i, cnmf_obj, num_workers) for i in range(num_workers)])


#----- Combine results and create k-selection plot
print('Combine results and create k-selection plot')

cnmf_obj.combine()
cnmf_obj.k_selection_plot(close_fig=False)

print('All done! Outputs are in ' + output_directory)