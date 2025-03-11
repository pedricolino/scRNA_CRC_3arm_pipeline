import scanpy as sc
import argparse

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Strip an anndata file saved in h5ad from most of its layers.")

    # Add arguments for input files and output file
    parser.add_argument('-i', '--input', required=True, help='Input file')
    parser.add_argument('-o', '--output', required=True, help='Output file')

    # Parse the arguments
    args = parser.parse_args()
    
    adata = sc.read(args.input)
    
    adata.X = adata.layers['counts'].copy()
    
    # which layers are there that are not 'counts', 'log1p_norm_of_counts', 'scran_normalization_of_counts'
    layers_to_del = [layer for layer in adata.layers.keys() if layer not in ['counts', 'log1p_norm_of_counts', 'scran_normalization_of_counts']]

    for l in layers_to_del:
        del adata.layers[l]

    del adata.raw

    adata.write(args.output)

if __name__ == "__main__":
    main()