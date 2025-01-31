import argparse
import glob
import re
import anndata as ad
import scanpy as sc

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Concatenate multiple anndata input files and save to a merged output file (no data integration / batch correction).")

    # Add arguments for input files and output file
    parser.add_argument('-i', '--inputs', nargs='+', required=True, help='Input file(s)')
    parser.add_argument('-o', '--output', required=True, help='Output file')

    # Parse the arguments
    args = parser.parse_args()
    
    sample_paths = args.inputs
    print("Input files:", sample_paths)
    # Get input files from snakemake. 
    # Not done here bc module imports will fail despite working normally in outside of snakemake in the same conda environment.
    # input_files = snakemake.input

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
    print("Samples dictionary created.")
    adatas = {}

    for idx, (sample_id, (sample_path, treatment, week)) in enumerate(samples.items(), start=1):
        print(f"Processing file {idx} of {len(samples)}: {sample_path}")
        try:
            with open(sample_path, 'r') as file:
                sample_adata = sc.read_h5ad(sample_path)
                sample_adata.obs['treatment'] = treatment
                sample_adata.obs['week'] = week
                adatas[sample_id] = sample_adata
                print(f"Successfully read file: {sample_path}")
        except FileNotFoundError:
            print(f"Error: File '{sample_path}' not found.")
        except Exception as e:
            print(f"Error reading file '{sample_path}': {e}")
        
    print("Concatenating anndata objects.")
    adata = ad.concat(adatas, label="sample")
    adata.obs_names_make_unique()

    # write the merged anndata object to file
    try:
        print(f"Writing merged anndata to: {args.output}")
        with open(args.output, 'w') as output_file:
            adata.write(args.output)
        print(f"Combined content written to '{args.output}' successfully.")
    except Exception as e:
        print(f"Error writing to file '{args.output}': {e}")

if __name__ == "__main__":
    main()
