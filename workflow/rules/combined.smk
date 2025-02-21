# if this rule does not work with the script, then do it with the merge_anndata_samples.ipynb notebook instead
rule merge_anndata_samples:
    input: 
        samples = expand('results/per_sample/{sample}/adata.h5ad', sample=samples.index),
        checkpoint = expand('results/per_sample/{sample}/checkpoints/annotation_w_{count_layer}.done', sample=wc, count_layer=config['count_layer'])
    output: 'results/merged/adata.h5ad'
    benchmark: 'benchmarks/merge_anndata_samples.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (160 * attempt**0.5), # square root to avoid overdoing
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "per_sample" + env_suffix
    shell:
        "python workflow/scripts/merge_anndata_samples.py -o {output} -i {input.samples}"


rule subsample:
    input: 'results/merged/adata.h5ad'
    output: 'results/subset/adata.h5ad'
    benchmark: 'benchmarks/subsample.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (100 * attempt**0.5),
        runtime=lambda wildcards, attempt: 1*30 if attempt == 1 else 4*60
    conda: env_prefix + "per_sample" + env_suffix
    params: fraction = config['subset']['fraction'] # what fraction of the data to keep
    shell:
        "python -c 'import scanpy as sc; adata = sc.read_h5ad(\"{input}\"); sc.pp.subsample(adata, {params.fraction}); adata.write(\"{output}\")'"


rule dim_reduc_concatenated_ds:
    input: 'results/' + filename + '/adata.h5ad'
    output: 'results/' + filename + 
            '/dim_reduc_potential_cc_removal' +
            '/count_layer__{count_layer}' +
            '/normalization__{normalization}' +
            '/scale_data_before_pca__{scale_data_before_pca}' +
            '/genes_for_pca__{genes_for_pca}' +
            '/cc_method__{cc_method}' +
            '/pca_n_components__{pca_n_components}' +
            '/umap_n_neighbors__{umap_n_neighbors}' +
            '/checkpoint.done'
    benchmark: 'benchmarks/dim_reduc_potential_cc_removal/' + filename +
                '/count_layer__{count_layer}' +
                '__{normalization}' +
                '__{scale_data_before_pca}' +
                '__{genes_for_pca}' +
                '__{cc_method}' +
                '__{pca_n_components}' +
                '__{umap_n_neighbors}' +
                '.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (100 * attempt),
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "scpa" + env_suffix
    shell:
        "papermill "
        "workflow/notebooks/dim_reduc_cc_removal.ipynb "
        'results/' + filename + 
            '/dim_reduc_potential_cc_removal' +
            '/count_layer__{wildcards.count_layer}' +
            '/normalization__{wildcards.normalization}' +
            '/scale_data_before_pca__{wildcards.scale_data_before_pca}' +
            '/genes_for_pca__{wildcards.genes_for_pca}' +
            '/cc_method__{wildcards.cc_method}' +
            '/pca_n_components__{wildcards.pca_n_components}' +
            '/umap_n_neighbors__{wildcards.umap_n_neighbors}' +
            '/notebook.ipynb '
        "-p input_file {input} "
        "-p scale_data_before_pca {wildcards.scale_data_before_pca} "
        "-p genes_for_pca {wildcards.genes_for_pca} "
        "-p normalization {wildcards.normalization} "
        "-p cc_method {wildcards.cc_method} "
        "-p pca_n_components {wildcards.pca_n_components} "
        "-p tsv_file results/" + filename + "/dim_reduc_potential_cc_removal/distances_between_week1_samples.tsv "
        "-p figures_folder results/" + filename + "/dim_reduc_potential_cc_removal/figures/ "
        '-p output_file results/' + filename + 
            '/dim_reduc_potential_cc_removal' +
            '/count_layer__{wildcards.count_layer}' +
            '/normalization__{wildcards.normalization}' +
            '/scale_data_before_pca__{wildcards.scale_data_before_pca}' +
            '/genes_for_pca__{wildcards.genes_for_pca}' +
            '/cc_method__{wildcards.cc_method}' +
            '/pca_n_components__{wildcards.pca_n_components}' +
            '/umap_n_neighbors__{wildcards.umap_n_neighbors}' +
            '/adata_without_layers.h5ad '
        "-p umap_n_neighbors {wildcards.umap_n_neighbors} && "
        "touch {output}"
