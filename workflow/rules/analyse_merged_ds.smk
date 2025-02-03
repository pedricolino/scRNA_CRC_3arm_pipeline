rule dim_reduc_concatenated_ds:
    input: 'results/preprocessing/' + filename + '.h5ad'
    output: 'results/analysis/' + filename + '/dim_reduc/checkpoints/{layer_to_use}__{scale_data_before_pca}__{genes_for_pca}__{cc_method}__{pca_n_components}__{umap_n_neighbors}.done'
    benchmark: 'benchmarks/dim_reduc/' + filename + '{layer_to_use}_{scale_data_before_pca}_{genes_for_pca}_{cc_method}_{pca_n_components}_{umap_n_neighbors}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "annotation_no_versions" + env_suffix
    shell:
        "papermill "
        "workflow/notebooks/merged_dim_reduc_diff_params.ipynb "
        "results/analysis/" + filename + "/dim_reduc/{wildcards.layer_to_use}__{wildcards.scale_data_before_pca}__{wildcards.genes_for_pca}__{wildcards.cc_method}__{wildcards.pca_n_components}__{wildcards.umap_n_neighbors}.ipynb "
        "-p input_file {input} "
        "-p scale_data_before_pca {wildcards.scale_data_before_pca} "
        "-p genes_for_pca {wildcards.genes_for_pca} "
        "-p layer_to_use {wildcards.layer_to_use} "
        "-p cc_method {wildcards.cc_method} "
        "-p pca_n_components {wildcards.pca_n_components} "
        "-p umap_n_neighbors {wildcards.umap_n_neighbors} && "
        "touch {output}"