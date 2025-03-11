# if this rule does not work with the script, then do it with the merge_anndata_samples.ipynb notebook instead
rule merge_anndata_samples:
    input: 
        samples = expand("results/per_sample/{sample}/adata_ready_for_merge_{{count_layer}}__{{qc_method}}.h5ad", sample=samples.index),
    output: 'results/merged/adata_{count_layer}_{qc_method}.h5ad'
    benchmark: 'benchmarks/merge/anndata_samples_w_{count_layer}_{qc_method}.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (95 * attempt), # square root to avoid overdoing
        runtime=lambda wildcards, attempt: 60 * attempt ** 2,
    conda: env_prefix + "preprocessing" + env_suffix
    shell:
        "python workflow/scripts/merge_anndata_samples.py -o {output} -i {input.samples}"


rule subsample_merged:
    input: 'results/merged/adata_{count_layer}_{qc_method}.h5ad'
    output: 'results/subset/adata_{count_layer}_{qc_method}.h5ad'
    benchmark: 'benchmarks/subsample/{count_layer}_{qc_method}.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (50 * attempt),
        runtime=lambda wildcards, attempt: 60 * attempt ** 2,
    conda: env_prefix + "preprocessing" + env_suffix
    params: fraction = config['subset']['fraction'] # what fraction of the data to keep
    shell:
        "python -c 'import scanpy as sc; adata = sc.read_h5ad(\"{input}\"); sc.pp.subsample(adata, {params.fraction}); adata.write(\"{output}\")'"


rule dim_reduc_concatenated_ds:
    input: 'results/' + filename + '/adata_{count_layer}_{qc_method}.h5ad'
    output: 'results/' + filename + 
            '/dim_reduc_potential_cc_removal' +
            '/count_layer__{count_layer}' +
            '/normalization__{normalization}' +
            '/scale_data_before_pca__{scale_data_before_pca}' +
            '/genes_for_pca__{genes_for_pca}' +
            '/cc_method__{cc_method}' +
            '/pca_n_components__{pca_n_components}' +
            '/umap_n_neighbors__{umap_n_neighbors}' +
            '/qc_method__{qc_method}' +
            '/notebook.done'
    benchmark: 'benchmarks/dim_reduc_potential_cc_removal/' + filename +
                '/count_layer__{count_layer}' +
                '__{normalization}' +
                '__{scale_data_before_pca}' +
                '__{genes_for_pca}' +
                '__{cc_method}' +
                '__{pca_n_components}' +
                '__{umap_n_neighbors}' +
                '__{qc_method}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (100 * attempt),
        runtime=lambda wildcards, attempt: 60 * attempt ** 2,
    conda: env_prefix + "scpca" + env_suffix
    shell:
        "papermill "
        "workflow/notebooks/dim_reduc_cc_removal.ipynb "
        '{output} '
        "-p input_file {input} "
        "-p scale_data_before_pca {wildcards.scale_data_before_pca} "
        "-p genes_for_pca {wildcards.genes_for_pca} "
        "-p normalization {wildcards.normalization} "
        "-p cc_method {wildcards.cc_method} "
        "-p pca_n_components {wildcards.pca_n_components} "
        "-p tsv_file results/" + filename + "/dim_reduc_potential_cc_removal/distances_between_week1_samples.tsv "
        "-p qc_method {wildcards.qc_method} "
        "-p count_layer {wildcards.count_layer} "
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
            '/qc_method__{wildcards.qc_method}' +
            '/adata_without_layers.h5ad '
        "-p umap_n_neighbors {wildcards.umap_n_neighbors}"

rule compare_parameter_options:
    input: 
        expand('results/' + filename + 
                '/dim_reduc_potential_cc_removal' +
                '/count_layer__{count_layer}' +
                '/normalization__{normalization}' +
                '/scale_data_before_pca__{scale_data_before_pca}' +
                '/genes_for_pca__{genes_for_pca}' +
                '/cc_method__{cc_method}' +
                '/pca_n_components__{pca_n_components}' +
                '/umap_n_neighbors__{umap_n_neighbors}' +
                '/qc_method__{qc_method}' +
                '/notebook.done', 
                count_layer=config['count_layer'],
                normalization=config['dim_reduc']['normalization'],
                scale_data_before_pca=config['dim_reduc']['scale_data_before_pca'],
                genes_for_pca=config['dim_reduc']['genes_for_pca'],
                cc_method=config['dim_reduc']['cc_method'],
                pca_n_components=config['dim_reduc']['pca_n_components'],
                umap_n_neighbors=config['dim_reduc']['umap_n_neighbors'],
                qc_method=config['qc_method'])
    output: 'results/' + filename + '/dim_reduc_potential_cc_removal/compare_parameter_options.ipynb'
    conda: env_prefix + "scpca" + env_suffix
    resources:
        mem=lambda wildcards, attempt: '%dG' % (4 * attempt), # reduce this before committing
        runtime=lambda wildcards, attempt: 60 * attempt ** 2,
    shell:
        "papermill "
        "workflow/notebooks/interactive_figure_viewer.ipynb "
        "-p input_folder results/" + filename + "/dim_reduc_potential_cc_removal/figures/ "
        "{output} "

rule save_adata_chosen_branch:
    input: chosen_branch + "/notebook.ipynb"
    output: 'results/chosen_branch/adata.h5ad'
    benchmark: 'benchmarks/save_adata_chosen_branch.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (100 * attempt),
        runtime=lambda wildcards, attempt: 60 * attempt ** 2,
    conda: env_prefix + "scpca" + env_suffix
    shell:
        "papermill " + chosen_branch + "/notebook.ipynb "
        "results/chosen_branch/notebook.ipynb "
        "-p output_file {output} "
        "-p save_adata True"

rule subsample_chosen_branch:
    input: 'results/chosen_branch/adata.h5ad'
    output: 'results/chosen_branch/adata_subset.h5ad'
    benchmark: 'benchmarks/subsample_chosen_branch.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (50 * attempt),
        runtime=lambda wildcards, attempt: 60 * attempt ** 2,
    conda: env_prefix + "preprocessing" + env_suffix
    params: fraction = config['subset']['fraction'] # what fraction of the data to keep
    shell:
        "python -c 'import scanpy as sc; adata = sc.read_h5ad(\"{input}\"); sc.pp.subsample(adata, {params.fraction}); adata.write(\"{output}\")'"

rule strip_adata:
    input: 'results/chosen_branch/adata.h5ad'
    output: 'results/chosen_branch/adata_stripped.h5ad'
    benchmark: 'benchmarks/strip_adata.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (50 * attempt),
        runtime=lambda wildcards, attempt: 60 * attempt ** 2,
    conda: env_prefix + "preprocessing" + env_suffix
    shell:
        "python workflow/scripts/strip_adata.py -i {input} -o {output}"

rule strip_adata_subset:
    input: 'results/chosen_branch/adata_subset.h5ad'
    output: 'results/chosen_branch/adata_subset_stripped.h5ad'
    benchmark: 'benchmarks/strip_adata_subset.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (20 * attempt),
        runtime=lambda wildcards, attempt: 60 * attempt ** 2,
    conda: env_prefix + "preprocessing" + env_suffix
    shell:
        "python workflow/scripts/strip_adata.py -i {input} -o {output}"

###-------------------- CLUSTERING ------------------------###

rule clustering_monocle3:
    input: 'results/chosen_branch/adata_stripped.h5ad'
    output: 'results/chosen_branch/clusters_leiden_monocle3.csv'
    benchmark: 'benchmarks/clustering_monocle3.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (200 * attempt),
        runtime=lambda wildcards, attempt: 60 * attempt ** 2,
    conda: env_prefix + "R" + env_suffix
    shell:
        "papermill "
        "workflow/notebooks/clustering_monocle3.ipynb "
        "results/chosen_branch/clustering_monocle3.ipynb "
        "-p input_file {input} "
        "-p output_csv {output} "

rule clustering_schist:
    input: 'results/chosen_branch/adata_stripped.h5ad'
    output: 'results/chosen_branch/clusters_schist.csv'
    benchmark: 'benchmarks/clustering_schist.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (400 * attempt),
        runtime=lambda wildcards, attempt: 24 * 60 * attempt,
    conda: env_prefix + "scpca" + env_suffix
    shell:
        "papermill "
        "workflow/notebooks/clustering_schist.ipynb "
        "results/chosen_branch/clustering_schist.ipynb "
        "-p input_file {input} "
        "-p output_csv {output} "

# rule clustering_scMiko_silhouette:
#     input: 'results/chosen_branch/adata_stripped.h5ad'
#     output: 
#         clusters_scmiko = 'results/chosen_branch/clusters_scMiko.csv',
#         clusters_silhouette = 'results/chosen_branch/clustering_silhouette.csv'
#     benchmark: 'benchmarks/clustering_scMiko_silhouette.tsv'
#     resources:
#         mem=lambda wildcards, attempt: '%dG' % (400 * attempt),
#         runtime=lambda wildcards, attempt: 24 * 60 * attempt,
#     conda: env_prefix + "R" + env_suffix
#     shell:
#         "papermill "
#         "workflow/notebooks/clustering_scMiko_silhouette.ipynb "
#         "results/chosen_branch/clustering_scMiko_silhouette.ipynb "
#         "-p input_file {input} "
#         "-p output_csv {output.clusters_scmiko} "

rule clustering_scMiko:
    input: 'results/chosen_branch/adata_stripped.h5ad'
    output: 'results/chosen_branch/clusters_scMiko.csv',
    benchmark: 'benchmarks/clustering_scMiko.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (400 * attempt),
        runtime=lambda wildcards, attempt: 24 * 60 * attempt,
    conda: env_prefix + "R" + env_suffix
    shell:
        "papermill "
        "workflow/notebooks/clustering_scMiko.ipynb "
        "results/chosen_branch/clustering_scMiko.ipynb "
        "-p input_file {input} "
        "-p output_csv {output} "

rule clustering_silhouette:
    input: 'results/chosen_branch/adata_stripped.h5ad'
    output: 'results/chosen_branch/clustering_silhouette.csv'
    benchmark: 'benchmarks/clustering_silhouette.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (400 * attempt),
        runtime=lambda wildcards, attempt: 24 * 60 * attempt,
    conda: env_prefix + "R" + env_suffix
    shell:
        "papermill "
        "workflow/notebooks/clustering_silhouette.ipynb "
        "results/chosen_branch/clustering_silhouette.ipynb "
        "-p input_file {input} "
        "-p output_csv {output} "

rule which_clustering_to_choose:
    input:
        adata = 'results/chosen_branch/adata_stripped.h5ad',
        clusters_monocle3 = 'results/chosen_branch/clusters_leiden_monocle3.csv',
        clusters_scmiko = 'results/chosen_branch/clusters_scMiko.csv',
        clusters_silhouette = 'results/chosen_branch/clustering_silhouette.csv'
    output: 'results/chosen_branch/checkpoint/which_clustering_to_choose.done'
    benchmark: 'benchmarks/which_clustering_to_choose.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (200 * attempt),
        runtime=lambda wildcards, attempt: 60 * attempt ** 2,
    conda: env_prefix + "scpca" + env_suffix
    shell:
        "papermill "
        "workflow/notebooks/which_clustering_to_choose.ipynb "
        "results/chosen_branch/which_clustering_to_choose.ipynb "
        "-p input_file {input.adata} "
        "-p clusters_monocle3 {input.clusters_monocle3} "
        "-p clusters_scmiko {input.clusters_scmiko} "
        "-p clusters_silhouette {input.clusters_silhouette} && "
        "touch {output}"

###-------------------- ANALYSES GROUPED BY CONDA ENVIRONMENT ------------------------###

rule analysis_in_R:
    input:
        adata = 'results/chosen_branch/adata_stripped.h5ad',
        clusters = config['chosen_clustering'],
        notebook = 'workflow/notebooks/{in_R}.ipynb'
    output:
        'results/chosen_branch/checkpoint/in_R/{in_R}.done'
    benchmark: 'benchmarks/{in_R}.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (200 * attempt),
        runtime=lambda wildcards, attempt: 60 * attempt ** 2
    threads: 14
    conda: env_prefix + "R" + env_suffix
    shell:
        "papermill "
        "{input.notebook} "
        "results/chosen_branch/{wildcards.in_R}.ipynb "
        "-p input_file {input.adata} "
        "-p clusters_file {input.clusters} && "
        "touch {output}"

rule analysis_with_decoupler:
    input:
        adata = 'results/chosen_branch/adata.h5ad',
        clusters = config['chosen_clustering'],
        notebook = 'workflow/notebooks/{with_decoupler}.ipynb'
    output:
        'results/chosen_branch/checkpoint/with_decoupler/{with_decoupler}.done'
    benchmark: 'benchmarks/{with_decoupler}.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (200 * attempt),
        runtime=lambda wildcards, attempt: 60 * attempt ** 2,
    conda: env_prefix + "decoupler" + env_suffix
    shell:
        "papermill "
        "{input.notebook} "
        "results/chosen_branch/{wildcards.with_decoupler}.ipynb "
        "-p input_file {input.adata} "
        "-p clusters_file {input.clusters} && "
        "touch {output}"

rule analysis_with_scpca:
    input:
        adata = 'results/chosen_branch/adata_stripped.h5ad',
        clusters = config['chosen_clustering'],
        notebook = 'workflow/notebooks/{with_scpca}.ipynb'
    output:
        'results/chosen_branch/checkpoint/with_scpca/{with_scpca}.done'
    benchmark: 'benchmarks/{with_scpca}.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (200 * attempt),
        runtime=lambda wildcards, attempt: 60 * attempt ** 2,
    conda: env_prefix + "scpca" + env_suffix
    shell:
        "papermill "
        "{input.notebook} "
        "results/chosen_branch/{wildcards.with_scpca}.ipynb "
        "-p input_file {input.adata} "
        "-p clusters_file {input.clusters} && "
        "touch {output}"

rule analysis_with_pathway_mod:
    input:
        adata = 'results/chosen_branch/adata_stripped.h5ad',
        clusters = config['chosen_clustering'],
        notebook = 'workflow/notebooks/{with_pathway_mod}.ipynb'
    output:
        'results/chosen_branch/checkpoint/with_pathway_mod/{with_pathway_mod}.done'
    benchmark: 'benchmarks/{with_pathway_mod}.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (200 * attempt),
        runtime=lambda wildcards, attempt: 60 * attempt ** 2,
    conda: env_prefix + "pathway_mod" + env_suffix
    shell:
        "papermill "
        "{input.notebook} "
        "results/chosen_branch/{wildcards.with_pathway_mod}.ipynb "
        "-p input_file {input.adata} "
        "-p clusters_file {input.clusters} && "
        "touch {output}"
