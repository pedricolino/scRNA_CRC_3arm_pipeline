rule quality_control:
    input: lambda wildcards: samples.at[wildcards.sample, 'path']
    output: 
        adata = 'results/preprocessing/{sample}/adata.h5ad',
        checkpoint = 'results/preprocessing/{sample}/checkpoints/quality_control.done'
        # nb = 'results/preprocessing/{sample}/quality_control.ipynb' # excluded bc it causes issues with pipeline
    benchmark: 'benchmarks/quality_control/{sample}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "preprocessing" + env_suffix
    shell:
        "papermill "
            "workflow/notebooks/quality_control.ipynb "
            "results/preprocessing/{wildcards.sample}/quality_control.ipynb "
            "-p input_file {input} "
            "-p output_dir results/preprocessing/{wildcards.sample}/ && "
        "touch {output.checkpoint}"


rule normalize_soupX_counts:
    input: 
        checkpoint = 'results/preprocessing/{sample}/checkpoints/quality_control.done'
    output: 'results/preprocessing/{sample}/checkpoints/normalize_soupX_counts.done'
    benchmark: 'benchmarks/normalize_soupX_counts/{sample}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "preprocessing" + env_suffix
    params: 
        count_layer_to_use = "soupX_counts",
        adata = 'results/preprocessing/{sample}/adata.h5ad' # should not be input bc if SnakeMake's detection of change in input files
    shell:
        "papermill "
            "workflow/notebooks/normalization.ipynb "
            "results/preprocessing/{wildcards.sample}/normalize_{params.count_layer_to_use}.ipynb "
            "-p input_file {params.adata} "
            "-p count_layer_to_use {params.count_layer_to_use} && "
        "touch {output}"


use rule normalize_soupX_counts as normalize_raw_counts with:
    input: 
        checkpoint = 'results/preprocessing/{sample}/checkpoints/normalize_soupX_counts.done'
    output: 'results/preprocessing/{sample}/checkpoints/normalize_raw_counts.done'
    benchmark: 'benchmarks/normalize_raw_counts/{sample}.tsv'
    params: 
        count_layer_to_use = "counts",
        adata = 'results/preprocessing/{sample}/adata.h5ad'


rule feature_selection:
    input: 
        checkpoint1 = 'results/preprocessing/{sample}/checkpoints/normalize_raw_counts.done',
        checkpoint2 = 'results/preprocessing/{sample}/checkpoints/normalize_soupX_counts.done'
    output: 'results/preprocessing/{sample}/checkpoints/feature_selection_w_' + config["count_layer_to_use"] + '.done'
    benchmark: 'benchmarks/feature_selection/{sample}_w_' + config["count_layer_to_use"] + '.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (48 * attempt), # 32 is not enough for sample Conti_1
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "preprocessing" + env_suffix
    params: 
        count_layer_to_use = "soupX_counts",
        adata = 'results/preprocessing/{sample}/adata.h5ad'
    shell:
        "papermill "
            "workflow/notebooks/feature_selection.ipynb "
            "results/preprocessing/{wildcards.sample}/feature_selection.ipynb "
            "-p input_file {params.adata} "
            "-p count_layer_to_use {params.count_layer_to_use} && "
        "touch {output}"


rule dimensionality_reduction:
    input:
        checkpoint = 'results/preprocessing/{sample}/checkpoints/feature_selection_w_' + config["count_layer_to_use"] + '.done'
    output: 'results/preprocessing/{sample}/checkpoints/dimensionality_reduction_w_' + config["count_layer_to_use"] + '.done'
    benchmark: 'benchmarks/dimensionality_reduction/{sample}_w_' + config["count_layer_to_use"] + '.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "preprocessing" + env_suffix
    params: 
        count_layer_to_use = "soupX_counts",
        adata = 'results/preprocessing/{sample}/adata.h5ad'
    shell:
        "papermill "
            "workflow/notebooks/dimensionality_reduction.ipynb "
            "results/preprocessing/{wildcards.sample}/dimensionality_reduction_w_" + config["count_layer_to_use"] + ".ipynb "
            "-p input_file {params.adata} "
            "-p count_layer_to_use {params.count_layer_to_use} && "
        "touch {output}"


rule clustering_per_sample:
    input: 
        checkpoint = 'results/preprocessing/{sample}/checkpoints/dimensionality_reduction_w_' + config["count_layer_to_use"] + '.done'
    output: 'results/preprocessing/{sample}/checkpoints/clustering_w_' + config["count_layer_to_use"] + '.done'
    benchmark: 'benchmarks/clustering/{sample}_w_' + config["count_layer_to_use"] + '.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "preprocessing" + env_suffix
    params: 
        count_layer_to_use = "soupX_counts",
        adata = 'results/preprocessing/{sample}/adata.h5ad'
    shell:
        "papermill "
            "workflow/notebooks/clustering.ipynb "
            "results/preprocessing/{wildcards.sample}/clustering_w_" + config["count_layer_to_use"] + ".ipynb "
            "-p input_file {params.adata} "
            "-p count_layer_to_use {params.count_layer_to_use} && "
        "touch {output}"


rule annotate_per_sample:
    input: 
        checkpoint = 'results/preprocessing/{sample}/checkpoints/clustering_w_' + config["count_layer_to_use"] + '.done'
    output: 'results/preprocessing/{sample}/checkpoints/annotation_w_' + config["count_layer_to_use"] + '.done'
    benchmark: 'benchmarks/annotation/{sample}_w_' + config["count_layer_to_use"] + '.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "preprocessing" + env_suffix
    params: 
        count_layer_to_use = "soupX_counts",
        adata = 'results/preprocessing/{sample}/adata.h5ad'
    shell:
        "papermill "
            "workflow/notebooks/annotation.ipynb "
            "results/preprocessing/{wildcards.sample}/annotation_w_" + config["count_layer_to_use"] + ".ipynb "
            "-p input_file {params.adata} "
            "-p count_layer_to_use {params.count_layer_to_use} && "
        "touch {output}"


# if this rule does not work with the script, then do it with the merge_anndata_samples.ipynb notebook instead
rule merge_anndata_samples:
    input: 
        samples = expand('results/preprocessing/{sample}/adata.h5ad', sample=samples.index),
        checkpoint = expand('results/preprocessing/{sample}/checkpoints/annotation_w_' + config["count_layer_to_use"] + '.done', sample=wc)
    output: 'results/preprocessing/merged.h5ad'
    benchmark: 'benchmarks/merge_anndata_samples.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (160 * attempt**0.5), # square root to avoid overdoing
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "preprocessing" + env_suffix
    shell:
        "python workflow/scripts/merge_anndata_samples.py -o {output} -i {input.samples}"

rule subsample:
    input: 'results/preprocessing/merged.h5ad'
    output: 'results/preprocessing/subset.h5ad'
    benchmark: 'benchmarks/subsample.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (100 * attempt**0.5),
        runtime=lambda wildcards, attempt: 1*30 if attempt == 1 else 4*60
    conda: env_prefix + "preprocessing" + env_suffix
    params: fraction = config['subset']['fraction'] # what fraction of the data to keep
    shell:
        "python -c 'import scanpy as sc; adata = sc.read_h5ad(\"{input}\"); sc.pp.subsample(adata, {params.fraction}); adata.write(\"{output}\")'"



###-------------------- experimental, different quality control methods ------------------------###

rule scAutoQC:
    input: lambda wildcards: samples.at[wildcards.sample, 'path']
    output: 'results/preprocessing/{sample}/scAutoQC.h5ad'
        # Excluded because the pipeline should not choke on this side product:
        # nb = 'results/preprocessing/{sample}/scAutoQC.ipynb'
    benchmark: 'benchmarks/scAutoQC/{sample}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "sctk" + env_suffix
    shell:
        "papermill "
            "workflow/notebooks/sctk_scAutoQC.ipynb "
            "results/preprocessing/{wildcards.sample}/scAutoQC.ipynb "
            "-p input_file {input} "
            "-p output_dir results/preprocessing/{wildcards.sample}/"