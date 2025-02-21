rule quality_control:
    input: lambda wildcards: samples.at[wildcards.sample, 'path']
    output: 
        adata = 'results/per_sample/{sample}/adata.h5ad',
        checkpoint = 'results/per_sample/{sample}/checkpoints/quality_control.done'
        # nb = 'results/per_sample/{sample}/quality_control.ipynb' # excluded bc it causes issues with pipeline
    benchmark: 'benchmarks/quality_control/{sample}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "preprocessing" + env_suffix
    shell:
        "papermill "
            "workflow/notebooks/quality_control.ipynb "
            "results/per_sample/{wildcards.sample}/quality_control.ipynb "
            "-p input_file {input} "
            "-p output_dir results/per_sample/{wildcards.sample}/ && "
        "touch {output.checkpoint}"


rule normalize:
    input: 
        checkpoint = 'results/per_sample/{sample}/checkpoints/quality_control.done'
    output: 'results/per_sample/{sample}/checkpoints/normalize_{count_layer}.done'
    benchmark: 'benchmarks/normalize/{sample}_{count_layer}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "preprocessing" + env_suffix
    params: 
        adata = 'results/per_sample/{sample}/adata.h5ad' # should not be input bc if SnakeMake's detection of change in input files
    shell:
        "papermill "
            "workflow/notebooks/normalization.ipynb "
            "results/per_sample/{wildcards.sample}/normalize_{wildcards.count_layer}.ipynb "
            "-p input_file {params.adata} "
            "-p count_layer {wildcards.count_layer} && "
        "touch {output}"


rule feature_selection:
    input: 'results/per_sample/{sample}/checkpoints/normalize_{count_layer}.done'
    output: 'results/per_sample/{sample}/checkpoints/feature_selection_w_{count_layer}.done'
    benchmark: 'benchmarks/feature_selection/{sample}_w_{count_layer}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (48 * attempt), # 32 is not enough for sample Conti_1
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "preprocessing" + env_suffix
    params: 
        adata = 'results/per_sample/{sample}/adata.h5ad'
    shell:
        "papermill "
            "workflow/notebooks/feature_selection.ipynb "
            "results/per_sample/{wildcards.sample}/feature_selection.ipynb "
            "-p input_file {params.adata} "
            "-p count_layer {wildcards.count_layer} && "
        "touch {output}"


rule dimensionality_reduction:
    input: 'results/per_sample/{sample}/checkpoints/feature_selection_w_{count_layer}.done'
    output: 'results/per_sample/{sample}/checkpoints/dimensionality_reduction_w_{count_layer}.done'
    benchmark: 'benchmarks/dimensionality_reduction/{sample}_w_{count_layer}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "preprocessing" + env_suffix
    params: 
        adata = 'results/per_sample/{sample}/adata.h5ad'
    shell:
        "papermill "
            "workflow/notebooks/dimensionality_reduction.ipynb "
            "results/per_sample/{wildcards.sample}/dimensionality_reduction_w_{wildcards.count_layer}.ipynb "
            "-p input_file {params.adata} "
            "-p count_layer {wildcards.count_layer} && "
        "touch {output}"


rule clustering_per_sample:
    input: 
        checkpoint = 'results/per_sample/{sample}/checkpoints/dimensionality_reduction_w_{count_layer}.done'
    output: 'results/per_sample/{sample}/checkpoints/clustering_w_{count_layer}.done'
    benchmark: 'benchmarks/clustering/{sample}_w_{count_layer}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "preprocessing" + env_suffix
    params: 
        adata = 'results/preprocessing/{sample}/adata.h5ad'
    shell:
        "papermill "
            "workflow/notebooks/clustering.ipynb "
            "results/per_sample/{wildcards.sample}/clustering_w_{wildcards.count_layer}.ipynb "
            "-p input_file {params.adata} "
            "-p count_layer {wildcards.count_layer} && "
        "touch {output}"


rule annotate_per_sample:
    input: 
        checkpoint = 'results/per_sample/{sample}/checkpoints/clustering_w_{count_layer}.done'
    output: 'results/per_sample/{sample}/checkpoints/annotation_w_{count_layer}.done'
    benchmark: 'benchmarks/annotation/{sample}_w_{count_layer}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "preprocessing" + env_suffix
    params: 
        adata = 'results/per_sample/{sample}/adata.h5ad'
    shell:
        "papermill "
            "workflow/notebooks/annotation.ipynb "
            "results/per_sample/{wildcards.sample}/annotation_w_{wildcards.count_layer}.ipynb "
            "-p input_file {params.adata} "
            "-p count_layer {wildcards.count_layer} && "
        "touch {output}"


###-------------------- experimental, different quality control methods ------------------------###

rule scAutoQC:
    input: lambda wildcards: samples.at[wildcards.sample, 'path']
    output: 'results/per_sample/{sample}/scAutoQC.h5ad'
        # Excluded because the pipeline should not choke on this side product:
        # nb = 'results/per_sample/{sample}/scAutoQC.ipynb'
    benchmark: 'benchmarks/scAutoQC/{sample}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "sctk" + env_suffix
    shell:
        "papermill "
            "workflow/notebooks/sctk_scAutoQC.ipynb "
            "results/per_sample/{wildcards.sample}/scAutoQC.ipynb "
            "-p input_file {input} "
            "-p output_dir results/per_sample/{wildcards.sample}/"