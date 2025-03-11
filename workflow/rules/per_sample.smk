rule quality_control:
    input: lambda wildcards: samples.at[wildcards.sample, 'path']
    output: 
        adata = 'results/per_sample/{sample}/adata.h5ad',
    benchmark: 'benchmarks/quality_control/{sample}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 13 * attempt ** 2, # 1h on debug partition, then 4h on short partition, then medium
    conda: env_prefix + "preprocessing" + env_suffix
    shell:
        "papermill "
            "workflow/notebooks/quality_control.ipynb "
            "results/per_sample/{wildcards.sample}/quality_control.ipynb "
            "-p input_file {input} "
            "-p output_dir results/per_sample/{wildcards.sample}/ "


rule per_sample_analysis:
    input: lambda wildcards: 'results/per_sample/{sample}/adata.h5ad' if wildcards.qc_method == "theislab_tutorial" else 'results/per_sample/{sample}/adata_scAutoQC_{count_layer}.h5ad'
    output: 'results/per_sample/{sample}/adata_ready_for_merge_{count_layer}__{qc_method}.h5ad'
    benchmark: 'benchmarks/per_sample_analysis/{sample}_{count_layer}__{qc_method}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 81 * attempt ** 2,
    conda: env_prefix + "preprocessing" + env_suffix
    shell:
        "papermill "
            "workflow/notebooks/per_sample_analysis.ipynb "
            "results/per_sample/{wildcards.sample}/per_sample_analysis_{wildcards.count_layer}_{wildcards.qc_method}.ipynb "
            "-p input_file {input} "
            "-p output_file {output} "
            "-p qc_method {wildcards.qc_method} "
            "-p count_layer {wildcards.count_layer}"


###-------------------- experimental, different quality control methods ------------------------###

rule scAutoQC:
    input: 'results/per_sample/{sample}/adata.h5ad'
    output: 'results/per_sample/{sample}/adata_scAutoQC_{count_layer}.h5ad'
    benchmark: 'benchmarks/scAutoQC/{sample}_{count_layer}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (4 * attempt),
        runtime=lambda wildcards, attempt: 6 * attempt ** 2,
    conda: env_prefix + "sctk" + env_suffix
    shell:
        "papermill "
            "workflow/notebooks/sctk_scAutoQC.ipynb "
            "results/per_sample/{wildcards.sample}/scAutoQC.ipynb "
            "-p input_file {input} "
            "-p count_layer {wildcards.count_layer} "
            "-p output_file {output} "
