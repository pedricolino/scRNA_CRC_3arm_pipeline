rule quality_control:
    input: lambda wildcards: samples.at[wildcards.sample, 'path']
    output: 'results/preprocessing/{sample}/quality_control.h5ad'
        # Excluded because the pipeline should not choke on this side product:
        # nb = 'results/preprocessing/{sample}/quality_control.ipynb'
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
            "-p output_dir results/preprocessing/{wildcards.sample}/"


rule normalization:
    input: 'results/preprocessing/{sample}/quality_control.h5ad'
    output: 'results/preprocessing/{sample}/normalization.h5ad'
    benchmark: 'benchmarks/normalization/{sample}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "preprocessing" + env_suffix
    shell:
        "papermill "
            "workflow/notebooks/normalization.ipynb "
            "results/preprocessing/{wildcards.sample}/normalization.ipynb "
            "-p input_file {input} "
            "-p output_dir results/preprocessing/{wildcards.sample}/"


rule feature_selection:
    input: 'results/preprocessing/{sample}/normalization.h5ad'
    output: 'results/preprocessing/{sample}/feature_selection.h5ad'
    benchmark: 'benchmarks/feature_selection/{sample}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda: env_prefix + "preprocessing" + env_suffix
    shell:
        "papermill "
            "workflow/notebooks/feature_selection.ipynb "
            "results/preprocessing/{wildcards.sample}/feature_selection.ipynb "
            "-p input_file {input} "
            "-p output_dir results/preprocessing/{wildcards.sample}/"


rule dimensionality_reduction:
    input: 'results/preprocessing/{sample}/feature_selection.h5ad'
    output: 'results/preprocessing/{sample}/dimensionality_reduction.h5ad'
    benchmark: 'benchmarks/dimensionality_reduction/{sample}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 1*60 if attempt == 1 else 4*60,
    conda:
        env_prefix + "preprocessing" + env_suffix
    shell:
        "papermill "
            "workflow/notebooks/dimensionality_reduction.ipynb "
            "results/preprocessing/{wildcards.sample}/dimensionality_reduction.ipynb "
            "-p input_file {input} "
            "-p output_dir results/preprocessing/{wildcards.sample}/"