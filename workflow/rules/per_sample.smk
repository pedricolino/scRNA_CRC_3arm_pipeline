rule quality_control:
    input: lambda wildcards: samples.at[wildcards.sample, 'path']
    output: 
        adata = 'results/per_sample/{sample}/adata.h5ad',
    benchmark: 'benchmarks/quality_control/{sample}.tsv'
    threads: 8
    params:
        use_ensembl_ids = config['use_ensembl_ids'],
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=lambda wildcards, attempt: 13 * attempt ** 2, # 1h on debug partition, then 4h on short partition, then medium
    conda: env_prefix + "preprocessing" + env_suffix
    shell:
        "papermill "
            "workflow/notebooks/quality_control.ipynb "
            "results/per_sample/{wildcards.sample}/quality_control.ipynb "
            "-p input_file {input} "
            "-p use_ensembl_ids {params.use_ensembl_ids} "
            "-p output_dir results/per_sample/{wildcards.sample}/ "


rule per_sample_analysis:
    input: lambda wildcards: 'results/per_sample/{sample}/adata.h5ad' if wildcards.qc_method == "theislab_tutorial" else 'results/per_sample/{sample}/adata_scAutoQC_{count_layer}.h5ad'
    output: 'results/per_sample/{sample}/adata_ready_for_merge_{count_layer}__{qc_method}.h5ad'
    wildcard_constraints:
         # Exclude patterns starting with "metacells_", otherwise this and the metacell computation rule might have identical output files
        count_layer="(?!metacells_).*"
    benchmark: 'benchmarks/per_sample_analysis/{sample}_{count_layer}__{qc_method}.tsv'
    threads: 8
    params:
        use_ensembl_ids = config['use_ensembl_ids'],
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
            "-p use_ensembl_ids {params.use_ensembl_ids} "
            "-p qc_method {wildcards.qc_method} "
            "-p count_layer {wildcards.count_layer}"


rule scAutoQC:
    input: 'results/per_sample/{sample}/adata.h5ad'
    output: 'results/per_sample/{sample}/adata_scAutoQC_{count_layer}.h5ad'
    benchmark: 'benchmarks/scAutoQC/{sample}_{count_layer}.tsv'
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (16 * attempt),
        runtime=lambda wildcards, attempt: 6 * attempt ** 2,
    conda: env_prefix + "sctk" + env_suffix
    shell:
        "papermill "
            "workflow/notebooks/sctk_scAutoQC.ipynb "
            "results/per_sample/{wildcards.sample}/scAutoQC.ipynb "
            "-p input_file {input} "
            "-p count_layer {wildcards.count_layer} "
            "-p output_file {output} "

if config['use_metacells']:
    rule SEACells_metacell_computation:
        input: 'results/per_sample/{sample}/adata_ready_for_merge_{count_layer}__{qc_method}.h5ad'
        output: 'results/per_sample/{sample}' + metacell_suffix + '/adata_ready_for_merge_{count_layer}__{qc_method}.h5ad'
        benchmark: 'benchmarks/SEACells_metacell_computation/{sample}_{count_layer}_{qc_method}.tsv'
        threads: 8
        resources:
            mem=lambda wildcards, attempt: '%dG' % (16 * attempt),
            runtime=lambda wildcards, attempt: 240 * attempt ** 2
        conda: env_prefix + "seacells" + env_suffix
        shell:
            "papermill "
            "workflow/notebooks/SEACell_computation.ipynb "
            "results/per_sample/{wildcards.sample}" + metacell_suffix + "/SEACell_computation_{wildcards.count_layer}_{wildcards.qc_method}.ipynb "
            "-p input_file {input} "
            "-p output_file {output} "
            "-p count_layer {wildcards.count_layer} "

rule velocyto:
    input: lambda wildcards: samples.at[wildcards.sample, 'path'].replace('sample_filtered_feature_bc_matrix.h5', 'sample_alignments.bam')
    output: 'results/per_sample/{sample}/velocyto/.done'
    conda: env_prefix + "velocity" + env_suffix
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (50 * attempt),
        runtime=lambda wildcards, attempt: 240 * attempt ** 2,
    shell:
        "velocyto run {input} resources/genes.gtf -o results/per_sample/{wildcards.sample}/velocyto && touch {output}"

rule numbat_preprocess:
    input:
        script = config['numbat']['script'],
        panel = config['numbat']['panel'],
        genetic_map = config['numbat']['genetic_map'],
        snp_vcf = config['numbat']['snp_vcf'],
        bam = lambda wildcards: samples.at[wildcards.sample, 'path'].replace('sample_filtered_feature_bc_matrix.h5', 'sample_alignments.bam'),
        barcodes = lambda wildcards: samples.at[wildcards.sample, 'path'].replace('sample_filtered_feature_bc_matrix.h5', 'sample_raw_feature_bc_matrix/barcodes.tsv.gz'),
    output: 'results/per_sample/{sample}/numbat/preprocessing.done'
    params:
        outdir = lambda wildcards: 'results/per_sample/{}/numbat'.format(wildcards.sample),
        label = 'exp',
        sample = lambda wildcards: wildcards.sample,
    conda: env_prefix + "R" + env_suffix
    threads: 32
    resources:
        mem=lambda wildcards, attempt: '%dG' % (50 * attempt),
        runtime=lambda wildcards, attempt: 60 * 24 * 4 * attempt ** 2,
    shell:
        "Rscript {input.script} "
        "--label {params.label} "
        "--samples {params.sample} "
        "--barcodes {input.barcodes} "
        "--bams {input.bam} "
        "--outdir {params.outdir} "
        "--gmap {input.genetic_map} "
        "--snpvcf {input.snp_vcf} "
        "--paneldir {input.panel} "
        "--ncores {threads} && "
        "touch {output}"
