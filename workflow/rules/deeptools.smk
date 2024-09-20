rule multiBigWigSummary:
    input:
        expand(f"results/bigwig/single/{bowtie2_dir}/{{sample}}.bw", sample=SAMPLES)
    output:
        f"results/deeptools/{bowtie2_dir}/bw_summary.npz",
    params:
        extra="",
    log:
        f"logs/deeptools/{bowtie2_dir}/multiBigWigSummary.log",
    threads: config["resources"]["deeptools"]["cpu"]
    resources: 
        runtime=config["resources"]["deeptools"]["time"]
    conda:
        "../envs/deeptools.yaml"
    shell:
        "multiBigwigSummary bins "
        "-p {threads} "
        "-b {input} "
        "-o {output} "
        "--smartLabels "
        "{params.extra}"


rule PCA:
    input:
        f"results/deeptools/{bowtie2_dir}/bw_summary.npz"
    output:
        f"results/deeptools/{bowtie2_dir}/pca.tab",
    params:
        extra="",
    log:
        f"logs/deeptools/{bowtie2_dir}/pca.log",
    threads: config["resources"]["deeptools"]["cpu"]
    resources: 
        runtime=config["resources"]["deeptools"]["time"]
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotPCA "
        "--corData {input} "
        "--outFileNameData {output} "
        "--transpose "
        "{params.extra} "
        "> {log} 2>&1"

#do this in for loop: value of mode can also be both scale-regions and reference-point
rule computeMatrix:
    input:
        bigwig=expand(f"results/bigwig/average_bw/{bowtie2_dir}/{{condition}}.bw", condition=CONDITIONS),
        bed=resources.gtf,
    output:
        matrix_gz=f"results/deeptools/{bowtie2_dir}/bw_matrix.gz",
    log:
        f"logs/deeptools/{bowtie2_dir}/computeMatrix.log",
    params:
        command=config["deeptools"]["computeMatrix"]["mode"],
        extra=computematrix_args(),
    threads: config["resources"]["deeptools"]["cpu"] * 4
    resources: 
        runtime=config["resources"]["deeptools"]["time"] * 2
    wrapper:
        f"{wrapper_version}/bio/deeptools/computematrix"
'''
#also a computeMatrix rule with regions being the peaks (for each condition)
#use mode center
rule computeMatrix_peaks:
    input:
        bed=expand(f"results/{peak_mode}/{bowtie2_dir}/fdr{fdr}/consensus_peaks/{{condition}}.bed", condition=CONDITIONS),
        bw=expand(f"results/bigwig/average_bw/{bowtie2_dir}/{{condition}}.bw", condition=CONDITIONS),
    output:
        matrix_gz=f"results/deeptools/{bowtie2_dir}/bw_matrix_peaks.gz",
    log:
        f"logs/deeptools/{bowtie2_dir}/computeMatrix_peaks.log",
    params:
''' 
