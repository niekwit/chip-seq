rule multiBigWigSummary:
    input:
        expand("results/bigwig/single/{sample}.bw", sample=SAMPLES)
    output:
        "results/deeptools/bw_summary.npz",
    params:
        extra="",
    log:
        "logs/deeptools/multiBigWigSummary.log",
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
        "results/deeptools/bw_summary.npz"
    output:
        "results/deeptools/pca.tab",
    params:
        extra="",
    log:
        "logs/deeptools/pca.log",
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


rule computeMatrix:
    input:
        bigwig=expand("results/bigwig/average_bw/{condition}.bw", condition=CONDITIONS),
        bed=resources.gtf,
    output:
        matrix_gz="results/deeptools/bw_matrix.gz",
    log:
        "logs/deeptools/computeMatrix.log",
    params:
        command=config["deeptools"]["computeMatrix"]["mode"],
        extra=computematrix_args(),
    threads: config["resources"]["deeptools"]["cpu"] * 4
    resources: 
        runtime=config["resources"]["deeptools"]["time"] * 2
    wrapper:
        f"{wrapper_version}/bio/deeptools/computematrix"