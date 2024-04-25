rule plot_readcounts: # plot read counts pre and post deduplication
    input:
        "results/qc/read_counts.csv",
    output:
        report("results/plots/readcounts_plot.pdf", caption="workflow/report/readcounts_plot.rst", category="Read counts"),
    threads: 1
    resources: 
        runtime=10
    conda:
        "../envs/R.yaml"
    log:
        "logs/qc/readcount_plot.log"
    script:
        "../scripts/readcount_plot.R"


rule plot_alignment_rates:
    input:
        expand("logs/bowtie2_align/{sample}.log", sample=SAMPLES)
    output:
        report("results/plots/alignment-rates.pdf", caption="workflow/report/alignment-rates.rst", category="Alignment rates")
    log:
        "logs/qc/plot-alignment-rate.log"
    threads: 1
    resources:
        runtime=10
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_alignment_rate.R"


rule plot_PCA:
    input:
        "results/deeptools/pca.tab",
    output:
        pca=report("results/plots/PCA.pdf", caption="../report/pca.rst", category="PCA"),
        scree=report("results/plots/scree.pdf", caption="../report/scree.rst", category="PCA"),
    params:
        extra=""
    threads: 1
    resources:
        runtime=10
    log:
        "logs/plotting/plotPCA.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_PCA.R"


rule plot_bam_fingerprint:
    input:
        bam_files=expand("results/mapped/{sample}.dedup.bam", sample=SAMPLES),
        bam_idx=expand("results/mapped/{sample}.dedup.bam.bai", sample=SAMPLES),
    output:
        fingerprint=report("results/plots/bam_fingerprint.pdf", caption="../report/bam_fingerprint.rst", category="BAM Fingerprint"),
        qc_metrics="logs/deeptools/plot_fingerprint_qc_metrics.txt",
    log:
        "logs/deeptools/bam_fingerprint.log"
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    wrapper:
        f"{wrapper_version}/bio/deeptools/plotfingerprint"


rule plot_correlation:
    input:
        "results/deeptools/bw_summary.npz",
    output:
        tab="results/deeptools/correlation.tab",
        pdf=report("results/plots/sample_correlation.pdf", caption="../report/correlation.rst", category="Sample correlation"),
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/plotting/plotCorrelation.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotCorrelation "
        "--corData {input} "
        "--corMethod spearman "
        "--whatToPlot heatmap "
        "--removeOutliers "
        "--plotFile {output.pdf} "
        "--outFileCorMatrix {output.tab} "
        "--colorMap viridis "
        "--skipZeros "
        "> {log} 2>&1"


rule plot_profile:
    input:
        "results/deeptools/bw_matrix.gz"
    output:
        plot_img=report("results/plots/profile.pdf", caption="../report/profile.rst", category="Profile"),
        regions="results/deeptools/plot_profile_regions.bed",
    log:
        "logs/deeptools/plot_profile.log"
    params:
        "--plotType=lines "
        "--perGroup "
        "--legendLocation=upper-right  "
    wrapper:
        f"{wrapper_version}/bio/deeptools/plotprofile"  


rule plot_heatmap:
    input:
        "results/deeptools/bw_matrix.gz"
    output:
        heatmap_img=report("results/plots/heatmap.pdf", caption="../report/heatmap.rst", category="Heatmap"),
    log:
        "logs/deeptools/plot_heatmap.log"
    params:
        "--colorMap viridis "
        "--whatToShow 'plot, heatmap and colorbar' "
    wrapper:
        f"{wrapper_version}/bio/deeptools/plotheatmap"