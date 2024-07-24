rule plot_readcounts: # plot read counts pre and post deduplication
    input:
        f"results/qc/{bowtie2_dir}/read_counts.csv",
    output:
        report(f"results/plots/{bowtie2_dir}/readcounts_plot.pdf", caption="workflow/report/readcounts_plot.rst", category="Read counts"),
    threads: 1
    resources: 
        runtime=10
    conda:
        "../envs/R.yaml"
    log:
        f"logs/qc/{bowtie2_dir}/readcount_plot.log"
    script:
        "../scripts/readcount_plot.R"


rule plot_alignment_rates:
    input:
        expand(f"logs/bowtie2_align/{bowtie2_dir}/{{sample}}.log", sample=SAMPLES)
    output:
        report(f"results/plots/{bowtie2_dir}/alignment-rates.pdf", caption="workflow/report/alignment-rates.rst", category="Alignment rates")
    log:
        f"logs/qc/{bowtie2_dir}/plot-alignment-rate.log"
    threads: 1
    resources:
        runtime=10
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_alignment_rate.R"


rule plot_PCA:
    input:
        f"results/deeptools/{bowtie2_dir}/pca.tab",
    output:
        pca=report(f"results/plots/{bowtie2_dir}/PCA.pdf", caption="../report/pca.rst", category="PCA"),
        scree=report(f"results/plots/{bowtie2_dir}/scree.pdf", caption="../report/scree.rst", category="PCA"),
    params:
        extra=""
    threads: 1
    resources:
        runtime=10
    log:
        f"logs/plotting/{bowtie2_dir}/plotPCA.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_PCA.R"


rule plot_bam_fingerprint:
    input:
        bam_files=expand(f"results/{bowtie2_dir}/{{sample}}.dedup.bam", sample=SAMPLES),
        bam_idx=expand(f"results/{bowtie2_dir}/{{sample}}.dedup.bam.bai", sample=SAMPLES),
    output:
        fingerprint=report(f"results/plots/{bowtie2_dir}/bam_fingerprint.pdf", caption="../report/bam_fingerprint.rst", category="BAM Fingerprint"),
        qc_metrics=f"logs/deeptools/{bowtie2_dir}/plot_fingerprint_qc_metrics.txt",
    log:
        f"logs/deeptools/{bowtie2_dir}/bam_fingerprint.log"
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    wrapper:
        f"{wrapper_version}/bio/deeptools/plotfingerprint"


rule plot_correlation:
    input:
        f"results/deeptools/{bowtie2_dir}/bw_summary.npz",
    output:
        tab=f"results/deeptools/{bowtie2_dir}/correlation.tab",
        pdf=report(f"results/plots/{bowtie2_dir}/sample_correlation.pdf", caption="../report/correlation.rst", category="Sample correlation"),
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        f"logs/plotting/{bowtie2_dir}/plotCorrelation.log"
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
        f"results/deeptools/{bowtie2_dir}/bw_matrix.gz"
    output:
        plot_img=report(f"results/plots/{bowtie2_dir}/profile.pdf", caption="../report/profile.rst", category="Profile"),
        regions=f"results/deeptools/{bowtie2_dir}/plot_profile_regions.bed",
    log:
        f"logs/deeptools/{bowtie2_dir}/plot_profile.log"
    params:
        "--plotType=lines "
        "--perGroup "
        "--legendLocation=upper-right  "
    wrapper:
        f"{wrapper_version}/bio/deeptools/plotprofile"  


rule plot_heatmap:
    input:
        f"results/deeptools/{bowtie2_dir}/bw_matrix.gz"
    output:
        heatmap_img=report(f"results/plots/{bowtie2_dir}/heatmap.pdf", caption="../report/heatmap.rst", category="Heatmap"),
    log:
        f"logs/deeptools/{bowtie2_dir}/plot_heatmap.log"
    params:
        "--colorMap viridis "
        "--whatToShow 'plot, heatmap and colorbar' "
    wrapper:
        f"{wrapper_version}/bio/deeptools/plotheatmap"