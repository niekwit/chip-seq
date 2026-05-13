rule plot_readcounts:  # plot read counts pre and post deduplication
    input:
        f"results/qc/{bowtie2_dir}/read_counts.csv",
    output:
        report(
            f"results/plots/{bowtie2_dir}/readcounts_plot.pdf",
            caption="workflow/report/readcounts_plot.rst",
            category="Read counts",
        ),
    threads: 1
    resources:
        runtime=10,
    conda:
        "../envs/R.yaml"
    log:
        f"logs/qc/{bowtie2_dir}/readcount_plot.log",
    script:
        "../scripts/readcount_plot.R"


rule plot_alignment_rates:
    input:
        expand(f"logs/bowtie2_align/{bowtie2_dir}/{{sample}}.log", sample=SAMPLES),
    output:
        report(
            f"results/plots/{bowtie2_dir}/alignment-rates.pdf",
            caption="workflow/report/alignment-rates.rst",
            category="Alignment rates",
        ),
    log:
        f"logs/qc/{bowtie2_dir}/plot-alignment-rate.log",
    threads: 1
    resources:
        runtime=10,
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_alignment_rate.R"


rule plot_PCA:
    input:
        f"results/deeptools/{bowtie2_dir}/pca.tab",
    output:
        pca=report(
            f"results/plots/{bowtie2_dir}/PCA.pdf",
            caption="../report/pca.rst",
            category="PCA",
        ),
        scree=report(
            f"results/plots/{bowtie2_dir}/scree.pdf",
            caption="../report/scree.rst",
            category="PCA",
        ),
    params:
        extra="",
    threads: 1
    resources:
        runtime=10,
    log:
        f"logs/plotting/{bowtie2_dir}/plotPCA.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_PCA.R"


rule plot_correlation:
    input:
        f"results/deeptools/{bowtie2_dir}/bw_summary.npz",
    output:
        tab=f"results/deeptools/{bowtie2_dir}/correlation.tab",
        pdf=report(
            f"results/plots/{bowtie2_dir}/sample_correlation.pdf",
            caption="../report/correlation.rst",
            category="Sample correlation",
        ),
    params:
        extra="",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"],
    log:
        f"logs/plotting/{bowtie2_dir}/plotCorrelation.log",
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


rule plot_profile_genome:
    input:
        f"results/deeptools/{bowtie2_dir}/bw_matrix_genome.gz",
    output:
        plot_img=report(
            f"results/plots/{bowtie2_dir}/profile_genome.pdf",
            caption="../report/profile.rst",
            category="Profile",
        ),
        regions=f"results/deeptools/{bowtie2_dir}/plot_profile_regions_genome.bed",
    log:
        f"logs/deeptools/{bowtie2_dir}/plot_profile_genome.log",
    params:
        "--plotType=lines " "--perGroup " "--legendLocation=upper-right  ",
    wrapper:
        f"{wrapper_version}/bio/deeptools/plotprofile"


use rule plot_profile_genome as plot_profile_peaks with:
    input:
        f"results/deeptools/{bowtie2_dir}/bw_matrix_peaks.gz",
    output:
        plot_img=report(
            f"results/plots/{bowtie2_dir}/profile_peaks.pdf",
            caption="../report/profile.rst",
            category="Profile",
        ),
        regions=f"results/deeptools/{bowtie2_dir}/plot_profile_regions_peaks.bed",
    log:
        f"logs/deeptools/{bowtie2_dir}/plot_profile_peaks.log",


rule plot_heatmap_genome:
    input:
        f"results/deeptools/{bowtie2_dir}/bw_matrix_genome.gz",
    output:
        heatmap_img=report(
            f"results/plots/{bowtie2_dir}/heatmap_genome.pdf",
            caption="../report/heatmap.rst",
            category="Heatmap",
        ),
    log:
        f"logs/deeptools/{bowtie2_dir}/plot_heatmap_genome.log",
    params:
        "--colorMap viridis " "--whatToShow 'plot, heatmap and colorbar' ",
    wrapper:
        f"{wrapper_version}/bio/deeptools/plotheatmap"


use rule plot_heatmap_genome as plot_heatmap_peaks with:
    input:
        f"results/deeptools/{bowtie2_dir}/bw_matrix_peaks.gz",
    output:
        heatmap_img=report(
            f"results/plots/{bowtie2_dir}/heatmap_peaks.pdf",
            caption="../report/heatmap.rst",
            category="Heatmap",
        ),
    log:
        f"logs/deeptools/{bowtie2_dir}/plot_heatmap_peaks.log",
