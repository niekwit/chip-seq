rule fastqc:
    input:
        "reads/{sample}{end}.fastq.gz"
    output:
        html="results/qc/fastqc/{sample}{end}.html",
        zip="results/qc/fastqc/{sample}{end}_fastqc.zip"
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}{end}.log"
    threads: config["resources"]["fastqc"]["cpu"]
    wrapper:
        "v2.0.0/bio/fastqc"


rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}{end}_fastqc.zip", sample=SAMPLES, end=["_R1","_R2"])
    output:
        "results/qc/multiqc.html",
    params:
        extra="",  # Optional: extra parameters for multiqc
    log:
        "logs/multiqc/multiqc.log"
    conda:
        "../envs/read-processing.yaml"
    shell:
        "multiqc " 
        "--force "
        "-o results/qc "
        "-n multiqc.html "
        "{params.extra} "
        "{input} "
        "2> {log}"


rule readcounts: # read count pre and post duplication
    input:
        pre=expand("results/mapped/{sample}_bl-sorted.bam", sample=SAMPLES),
        post=expand("results/mapped/{sample}_dedup.bam", sample=SAMPLES),
    output:
        "results/qc/readcounts.csv",
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    conda:
        "../envs/plotting.yaml"
    log:
        "logs/qc/readcounts.log"
    script:
        "../scripts/readcounts.py"


rule readcount_plot: # plot read counts pre and post deduplication
    input:
        "results/qc/readcounts.csv",
    output:
        "results/qc/readcounts_plot.pdf",
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    conda:
        "../envs/plotting.yaml"
    log:
        "logs/qc/readcount_plot.log"
    script:
        "../scripts/readcount_plot.R"


rule plot_alignment_rate:
    input:
        expand("logs/hisat2_align/{sample}.log", sample=SAMPLES)
    output:
        report("results/qc/alignment-rates.pdf", caption="workflow/report/alignment-rates.rst", category="Alignment rates")
    log:
        "logs/qc/plot-alignment-rate.log"
    conda:
        "../envs/plotting.yaml"
    script:
        "../scripts/plot_alignment_rate.R"