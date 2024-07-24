rule fastqc_post_trim:
    input:
        "results/trimmed/{sample}_R{end}.fq.gz",
    output:
        html="results/qc/fastqc_trimmed/{sample}_R{end}.html",
        zip="results/qc/fastqc_trimmed/{sample}_R{end}_fastqc.zip"
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}{end}.log"
    threads: config["resources"]["fastqc"]["cpu"]
    resources:
        mem_mb = 2048
    wrapper:
        f"{wrapper_version}/bio/fastqc"


rule fastqc_pre_trim:
    input:
        "reads/{sample}_R{end}_001.fastq.gz",
    output:
        html="results/qc/fastqc/{sample}_R{end}.html",
        zip="results/qc/fastqc/{sample}_R{end}_fastqc.zip"
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}{end}.log"
    threads: config["resources"]["fastqc"]["cpu"]
    resources:
        mem_mb = 2048
    wrapper:
        f"{wrapper_version}/bio/fastqc"


rule multiqc:
    input:
        expand("results/qc/{dir}/{sample}_R{end}_fastqc.zip", sample=SAMPLES, end=["1","2"], dir=["fastqc", "fastqc_trimmed"]),
    output:
        "results/qc/multiqc.html",
        directory("results/qc/multiqc_data"),
    params:
        extra="",  # Optional: extra parameters for multiqc
    log:
        "logs/multiqc.log"
    wrapper:
        f"{wrapper_version}/bio/multiqc"


rule readcounts: # Read count pre and post duplication
    input:
        pre=expand("results/mapped/{sample}.bl.sorted.bam", sample=SAMPLES),
        post=expand("results/mapped/{sample}.dedup.bam", sample=SAMPLES),
    output:
        "results/qc/read_counts.csv",
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    conda:
        "../envs/plotting.yaml"
    log:
        "logs/qc/readcounts.log"
    script:
        "../scripts/readcounts.py"
