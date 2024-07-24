rule fastqc_post_trim:
    input:
        "results/trimmed/{sample}_R{end}.fq.gz",
    output:
        html="results/qc/fastqc/post-trim/{sample}_R{end}.html",
        zip="results/qc/fastqc/post-trim/{sample}_R{end}_fastqc.zip"
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/post-trim/{sample}{end}.log"
    threads: config["resources"]["fastqc"]["cpu"]
    resources:
        mem_mb = 2048
    wrapper:
        f"{wrapper_version}/bio/fastqc"


rule fastqc_pre_trim:
    input:
        "reads/{sample}_R{end}_001.fastq.gz",
    output:
        html="results/qc/fastqc/pre-trim/{sample}_R{end}.html",
        zip="results/qc/fastqc/pre-trim/{sample}_R{end}_fastqc.zip"
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/pre-trim/{sample}{end}.log"
    threads: config["resources"]["fastqc"]["cpu"]
    resources:
        mem_mb = 2048
    wrapper:
        f"{wrapper_version}/bio/fastqc"


rule multiqc:
    input:
        expand("results/qc/fastqc/{dir}/{sample}_R{end}_fastqc.zip", sample=SAMPLES, end=["1","2"], dir=["pre-trim", "post-trim"]),
    output:
        "results/qc/multiqc.html",
        directory("results/qc/multiqc_data"),
    params:
        extra="--dirs --dirs-depth 1",  # Optional: extra parameters for multiqc
    log:
        "logs/multiqc.log"
    wrapper:
        f"{wrapper_version}/bio/multiqc"


rule readcounts: # Read count pre and post duplication
    input:
        pre=expand(f"results/{bowtie2_dir}/{{sample}}.bl.sorted.bam", sample=SAMPLES),
        post=expand(f"results/{bowtie2_dir}/{{sample}}.dedup.bam", sample=SAMPLES),
    output:
        f"results/qc/{bowtie2_dir}/read_counts.csv",
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    conda:
        "../envs/plotting.yaml"
    log:
        f"logs/qc/{bowtie2_dir}/readcounts.log"
    script:
        "../scripts/readcounts.py"
