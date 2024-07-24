rule bigwig:
    input:
        bam=f"results/{bowtie2_dir}/{{sample}}.dedup.bam",
        bai=f"results/{bowtie2_dir}/{{sample}}.dedup.bam.bai",
    output:
        f"results/bigwig/single/{bowtie2_dir}/{{sample}}.bw",
    params:
        effective_genome_size=effective_genome_size(),
        extra=f"--extendReads -bs {config['deeptools']['bigwig']['binSize']} --normalizeUsing {config['deeptools']['bigwig']['normalizeUsing']} {config['deeptools']['bigwig']['extra']}",
    log:
        f"logs/deeptools/bigwig/{bowtie2_dir}/{{sample}}.log",
    threads: config["resources"]["deeptools"]["cpu"]
    resources: 
        runtime=config["resources"]["deeptools"]["time"]
    wrapper:
        f"{wrapper_version}/bio/deeptools/bamcoverage"


rule average_wig:
    input:
        expand(f"results/bigwig/single/{bowtie2_dir}/{{ip_sample}}.bw", ip_sample=IP_SAMPLES),
    output:
        wig=temp(f"results/bigwig/average_bw/{bowtie2_dir}/{{condition}}.wig"),
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/wiggletools/wig_average_{condition}.log"
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/average_wig.py"


rule wig2bigwig:
    input:
        wig=f"results/bigwig/average_bw/{bowtie2_dir}/{{condition}}.wig",
        cs=f"resources/{resources.genome}_chrom.sizes",
    output:
        f"results/bigwig/average_bw/{bowtie2_dir}/{{condition}}.bw",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        f"logs/wigToBigWig/{bowtie2_dir}/{{condition}}.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "wigToBigWig {input.wig} {input.cs} {output}"
