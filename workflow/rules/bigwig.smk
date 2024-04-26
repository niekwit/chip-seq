rule bigwig:
    input:
        bam="results/mapped/{sample}.dedup.bam",
        bai="results/mapped/{sample}.dedup.bam.bai",
    output:
        "results/bigwig/single/{sample}.bw",
    params:
        effective_genome_size=effective_genome_size(),
        extra=f"--extendReads -bs {config['deeptools']['bigwig']['binSize']} --normalizeUsing {config['deeptools']['bigwig']['normalizeUsing']} {config['deeptools']['bigwig']['extra']}",
    log:
        "logs/deeptools/bigwig/{sample}.log",
    threads: config["resources"]["deeptools"]["cpu"]
    resources: 
        runtime=config["resources"]["deeptools"]["time"]
    wrapper:
        f"{wrapper_version}/bio/deeptools/bamcoverage"


rule average_wig:
    input:
        expand("results/bigwig/single/{ip_sample}.bw", ip_sample=IP_SAMPLES),
    output:
        wig=temp("results/bigwig/average_bw/{condition}.wig"),
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
        wig="results/bigwig/average_bw/{condition}.wig",
        cs=f"resources/{resources.genome}_chrom.sizes",
    output:
        "results/bigwig/average_bw/{condition}.bw",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/wigToBigWig/{condition}.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "wigToBigWig {input.wig} {input.cs} {output}"
