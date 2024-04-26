rule bowtie2_align:
    input:
        sample=["results/trimmed/{sample}_R1.fq.gz",
                "results/trimmed/{sample}_R2.fq.gz"],
        idx=multiext(
            f"resources/{resources.genome}_{resources.build}_index/index",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        temp("results/mapped/{sample}.bam"),
    log:
        "logs/bowtie2_align/{sample}.log",
    params:
        extra="",
    threads: config["resources"]["mapping"]["cpu"]
    resources: 
        runtime=config["resources"]["mapping"]["time"]
    wrapper:
        f"{wrapper_version}/bio/bowtie2/align"


rule bam_filtering:
    input:
        bam="results/mapped/{sample}.bam",
    output:
        temp("results/mapped/{sample}.filtered.bam"),
    params:
        extra=f"--min-MQ {config['samtools']['mapq']} -f 3",
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    log:
        "logs/mapq_filtering/{sample}.log"
    wrapper:
        f"{wrapper_version}/bio/samtools/view"


rule remove_blacklisted_regions:
    input:
        left="results/mapped/{sample}.filtered.bam",
        right=resources.blacklist,
    output:
        temp("results/mapped/{sample}.bl.bam"),
    params:
        extra="-v -nonamecheck",
    threads: config["resources"]["deeptools"]["cpu"]
    resources: 
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/remove_blacklisted_regions/{sample}.log"
    wrapper:
        f"{wrapper_version}/bio/bedtools/intersect"


rule sort:
    input:
        "results/mapped/{sample}.bl.bam",
    output:
        "results/mapped/{sample}.bl.sorted.bam",
    params:
        extra="-m 2G",
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    log:
        "logs/sort_bam/{sample}.log"
    wrapper:
        f"{wrapper_version}/bio/samtools/sort"


rule deduplication:
    input:
        bams="results/mapped/{sample}.bl.sorted.bam",
    output:
        bam="results/mapped/{sample}.dedup.bam",
        metrics="logs/deduplication/metrics.{sample}.txt",
    params:
        extra="--REMOVE_DUPLICATES true",
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    log:
        "logs/deduplication/{sample}.log"
    wrapper:
        f"{wrapper_version}/bio/picard/markduplicates"


rule bam_index:
    input:
        "results/mapped/{sample}.dedup.bam"
    output:
        "results/mapped/{sample}.dedup.bam.bai",
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    log:
        "logs/index_bam/{sample}.log"
    wrapper:
        f"{wrapper_version}/bio/samtools/index"