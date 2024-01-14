rule hisat2_index:
    input:
        fasta=fasta
    output:
        dir=directory("resources/index"),
    params:
        prefix=f"resources/index/{genome}"
    log:
        "logs/hisat2_index/index.log"
    threads: config["resources"]["mapping"]["cpu"]
    resources: 
        runtime=config["resources"]["mapping"]["time"]
    wrapper:
        "v2.10.0/bio/hisat2/index"


rule hisat2_align:
    input:
        r1="results/trimmed/{sample}_val_1.fq.gz",
        r2="results/trimmed/{sample}_val_2.fq.gz",
        idx="resources/index",
    output:
        temp("results/mapped/{sample}.bam"),
    log:
        "logs/hisat2_align/{sample}.log",
    params:
        extra="",
        prefix=f"resources/index/{genome}"
    threads: config["resources"]["mapping"]["cpu"]
    resources: 
        runtime=config["resources"]["mapping"]["time"]
    conda:
        "../envs/read-processing.yaml"
    shell:
        "hisat2 -t -p {threads} -x {params.prefix} -1 {input.r1} -2 {input.r2} 2> {log} | "
        "samtools view -F 260 -bS -@ {threads} > {output}"


rule remove_blacklisted_regions:
    input:
        bam="results/mapped/{sample}.bam",
        bl=resources.blacklist,
    output:
        temp("results/mapped/{sample}_bl.bam"),
    threads: config["resources"]["deeptools"]["cpu"]
    resources: 
        runtime=config["resources"]["deeptools"]["time"]
    conda:
        "../envs/read-processing.yaml"
    log:
        "logs/remove_blacklisted_regions/{sample}.log"
    shell:
        "bedtools intersect -v -a {input.bam} -b {input.bl} -nonamecheck > {output} 2> {log}"


rule sort:
    input:
        "results/mapped/{sample}_bl.bam",
    output:
        "results/mapped/{sample}_bl-sorted.bam",
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    conda:
        "../envs/read-processing.yaml"
    log:
        "logs/sort/{sample}.log"
    shell:
        "samtools sort -@ {threads} -o {output} {input} 2> {log}"


rule deduplication:
    input:
        "results/mapped/{sample}_bl-sorted.bam"
    output:
        "results/mapped/{sample}_dedup.bam"
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    conda:
        "../envs/read-processing.yaml"
    log:
        "logs/deduplication/{sample}.log"
    shell:
        "picard MarkDuplicates INPUT={input} OUTPUT={output} REMOVE_DUPLICATES=TRUE "
        "METRICS_FILE=logs/deduplication/picard.{wildcards.sample}.log 2> {log}"


rule bam_index:
    input:
        "results/mapped/{sample}_dedup.bam"
    output:
        "results/mapped/{sample}_dedup.bam.bai",
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    conda:
        "../envs/read-processing.yaml"
    log:
        "logs/index/{sample}.log"
    shell:
        "samtools index -@ {threads} {input} 2> {log}"


