rule bigwig:
    input:
        bam="mapped/{sample}_dedup.bam",
        index="mapped/{sample}_dedup.bam.bai",
    output:
        "results/bigwig/single/{sample}.bw",
    params:
        bs=config["bigwig"]["binSize"],
    threads: config["resources"]["deeptools"]["cpu"]
    resources: 
        runtime=config["resources"]["deeptools"]["time"]
    conda:
        "envs/deeptools.yaml"
    shell:
        "bamCoverage -p {threads} -b {input.bam} -o {output} --extendReads"