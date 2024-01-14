rule trim_galore_pe:
    input:
        r1="reads/{sample}_R1.fastq.gz", 
        r2="reads/{sample}_R2.fastq.gz",
    output:
        r1=temp("results/trimmed/{sample}_val_1.fq.gz"),
        r2=temp("results/trimmed/{sample}_val_2.fq.gz"),
    threads: config["resources"]["trim"]["cpu"],
    params:
        extra="--illumina -q 20",
    log:
        "logs/trim_galore/{sample}.log",
    conda:
        "../envs/read-processing.yaml",
    shell:
        "trim_galore -j {threads} "
        "--basename {wildcards.sample} "
        "-o results/trimmed "
        "--paired {input.r1} {input.r2} "
        "{params.extra} " 
        " 2> {log}"
