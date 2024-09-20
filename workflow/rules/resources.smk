rule get_fasta:
    output:
        resources.fasta,
    retries: 3
    params:
        url=resources.fasta_url,
    log:
        "logs/resources/get_fasta.log"
    cache: False
    conda:
        "../envs/trim_galore.yml"
    resources: 
        runtime=15
    script:
        "../scripts/get_resource.sh"


rule index_fasta:
    input:
        resources.fasta,
    output:
        f"{resources.fasta}.fai",
    log:
        "logs/resources/index_fasta.log"
    cache: False
    threads: config["resources"]["r"]["cpu"]
    resources: 
        runtime=config["resources"]["r"]["time"]
    wrapper:
        f"{wrapper_version}/bio/samtools/faidx"


rule chrom_sizes:
    input:
        fa=resources.fasta,
        fai=f"{resources.fasta}.fai",
    output:
        f"resources/{resources.genome}_chrom.sizes",
    log:
        "logs/resources/chrom_sizes.log"
    threads: config["resources"]["r"]["cpu"]
    resources: 
        runtime=config["resources"]["r"]["time"]
    conda:
        "../envs/deeptools.yaml"
    shell:
        "awk '{{print $1,$2}}' {input.fai} | "
        r"sed 's/ /\t/'  > {output}"


use rule get_fasta as get_black_list with:
    output:
        resources.blacklist,
    params:
        url=resources.blacklist_url,
    log:
        "logs/resources/get_black_list.log"


use rule get_fasta as get_gtf with:
        output:
            resources.gtf,
        params:
            url=resources.gtf_url,
        log:
            "logs/resources/get_gtf.log"


rule create_annotation_file:
    input:
        gtf=resources.gtf,
    output:
        rdata=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
    log:
        "logs/resources/create_annotation_file.log"
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/create_annotation_file.R"


rule bowtie2_build:
    input:
        ref=fasta
    output:
        multiext(
            f"resources/{resources.genome}_{resources.build}_index/index",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    params:
        extra=""
    cache: False
    log:
        "logs/bowtie2_build/index.log"
    threads: config["resources"]["mapping"]["cpu"] * 2
    resources: 
        runtime=config["resources"]["mapping"]["time"]
    wrapper:
        f"{wrapper_version}/bio/bowtie2/build"


if config["spike_in"]["apply"]:
    rule get_si_fasta:
        output:
            si_resources.fasta,
        retries: 3
        params:
            url=si_resources.fasta_url,
        log:
            "logs/resources/get_si_fasta.log"
        conda:
            "../envs/trim_galore.yml"
        cache: False
        resources: 
            runtime=15
        script:
            "../scripts/get_resource.sh"
    
    
    if not config["spike_in"]["downsample_only"]:
        rule bowtie2_build_spike_in:
            input:
                ref=si_resources.fasta
            output:
                multiext(
                    f"resources/{si_resources.genome}_{si_resources.build}_index/index",
                    ".1.bt2",
                    ".2.bt2",
                    ".3.bt2",
                    ".4.bt2",
                    ".rev.1.bt2",
                    ".rev.2.bt2",
                ),
            params:
                extra=""
            cache: False
            log:
                "logs/bowtie2_build_spike_in/index.log"
            threads: config["resources"]["mapping"]["cpu"] * 2
            resources: 
                runtime=config["resources"]["mapping"]["time"]
            wrapper:
                f"{wrapper_version}/bio/bowtie2/build"
