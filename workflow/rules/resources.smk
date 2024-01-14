rule get_fasta:
    output:
        fasta,
    retries: 3
    params:
        url=resources.fasta_url,
    log:
        "logs/resources/get_fasta.log"
    conda:
        "../envs/read-processing.yaml"
    shell:
        "wget -q {params.url} -O {output}.gz && gunzip -f {output}.gz 2> {log}"


rule get_black_list:
    output:
        resources.blacklist,
    retries: 3
    params:
        url=resources.blacklist_url,
    log:
        "logs/resources/get_black_list.log"
    conda:
        "../envs/read-processing.yaml"
    shell:
        "wget -q {params.url} -O {output}.gz && gunzip -f {output}.gz 2> {log}"



rule install_homer_genome:
    output:
        output=touch(f"resources/homer_{genome}_installed"),
    conda:
        "../envs/macs2.yaml"
    params:
        genome=genome,
    log:
        f"logs/homer/homer_install_{genome}.log"
    conda:
        "../envs/macs2.yaml"
    shell:
        "perl $CONDA_PREFIX/share/homer/configureHomer.pl -install human-p; "
        "perl $CONDA_PREFIX/share/homer/configureHomer.pl -install {params.genome} 2> {log}"