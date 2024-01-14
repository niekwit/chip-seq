if regions == "narrow":
    rule peak_calling_narrow_single:
        input:
            ip="results/mapped/{sample_ip}_dedup.bam",
            ipx="results/mapped/{sample_ip}_dedup.bam.bai",
            inpt="results/mapped/{sample_input}_dedup.bam",
            inptx="results/mapped/{sample_input}_dedup.bam.bai",
        output:
            multiext("results/peaks/narrow/{sample_ip}_vs_{sample_input}/{sample_ip}", 
            "_peaks.narrowPeak", 
            "_peaks.xls",
            "_summits.bed",),
        params:
            f=config["macs2"]["format"],
            g=re.sub(r"\d+$","",genome),
            q=config["macs2"]["qvalue"],
            extra=config["macs2"]["extra_params"],
            d="results/peaks/narrow/",
        threads: config["resources"]["macs2"]["cpu"]
        resources: 
            runtime=config["resources"]["macs2"]["time"]
        log:
            "logs/macs2/{sample_ip}_vs_{sample_input}.log"
        conda:
            "../envs/macs2.yaml"
        shell:
            "macs2 callpeak -t {input.ip} -c {input.inpt} -n {wildcards.sample_ip} --outdir {params.d}{wildcards.sample_ip}_vs_{wildcards.sample_input} -f {params.f} -g {params.g} -q {params.q} {params.extra} 2> {log}"

   
    rule peak_calling_narrow_all:
        input:
            ip=expand("results/mapped/{sample_ip}_dedup.bam", sample_ip=IP_SAMPLES),
            ipx=expand("results/mapped/{sample_ip}_dedup.bam.bai", sample_ip=IP_SAMPLES),
            inpt=expand("results/mapped/{sample_input}_dedup.bam", sample_input=INPUT_SAMPLES),
            inptx=expand("results/mapped/{sample_input}_dedup.bam.bai", sample_input=INPUT_SAMPLES),
        output:
            expand("results/peaks/narrow/{name}/{name}_peaks.narrowPeak", name=UNIQUE_IP_NAMES),
            expand("results/peaks/narrow/{name}/{name}_peaks.xls", name=UNIQUE_IP_NAMES),
            expand("results/peaks/narrow/{name}/{name}_summits.bed", name=UNIQUE_IP_NAMES),
        params:
            f=config["macs2"]["format"],
            g=genome,
            q=config["macs2"]["qvalue"],
            e=config["macs2"]["extra_params"],
            r=regions,
            d="results/peaks/narrow/",
        threads: config["resources"]["macs2"]["cpu"]
        resources: 
            runtime=config["resources"]["macs2"]["time"]
        conda:
            "../envs/macs2.yaml"
        log:
            expand("logs/macs2/{name}/{name}.log", name=UNIQUE_IP_NAMES)
        script:
            "../scripts/macs2_replicates.py"

elif regions == "broad":
    regions_param = "--broad --broadcutoff"

    rule peak_calling_broad_single:
        input:
            ip="results/mapped/{sample_ip}_dedup.bam",
            inpt="results/mapped/{sample_input}_dedup.bam",
        output:
            multiext("results/peaks/{regions}/{sample_ip}_vs_{sample_input}/{sample_ip}", 
            "_peaks.broadPeak", 
            "_peaks.xls",
            "_summits.bed",
            "_peaks.gappedPeak",),
        params:
            g=config["macs2"]["genome"],
            q=config["macs2"]["qvalue"],
            r=regions_param,
            rc=config["macs2"]["broad_cutoff"],
            extra=config["macs2"]["extra_params"],
            d="results/peaks/broad/",
        threads: config["resources"]["macs2"]["cpu"]
        resources: 
            runtime=config["resources"]["macs2"]["time"]
        log:
            "logs/macs2/{sample_ip}_vs_{sample_input}.log"
        conda:
            "../envs/macs2.yaml"
        shell:
            "macs2 callpeak -t {input.ip} -c {input.inpt} -n {wildcards.sample_ip} --outdir {params.d}{wildcards.sample_ip} "
            "-f BAMPE -g {params.g} -q {params.q} {params.r} {params.rc} {params.extra} 2> {log}"


rule convert_xls2bed:
    input:
        "results/peaks/{regions}/{sample_ip}_vs_{sample_input}/{sample_ip}_peaks.xls"
    output:
        "results/peaks/{regions}/{sample_ip}_vs_{sample_input}/{sample_ip}_peaks.bed"
    shell: #this command is split over several lines to avoid issues with quotes 
        r"tail -n +30 {input} | awk -v OFS='\t' '{{print $1,$2,$3,$10," #raw string escapes \t
        '".","+"}}'
        "' "
        "| sed -e 's/^/chr/' > {output}"


rule annotate_peaks:
        input:
            bed="results/peaks/{regions}/{sample_ip}_vs_{sample_input}/{sample_ip}_peaks.bed",
            f=f"resources/homer_{genome}_installed",
        output:
            "results/peaks/{regions}/{sample_ip}_vs_{sample_input}/{sample_ip}_peaks_annotated.txt"
        params:
            genome=config["genome"],
        threads: config["resources"]["macs2"]["cpu"]
        resources: 
            runtime=config["resources"]["macs2"]["time"]
        conda:
            "../envs/macs2.yaml"
        log: 
            "logs/homer/{regions}/{sample_ip}_vs_{sample_input}_annotatePeaks.log"
        shell:
            "annotatePeaks.pl {input.bed} {params.genome} > {output} 2> {log}"


rule convert_xls2bed_merged:
    input:
        "results/peaks/{regions}/{name}/{name}_peaks.xls"
    output:
        "results/peaks/{regions}/{name}/{name}_peaks.bed"
    log:
        "logs/xls2bed/{regions}/{name}.log"
    shell:
        r"tail -n +30 {input} | awk -v OFS='\t' '{{print $1,$2,$3,$10,"
        '".","+"}}'
        "' "
        "| sed -e 's/^/chr/' > {output} 2> {log}"


rule annotate_peaks_merged:
    input:
        bed="results/peaks/{regions}/{name}/{name}_peaks.bed",
        f=f"resources/homer_{genome}_installed",
    output:
        "results/peaks/{regions}/{name}/{name}_peaks_annotated.txt"
    params:
        genome=config["genome"],
    threads: config["resources"]["macs2"]["cpu"]
    resources: 
        runtime=config["resources"]["macs2"]["time"]
    conda:
        "../envs/macs2.yaml"
    log:
        "logs/homer/{regions}/{name}_annotatePeaks.log"
    shell:
        "annotatePeaks.pl {input.bed} {params.genome} > {output} 2> {log}"