if config["peak_calling"]["macs2"]["run"]:
    if regions == "narrow":
        fdr = peak_fdr("macs2_narrow")
        logger.info(f"MAC2S narrow peak calling selected with qvalue {fdr}...")
        rule peak_calling_narrow_single:
            input:
                treatment="results/mapped/{ip_sample}.dedup.bam",
                ipx="results/mapped/{ip_sample}.dedup.bam.bai",
                control="results/mapped/{control_sample}.dedup.bam",
                inptx="results/mapped/{control_sample}.dedup.bam.bai",
            output:
                multiext(f"results/macs2_narrow/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}", 
                "_peaks.narrowPeak", 
                "_peaks.xls",
                "_summits.bed",),
            params:
                macs2_params(),
            threads: config["resources"]["macs2"]["cpu"]
            resources: 
                runtime=config["resources"]["macs2"]["time"]
            log:
                f"logs/macs2_narrow/fdr{fdr}/{{ip_sample}}_vs_{{control_sample}}.log"
            wrapper:
                    f"{wrapper_version}/bio/macs2/callpeak"

        
        rule macs2_narrow_replicates:
                input:
                    bams=expand("results/mapped/{ip_sample}.dedup.bam", ip_sample=IP_SAMPLES),
                    bais=expand("results/mapped/{ip_sample}.dedup.bam.bai", ip_sample=IP_SAMPLES),
                    cbams=expand("results/mapped/{control_sample}.dedup.bam", control_sample=CONTROL_SAMPLES), 
                    cbais=expand("results/mapped/{control_sample}.dedup.bam.bai", control_sample=CONTROL_SAMPLES),
                output:
                    xls=expand(f"results/macs2_narrow/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.xls", conditions=CONDITIONS),
                    peak=expand(f"results/macs2_narrow/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.narrowPeak", conditions=CONDITIONS),
                    bed=expand(f"results/macs2_narrow/fdr{fdr}/{{conditions}}/{{conditions}}_summits.bed", conditions=CONDITIONS),
                params:
                    base_dir=lambda wildcards, output: os.path.dirname(os.path.dirname(output["xls"][0])),
                    genome=resources.genome,
                    mode="narrow",
                    qvalue=config["peak_calling"]["macs2"]["qvalue"],
                    extra=config["peak_calling"]["macs2"]["extra"],
                threads: config["resources"]["macs2"]["cpu"]
                resources:
                    runtime=config["resources"]["macs2"]["time"] * 4
                conda:
                    "../envs/macs2.yaml",
                log:
                    expand(f"logs/macs2_narrow/fdr{fdr}/{{conditions}}.log", conditions=CONDITIONS)
                script:
                    "../scripts/macs2_replicates.py"

        rule peak_annotation_plots:
            input:
                bed=expand(f"results/macs2_narrow/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.narrowPeak", conditions=CONDITIONS),
                gtf=resources.gtf,
            output:
                dt=f"results/plots/macs2_narrow/fdr{fdr}/peaks_distance_to_TSS.pdf",
                fd=f"results/plots/macs2_narrow/fdr{fdr}/peak_distributions.pdf",
            log: f"logs/plots/macs2_narrow/fdr{fdr}/peak_annotation_plots.log"
            threads: config["resources"]["r"]["cpu"]
            resources:
                runtime=config["resources"]["r"]["time"]
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/peak_annotation_plots.R"
    
    if regions == "broad":
        fdr = peak_fdr("macs2_broad")
        logger.info(f"MAC2S broad peak calling selected with qvalue {fdr}...")
        rule peak_calling_broad_single:
            input:
                treatment="results/mapped/{ip_sample}.dedup.bam",
                ipx="results/mapped/{ip_sample}.dedup.bam.bai",
                control="results/mapped/{control_sample}.dedup.bam",
                inptx="results/mapped/{control_sample}.dedup.bam.bai",
            output:
                multiext(f"results/macs2_broad/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}", 
                "_peaks.broadPeak", 
                "_peaks.xls",
                "_gappedPeak",),
            params:
                macs2_params(),
            threads: config["resources"]["macs2"]["cpu"]
            resources: 
                runtime=config["resources"]["macs2"]["time"]
            log:
                f"logs/macs2_broad/fdr{fdr}/{{ip_sample}}_vs_{{control_sample}}.log"
            wrapper:
                    f"{wrapper_version}/bio/macs2/callpeak"
        

        rule macs2_broad_replicates:
            input:
                bams=expand("results/mapped/{ip_sample}.dedup.bam", ip_sample=IP_SAMPLES),
                bais=expand("results/mapped/{ip_sample}.dedup.bam.bai", ip_sample=IP_SAMPLES),
                cbams=expand("results/mapped/{control_sample}.dedup.bam", control_sample=CONTROL_SAMPLES), 
                cbais=expand("results/mapped/{control_sample}.dedup.bam.bai", control_sample=CONTROL_SAMPLES),
            output:
                xls=expand(f"results/macs2_broad/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.xls", conditions=CONDITIONS),
                peak=expand(f"results/macs2_broad/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.narrowPeak", conditions=CONDITIONS),
                bed=expand(f"results/macs2_broad/fdr{fdr}/{{conditions}}/{{conditions}}_summits.bed", conditions=CONDITIONS),
            params:
                base_dir=lambda wildcards, output: os.path.dirname(os.path.dirname(output["xls"][0])),
                genome=resources.genome,
                mode="broad",
                qvalue=config["peak_calling"]["macs2"]["qvalue"],
                extra=config["peak_calling"]["macs2"]["extra"],
            threads: config["resources"]["macs2"]["cpu"]
            resources:
                runtime=config["resources"]["macs2"]["time"] * 4
            conda:
                "../envs/macs2.yaml",
            log:
                expand(f"logs/macs2_broad/fdr{fdr}/{{conditions}}.log", conditions=CONDITIONS)
            script:
                "../scripts/macs2_replicates.py"
        

        rule peak_annotation_plots:
            input:
                bed=expand(f"results/macs2_broad/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.broadPeak", conditions=CONDITIONS),
                gtf=resources.gtf,
            output:
                dt=f"results/plots/macs2_broad/fdr{fdr}/peaks_distance_to_TSS.pdf",
                fd=f"results/plots/macs2_broad/fdr{fdr}/peak_distributions.pdf",
            log: f"logs/plots/macs2_broad/fdr{fdr}/peak_annotation_plots.log"
            threads: config["resources"]["r"]["cpu"]
            resources:
                runtime=config["resources"]["r"]["time"]
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/peak_annotation_plots.R"

    rule annotate_peaks:
        input:
            xls=f"results/{peak_mode}/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.xls",
            adb=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
            gtf=resources.gtf,
        output:
            bed=f"results/{peak_mode}/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.bed",
            txt=f"results/{peak_mode}/fdr{fdr}/{{conditions}}/{{conditions}}_annotated.peaks.txt",
        params:
            pm=peak_mode,
            extra=""
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"]
        log:
            f"logs/annotate_peaks/{peak_mode}/fdr{fdr}/{{conditions}}.log"
        conda:
            "../envs/R.yaml"
        script:
            "../scripts/annotate_peaks.R"


    rule get_gene_names_macs2:
            input:
                txt=f"results/{peak_mode}/fdr{fdr}/{{conditions}}/{{conditions}}_annotated.peaks.txt",
            output:
                ids=f"results/{peak_mode}/fdr{fdr}/{{conditions}}.geneIDs.txt",
            threads: 1
            resources:
                runtime=5
            log:
                f"logs/geneIDs_peaks/{peak_mode}/fdr{fdr}/{{conditions}}.log"
            conda:
                "../envs/deeptools.yaml"
            shell:
                "sed '1d' {input.txt} | "
                "awk '{{print $(NF-4),$(NF-1)}}' |"
                " sort | "
                "uniq > {output.ids}"

if config["peak_calling"]["htseq_deseq2"]["run"]:
    logger.info("Peak calling with htseq-count/DESeq2 selected...")
    rule call_peaks_htseq_count:
        input:
            bam="results/mapped/{sample}.dedup.bam",
            bai="results/mapped/{sample}.dedup.bam.bai",
            gtf=resources.gtf,
        output:
            counts="results/htseq_count/{sample}.tsv",
        params:
            mode=config["peak_calling"]["htseq_count"]["mode"],
            f=config["peak_calling"]["htseq_count"]["feature"],
            mapq=config["bowtie2"]["MAPQ_cutoff"],
            extra=config["peak_calling"]["htseq_count"]["extra"],
        threads: config["resources"]["deeptools"]["cpu"] * 2
        resources:
            runtime=config["resources"]["deeptools"]["time"]
        log:
            "logs/peaks/htseq_count/{sample}.log"
        conda:
            "../envs/macs2.yaml"
        shell:
            "htseq-count "
            "-m {params.mode} "
            "-f bam " # data is bam format
            "-r pos " # bam is sorted on position
            "-s no " # Cut & Run data is not stranded
            "-t {params.f} " # get signal over whole gene, instead of just exons
            "-i gene_id " # use gene_id as identifier
            "-a {params.mapq} " # MAPQ cutoff
            "--additional-attr=gene_name " # use for annotation
            "--additional-attr=gene_biotype " # use for annotation
            "-n {threads} "
            "{params.extra} "
            "{input.bam} "
            "{input.gtf} "
            "2> {log} | "
            r"sed 's/\t\t/\tNA\t/g' > {output.counts}" # replace empty fields with NA (genes with no gene_name attributes)

    rule differential_peaks_DESeq2:
        input:
            counts=expand("results/htseq_count/{sample}.tsv", sample=SAMPLES),
        output:
            xlsx="results/htseq_count/DESeq2/differential_peaks.xlsx",
            rdata="results/htseq_count/DESeq2/dds.RData",
        params:
            alpha=config["peak_calling"]["htseq_count"]["DESeq2"]["alpha"],
            fc=config["peak_calling"]["htseq_count"]["DESeq2"]["fc"],
            control=config["peak_calling"]["htseq_count"]["DESeq2"]["deseq2_apply_control"],
            cfo=config["peak_calling"]["htseq_count"]["DESeq2"]["cumulative_filter_out"],
            sg=config["peak_calling"]["htseq_count"]["DESeq2"]["smallest_group"],
            extra="",
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"]
        log:
            "logs/peaks/DESeq2/differential_peaks.log"
        conda:
            "../envs/R.yaml"
        script:
            "../scripts/differential_peaks_DESeq2.R"

if config["peak_calling"]["genomic_blanket"]["run"]:
    logger.info("...")
    window_size = config["peak_calling"]["genomic_blanket"]["window_size"]
    rule create_genomic_windows_bed:
        input:
            cs=f"resources/{resources.genome}_chrom.sizes",
        output:
            bed=f"resources/{resources.genome}_windows_{window_size}kbp.bed",
        params:
            size=window_size,
            extra="",
        log:
            "logs/create_genomic_windows.log"
        conda:
            "../envs/deeptools.yaml"
        script:
            "../scripts/create_genomic_windows.py"


    rule read_coverage_bed:
        input:
            bam=expand("results/mapped/{ip_sample}.dedup.bam", ip_sample=IP_SAMPLES),
            bai=expand("results/mapped/{ip_sample}.dedup.bam.bai", ip_sample=IP_SAMPLES),
            bed=f"resources/{resources.genome}_windows_{window_size}kbp.bed",
        output:
            counts=f"results/multibamsummary_bed/counts.tab_{window_size}kbp",
            npz=temp("results/multibamsummary_bed/summary.npz"),
        params:
            labels=lambda wildcards, input: os.path.basename(input.bam).replace(".dedup.bam", ""),
            extra="",
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"]
        log:
            f"logs/multibamsummary_bed/{sample}_{window_size}kbp.log"
        conda:
            "../envs/deeptools.yaml"
        shell:
            "multiBamSummary BED-file "
            "-p {threads} "
            "--BED {input.bed} "
            "--bamfiles {input.bam} "
            "--labels {params.labels} "
            "--outRawCounts {output.counts} "
            "-o {output.npz} "
            "{params.extra} "
            "> {log} 2>&1"


    rule RPKM_genomic_windows:
        input:
            counts=f"results/multibamsummary_bed/counts.tab_{window_size}kbp",
        output:
            pdf=f"results/plots/RPKM_genomic_windows_{window_size}kbp.pdf",
