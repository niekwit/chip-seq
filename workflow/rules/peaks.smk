if config["peak_calling"]["macs2"]["run"]:
    if regions == "narrow":
        fdr = config["peak_calling"]["macs2"]["qvalue"]
        logger.info(f"MAC2S narrow peak calling selected with qvalue {fdr}")
        rule peak_calling_narrow:
            input:
                treatment=f"results/{bowtie2_dir}/{{ip_sample}}.dedup.bam",
                ipx=f"results/{bowtie2_dir}/{{ip_sample}}.dedup.bam.bai",
                control=f"results/{bowtie2_dir}/{{control_sample}}.dedup.bam",
                inptx=f"results/{bowtie2_dir}/{{control_sample}}.dedup.bam.bai",
            output:
                multiext(f"results/macs2_narrow/{bowtie2_dir}/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}", 
                "_peaks.narrowPeak", 
                "_peaks.xls",
                "_summits.bed",),
            params:
                macs2_params(),
            threads: config["resources"]["macs2"]["cpu"]
            resources: 
                runtime=config["resources"]["macs2"]["time"]
            log:
                f"logs/macs2_narrow/{bowtie2_dir}/fdr{fdr}/{{ip_sample}}_vs_{{control_sample}}.log"
            wrapper:
                f"{wrapper_version}/bio/macs2/callpeak"

        
        rule consensus_peaks:
            input:
                beds=expand(f"results/macs2_narrow/{bowtie2_dir}/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_peaks.narrowPeak", zip, ip_sample=IP_SAMPLES, control_sample=CONTROL_SAMPLES),
                chrom_sizes=f"resources/{resources.genome}_chrom.sizes",
            output:
                bed_out_intermediate=expand(f"results/macs2_narrow/{bowtie2_dir}/fdr{fdr}/consensus_peaks/{{condition}}.multi_intersect.bed", condition=CONDITIONS),
                bed_out=expand(f"results/macs2_narrow/{bowtie2_dir}/fdr{fdr}/consensus_peaks/{{condition}}.bed", condition=CONDITIONS),
            params:
                max_size=config["peak_calling"]["macs2"]["consensus_peaks"]["max_size"],
                extend_by=config["peak_calling"]["macs2"]["consensus_peaks"]["extend_by"],
                keep=config["peak_calling"]["macs2"]["consensus_peaks"]["keep"],
                conditions=CONDITIONS,
                extra=""
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                expand(f"logs/macs2_broad/{bowtie2_dir}/fdr{fdr}/consensus_peaks/{{condition}}.log", condition=CONDITIONS)
            conda:
                "../envs/deeptools.yaml"
            script:
                "../scripts/consensus_peaks.py"


        rule peak_annotation_plots:
            input:
                bed=expand(f"results/macs2_narrow/{bowtie2_dir}/fdr{fdr}/consensus_peaks/{{condition}}.bed", zip,  condition=CONDITIONS),
                gtf=resources.gtf,
            output:
                dt=f"results/plots/macs2_narrow/{bowtie2_dir}/fdr{fdr}/peaks_distance_to_TSS.pdf",
                fd=f"results/plots/macs2_narrow/{bowtie2_dir}/fdr{fdr}/peak_distributions.pdf",
            log: f"logs/plots/macs2_narrow/{bowtie2_dir}/fdr{fdr}/peak_annotation_plots.log"
            threads: config["resources"]["r"]["cpu"]
            resources:
                runtime=config["resources"]["r"]["time"]
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/peak_annotation_plots.R"
    
    elif regions == "broad":
        fdr = config["peak_calling"]["macs2"]["broad_cutoff"]
        logger.info(f"MAC2S broad peak calling selected with qvalue {fdr}")
        
        rule peak_calling_broad:
            input:
                treatment=f"results/{bowtie2_dir}/{{ip_sample}}.dedup.bam",
                ipx=f"results/{bowtie2_dir}/{{ip_sample}}.dedup.bam.bai",
                control=f"results/{bowtie2_dir}/{{control_sample}}.dedup.bam",
                inptx=f"results/{bowtie2_dir}/{{control_sample}}.dedup.bam.bai",
            output:
                multiext(f"results/macs2_broad/{bowtie2_dir}/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}", 
                "_peaks.broadPeak", 
                "_peaks.xls",),
            params:
                macs2_params(),
            threads: config["resources"]["macs2"]["cpu"]
            resources: 
                runtime=config["resources"]["macs2"]["time"]
            log:
                f"logs/macs2_broad/{bowtie2_dir}/fdr{fdr}/{{ip_sample}}_vs_{{control_sample}}.log"
            wrapper:
                f"{wrapper_version}/bio/macs2/callpeak"

        
        rule consensus_peaks:
            input:
                beds=expand(f"results/macs2_broad/{bowtie2_dir}/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_peaks.broadPeak", zip, ip_sample=IP_SAMPLES, control_sample=CONTROL_SAMPLES),
                chrom_sizes=f"resources/{resources.genome}_chrom.sizes",
            output:
                bed_out_intermediate=expand(f"results/macs2_broad/{bowtie2_dir}/fdr{fdr}/consensus_peaks/{{condition}}.multi_intersect.bed", condition=CONDITIONS),
                bed_out=expand(f"results/macs2_broad/{bowtie2_dir}/fdr{fdr}/consensus_peaks/{{condition}}.bed", condition=CONDITIONS),
            params:
                max_size=config["peak_calling"]["macs2"]["consensus_peaks"]["max_size"],
                extend_by=config["peak_calling"]["macs2"]["consensus_peaks"]["extend_by"],
                keep=config["peak_calling"]["macs2"]["consensus_peaks"]["keep"],
                conditions=CONDITIONS,
                extra=""
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                expand(f"logs/macs2_broad/{bowtie2_dir}/fdr{fdr}/consensus_peaks/{{condition}}.log", condition=CONDITIONS)
            conda:
                "../envs/deeptools.yaml"
            script:
                "../scripts/consensus_peaks.py"


        rule peak_annotation_plots:
            input:
                bed=expand(f"results/macs2_broad/{bowtie2_dir}/fdr{fdr}/consensus_peaks/{{condition}}.bed", condition=CONDITIONS),
                gtf=resources.gtf,
            output:
                dt=f"results/plots/macs2_broad/{bowtie2_dir}/fdr{fdr}/peaks_distance_to_TSS.pdf",
                fd=f"results/plots/macs2_broad/{bowtie2_dir}/fdr{fdr}/peak_distributions.pdf",
            log: f"logs/plots/macs2_broad/{bowtie2_dir}/fdr{fdr}/peak_annotation_plots.log"
            threads: config["resources"]["r"]["cpu"]
            resources:
                runtime=config["resources"]["r"]["time"]
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/peak_annotation_plots.R"

if config["peak_calling"]["htseq_deseq2"]["run"]:
    logger.info("Peak calling with htseq-count/DESeq2 selected")
    rule call_peaks_htseq_count:
        input:
            bam=f"results/{bowtie2_dir}/{{sample}}.dedup.bam",
            bai=f"results/{bowtie2_dir}/{{sample}}.dedup.bam.bai",
            gtf=resources.gtf,
        output:
            counts=f"results/htseq_count/{bowtie2_dir}/{{sample}}.tsv",
        params:
            mode=config["peak_calling"]["htseq_count"]["mode"],
            f=config["peak_calling"]["htseq_count"]["feature"],
            mapq=config["bowtie2"]["MAPQ_cutoff"],
            extra=config["peak_calling"]["htseq_count"]["extra"],
        threads: config["resources"]["deeptools"]["cpu"] * 2
        resources:
            runtime=config["resources"]["deeptools"]["time"]
        log:
            f"logs/peaks/htseq_count/{bowtie2_dir}/{{sample}}.log"
        conda:
            "../envs/macs2.yaml"
        shell:
            "htseq-count "
            "-m {params.mode} "
            "-f bam " # data is bam format
            "-r pos " # bam is sorted on position
            "-s no " # ChIP data is not stranded
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
            counts=expand(f"results/htseq_count/{bowtie2_dir}/{{sample}}.tsv", sample=SAMPLES),
        output:
            xlsx=f"results/htseq_count/DESeq2/{bowtie2_dir}/differential_peaks.xlsx",
            rdata=f"results/htseq_count/DESeq2/{bowtie2_dir}/dds.RData",
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
            f"logs/peaks/DESeq2/{bowtie2_dir}/differential_peaks.log"
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
            bam=expand(f"results/{bowtie2_dir}/{{ip_sample}}.dedup.bam", ip_sample=IP_SAMPLES),
            bai=expand(f"results/{bowtie2_dir}/{{ip_sample}}.dedup.bam.bai", ip_sample=IP_SAMPLES),
            bed=f"resources/{resources.genome}_windows_{window_size}kbp.bed",
        output:
            counts=f"results/multibamsummary_bed/{bowtie2_dir}/counts.tab_{window_size}kbp",
            npz=temp(f"results/multibamsummary_bed/{bowtie2_dir}/summary.npz"),
        params:
            labels=lambda wildcards, input: os.path.basename(input.bam).replace(".dedup.bam", ""),
            extra="",
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"]
        log:
            f"logs/multibamsummary_bed/{bowtie2_dir}/{sample}_{window_size}kbp.log"
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
            counts=f"results/multibamsummary_bed/{bowtie2_dir}/counts.tab_{window_size}kbp",
        output:
            pdf=f"results/plots/{bowtie2_dir}/RPKM_genomic_windows_{window_size}kbp.pdf",
