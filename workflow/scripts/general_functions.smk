def targets():
    # Base targets
    TARGETS = [
        f"results/plots/{bowtie2_dir}/readcounts_plot.pdf",
        f"results/plots/{bowtie2_dir}/alignment-rates.pdf",
        f"results/plots/{bowtie2_dir}/PCA.pdf",
        f"results/plots/{bowtie2_dir}/scree.pdf",
        f"results/plots/{bowtie2_dir}/bam_fingerprint.pdf",
        f"results/plots/{bowtie2_dir}/sample_correlation.pdf",
        f"results/plots/{bowtie2_dir}/profile_genome.pdf",
        f"results/plots/{bowtie2_dir}/heatmap_genome.pdf",
        expand(f"results/bigwig/single/{bowtie2_dir}/{{sample}}.bw", sample=SAMPLES),
        # "results/qc/multiqc.html",
    ]

    if config["peak_calling"]["macs2"]["run"]:
        TARGETS.extend(
            [
                f"results/plots/{bowtie2_dir}/profile_peaks.pdf",
                f"results/plots/{bowtie2_dir}/heatmap_peaks.pdf",
            ]
        )
        if regions == "narrow":
            TARGETS.extend(
                [
                    expand(
                        f"results/macs2_narrow/{bowtie2_dir}/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_peaks.narrowPeak",
                        zip,
                        ip_sample=IP_SAMPLES,
                        control_sample=CONTROL_SAMPLES,
                    ),
                    expand(
                        f"results/macs2_narrow/{bowtie2_dir}/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_peaks.xls",
                        zip,
                        ip_sample=IP_SAMPLES,
                        control_sample=CONTROL_SAMPLES,
                    ),
                    expand(
                        f"results/macs2_narrow/{bowtie2_dir}/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_summits.bed",
                        zip,
                        ip_sample=IP_SAMPLES,
                        control_sample=CONTROL_SAMPLES,
                    ),
                    f"results/plots/macs2_narrow/{bowtie2_dir}/fdr{fdr}/peaks_distance_to_TSS.pdf",
                    f"results/plots/macs2_narrow/{bowtie2_dir}/fdr{fdr}/peak_distributions.pdf",
                ]
            )
        elif regions == "broad":
            TARGETS.extend(
                [
                    expand(
                        f"results/macs2_broad/{bowtie2_dir}/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_peaks.broadPeak",
                        zip,
                        ip_sample=IP_SAMPLES,
                        control_sample=CONTROL_SAMPLES,
                    ),
                    expand(
                        f"results/macs2_broad/{bowtie2_dir}/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_peaks.xls",
                        zip,
                        ip_sample=IP_SAMPLES,
                        control_sample=CONTROL_SAMPLES,
                    ),
                    f"results/plots/macs2_broad/{bowtie2_dir}/fdr{fdr}/peaks_distance_to_TSS.pdf",
                    f"results/plots/macs2_broad/{bowtie2_dir}/fdr{fdr}/peak_distributions.pdf",
                ]
            )
    return TARGETS


def samples(paired_end):
    """
    Imports all samples from config/samples.csv
    """
    csv = pd.read_csv("config/samples.csv")
    SAMPLES = csv["sample"].tolist()
    SAMPLES.extend(csv["control"].unique().tolist())

    # Check if samples from samples.csv match fastq files (both R1 and R2) in reads folder
    not_found = []
    if paired_end:
        for sample in SAMPLES:
            r1 = f"reads/{sample}_R1_001.fastq.gz"
            r2 = f"reads/{sample}_R2_001.fastq.gz"

            if not os.path.isfile(r1):
                not_found.append(r1)
            if not os.path.isfile(r2):
                not_found.append(r2)

            if len(not_found) > 0:
                not_found = "\n".join(not_found)
                raise FileNotFoundError(
                    f"Following fastq files not found:\n{not_found}"
                )
    else:
        for sample in SAMPLES:
            r = f"reads/{sample}.fastq.gz"

            if not os.path.isfile(r):
                not_found.append(r)

            if len(not_found) > 0:
                not_found = "\n".join(not_found)
                raise FileNotFoundError(
                    f"Following fastq files not found:\n{not_found}"
                )

    assert len(SAMPLES) > 0, "No samples found in config/samples.csv"

    return SAMPLES


def ip_samples():
    """
    Get all IP samples from config/samples.csv
    """
    csv = pd.read_csv("config/samples.csv")
    return csv["sample"].tolist()


def control_samples():
    """
    Get all input samples from config/samples.csv
    """
    csv = pd.read_csv("config/samples.csv")
    return csv["control"].tolist()


def conditions():
    """
    Get all unique IP sample names from config/samples.csv
    """
    # Remove underscore and trailing numbers from IP sample names
    return list(set([re.sub(r"_*[0-9]+$", "", x) for x in IP_SAMPLES]))


def effective_genome_size():
    """
    Returns genome argument for bamCoverage based on genome specified in config file.
    """
    read_length = str(config["read_length"])

    if resources.genome == "hg19":
        genome = "GRCh37"
    elif resources.genome == "hg38":
        genome = "GRCh38"
    elif resources.genome == "mm38":
        genome = "GRCm38"
    elif resources.genome == "mm39":
        genome = "GRCm39"
    elif resources.genome == "dm3":
        genome = "dm3"
    elif resources.genome == "dm6":
        genome = "dm6"
    else:
        raise ValueError("Unsupported genome specified in config file...")

    default_effective_genome_size = {
        "dm3": {
            "50": 130428560,
            "75": 135004462,
            "100": 139647232,
            "150": 144307808,
            "200": 148524010,
        },
        "dm6": {
            "50": 125464728,
            "75": 127324632,
            "100": 129789873,
            "150": 129941135,
            "200": 132509163,
        },
        "GRCh37": {
            "50": 2685511504,
            "75": 2736124973,
            "100": 2776919808,
            "150": 2827437033,
            "200": 2855464000,
        },
        "GRCh38": {
            "50": 2701495761,
            "75": 2747877777,
            "100": 2805636331,
            "150": 2862010578,
            "200": 2887553303,
        },
        "GRCm37": {
            "50": 2304947926,
            "75": 2404646224,
            "100": 2462481010,
            "150": 2489384235,
            "200": 2513019276,
        },
        "GRCm38": {
            "50": 2308125349,
            "75": 2407883318,
            "100": 2467481108,
            "150": 2494787188,
            "200": 2520869189,
        },
        "GRCm39": {  # copy from GRCm38
            "50": 2308125349,
            "75": 2407883318,
            "100": 2467481108,
            "150": 2494787188,
            "200": 2520869189,
        },
    }
    return default_effective_genome_size[genome][read_length]


def computematrix_genome_args():
    """
    Returns computeMatrix arguments as string based on config file.
    """
    # Add mode argument
    mode = config["deeptools"]["computeMatrix"]["mode"]
    if mode == "scale-regions":
        rbl = config["deeptools"]["computeMatrix"]["regionBodyLength"]
        args = f"--regionBodyLength {rbl} "
    else:
        rp = config["deeptools"]["computeMatrix"]["referencePoint"]
        args = f"--referencePoint {rp} "

    # Add common arguments
    b = config["deeptools"]["computeMatrix"]["upstream"]
    a = config["deeptools"]["computeMatrix"]["downstream"]
    bs = config["deeptools"]["computeMatrix"]["binSize"]
    atb = config["deeptools"]["computeMatrix"]["averageTypeBins"]

    args = f"{args} --upstream {b} --downstream {a} --binSize {bs} --averageTypeBins {atb} --missingDataAsZero "

    # Add region argument
    r = f"--regionsFileName {resources.gtf}"
    args = f"{args} {r}"

    return args


def computematrix_peaks_args():
    """
    Returns computeMatrix arguments as string based on config file.
    """
    # Add mode argument
    args = "--referencePoint center"

    # Add common arguments
    b = config["deeptools"]["computeMatrix"]["upstream"]
    a = config["deeptools"]["computeMatrix"]["downstream"]
    bs = config["deeptools"]["computeMatrix"]["binSize"]
    atb = config["deeptools"]["computeMatrix"]["averageTypeBins"]

    args = f"{args} --upstream {b} --downstream {a} --binSize {bs} --averageTypeBins {atb} --missingDataAsZero --smartLabels "

    return args


def macs2_params(paired_end):
    """
    Returns MACS2 parameters based on the genome and mode
    (narrow or broad) specified in the config file.
    """
    if paired_end:
        format_ = "BAMPE"
    else:
        format_ = "BAM"

    if "hg" in resources.genome:
        genome = "hs"
    elif "mm" in resources.genome:
        genome = "mm"
    elif "dm" in resources.genome:
        genome = "dm"

    if config["peak_calling"]["macs2"]["regions"] == "broad":
        cutoff = config["peak_calling"]["macs2"]["broad_cutoff"]
        broad = f"--broad --broad-cutoff {cutoff} "
        qvalue = ""
    else:
        broad = ""
        qvalue = config["peak_calling"]["macs2"]["qvalue"]
        qvalue = f"-q {qvalue}"

    extra = config["peak_calling"]["macs2"]["extra"]

    return f"-f {format_} -g {genome} {qvalue} {broad} {extra}"


def peak_mode():
    """
    Returns MACS2 peak calling mode as string based on config file
    """
    if config["peak_calling"]["macs2"]["regions"] == "broad":
        return "macs2_broad"
    else:
        return "macs2_narrow"


def bowtie2_dir():
    if config["bowtie2"]["k_mode"] == 0:
        return f"mapped_q{config["samtools"]["mapq"]}"
    else:
        return f"mapped_k{config['bowtie2']['k_mode']}"


def paired_end():
    """
    Checks if paired-end reads are used
    """
    # Get one fastq file
    reads = glob.glob("reads/*fastq.gz")
    if len(reads) == 0:
        reads = glob.glob("reads/*fastq.gz")
    assert len(reads) != 0, "No fastq files found..."

    fastq = reads[0]

    # Check file extension to see if paired-end reads are used
    if fastq.endswith("_R1_001.fastq.gz"):
        logger.info("Paired-end reads detected...")
        return True
    elif fastq.endswith("_R2_001.fastq.gz"):
        logger.info("Paired-end reads detected...")
        return True
    else:
        logger.info("Single-end reads detected...")
        return False


def bamcoverage_args(paired_end):
    """
    Returns bamCoverage arguments as string based on config file.
    """
    extend_reads = (
        f"--extendReads {config['deeptools']['bigwig']['extendReads']}"
        if config["deeptools"]["bigwig"]["extendReads"]
        else ""
    )
    bin_size = f"-bs {config['deeptools']['bigwig']['binSize']}"
    normalize_using = (
        f"--normalizeUsing {config['deeptools']['bigwig']['normalizeUsing']}"
    )
    extra = config["deeptools"]["bigwig"]["extra"]

    return f"{extend_reads} {bin_size} {normalize_using} {extra}"
