import re
import os
import pandas as pd
from snakemake.shell import shell

# Get snakemake params
genome = snakemake.params["genome"]
qvalue = snakemake.params["qvalue"]
extra = snakemake.params["extra"]
peak_files = snakemake.output["xls"]
logs = snakemake.log
base_dir = snakemake.params["base_dir"]

csv = pd.read_csv("config/samples.csv")

# Set MACS2 genome size parameter 
if "hg" in genome:
    genome = "hs"
else:
    # Remove trailing numbers from genome name
    genome = re.sub(r"\d+$", "", genome)

# Get conditions from peak xls files
conditions = [os.path.basename(x).replace("_peaks.xls", "") for x in peak_files]

# Create MACS2 command for unique IP sample name with matching input files
for condition in conditions:
    # Get log file
    log = [x for x in logs if condition in x]
    
    # Get replicate IP BAM files
    ip = csv[csv["sample"].str.contains(condition)]["sample"].tolist()
    ip_bams = " ".join([f"results/mapped/{x}.dedup.bam" for x in ip])
    
    # Get replicate matching input samples
    # Some samples may have the same control samples
    inpt = csv[csv["sample"].str.contains(condition)]["control"].unique().tolist()
    input_bams = " ".join([f"results/mapped/{x}.dedup.bam" for x in inpt])
    
    # Create output dir string
    outdir = f"{base_dir}/{condition}"
    
    # Run MACS2 for replicate IP and matching input samples
    shell(
        "macs2 callpeak "
        "-t {ip_bams} "
        "-c {input_bams} "
        "-n {condition} "
        "--outdir {outdir} "
        "-f BAMPE "
        "-g {genome} "
        "-q {qvalue} "
        "{extra} "
        "> {log} 2>&1"
        )
