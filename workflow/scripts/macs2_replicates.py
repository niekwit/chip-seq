import pandas as pd
import re
from snakemake.shell import shell

# write all stdout and stderr to log file
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# get snakemake params
format= snakemake.params["f"]
genome = snakemake.params["g"]
qvalue = snakemake.params["q"]
extra = snakemake.params["e"]
d = snakemake.params["d"]

# get sample info
csv = pd.read_csv("config/samples.csv")

# set MACS2 genome size parameter (removes trailing numbers from genome name)
genome = re.sub(r"\d+$","",genome)

# create MACS2 command for unique IP sample name with matching input files
for n in snakemake.wildcards["name"]:
    # get replicate IP BAM files
    ip = csv[csv["Sample"].str.contains(n)]["Sample"].tolist()
    ip_bams = " ".join([f"results/mapped/{x}_dedup.bam" for x in ip])
    
    # get replicate matching input samples
    inpt = csv[csv["Sample"].str.contains(n)]["Input"].tolist()
    input_bams = " ".join([f"results/mapped/{x}_dedup.bam" for x in inpt])
    
    # create output dir string
    outdir = f"{d}/{n}"
    
    # run MACS2 for replicate IP and matching input samples
    shell(
        "macs2 callpeak "
        f"-t {ip_bams} "
        f"-c {input_bams} "
        f"-n {n} "
        f"--outdir {outdir} "
        f"-f {format} "
        f"-g {genome} "
        f"-q {qvalue} "
        f"{extra} {log}"
        )


