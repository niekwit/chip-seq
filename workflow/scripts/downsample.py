"""
Downsample all BAM files to the lowest number of reads in the dataset.
"""
import pysam

threads = snakemake.threads
bams = snakemake.input["bam"]
out_bams = snakemake.output["bam"]
read_count = [pysam.AlignmentFile(x).mapped for x in bams]

# Get the minimum number of reads
min_reads = min(read_count)

# Downsample all BAM files to the minimum number of reads
seed = 42
for i, bam in enumerate(bams):
    out = out_bams[i]
    FLOAT = min_reads / read_count[i]
    pysam.view("--threads", threads,
               "--subsample-seed", seed, 
               "--subsample", FLOAT, 
               "-b", "-o", out, bam)