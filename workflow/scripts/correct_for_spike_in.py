"""
Correct all (downsampled) BAM file for spike-in reads:
    correction factor:
        lowest number of spike-inreads in the dataset / number of spike-in reads in the BAM file
"""
import pysam

threads = snakemake.threads
bams = snakemake.input["bam"]
out_bams = snakemake.output["bam"]
si_bams = snakemake.input["si_bam"]
si_read_counts = [pysam.AlignmentFile(x).mapped for x in si_bams]

# Get the minimum number of spike-in reads
min_si_reads = min(si_read_counts)

# Correct all BAM files for spike-in reads
seed = 1234 # Should be different from the seed used in downsample.py
for i, bam in enumerate(bams):
    FLOAT = min_si_reads / si_read_counts[i]
    pysam.view("--threads", threads,
               "--subsample-seed", seed,
               "--subsample", FLOAT,
               "-b", "-o", out_bams[i], bam)
