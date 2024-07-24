import pysam
import pandas as pd
import os

# Get bam files
prededup = snakemake.input["pre"]
dedup = snakemake.input["post"]

threads = str(snakemake.threads)

# Lists to store read counts and convert to df later
samples = []
pre_counts = []
post_counts = []

# count reads pre and post deduplication
for i, f in enumerate(prededup):
    # Get sample name
    samples.append(os.path.basename(f).replace(".bl.sorted.bam", ""))
    
    # Count reads pre deduplication
    pre_counts.append(pysam.view("-@", threads, "-c", "-F", "260", f).strip())
    
    # Count reads post deduplication
    post_counts.append(pysam.view("-@", threads, "-c", "-F", "260", dedup[i]).strip())
    
# Create df
df = pd.DataFrame({
    "sample": samples,
    "pre.dedup_counts": pre_counts,
    "post.dedup_counts": post_counts
})

# save df to csv
df.to_csv(snakemake.output[0], index=False)
