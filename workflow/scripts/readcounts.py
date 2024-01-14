import pysam
import pandas as pd

# get bam files
prededup = snakemake.input["pre"]
dedup = snakemake.input["post"]

# nested list to store read counts and convert to df later
df = [] 

# count reads pre and post deduplication
for i in range(len(dedup)):
    # get sample name
    sample = snakemake.wildcards["sample"][i]
    
    # count reads pre deduplication
    pre_count = pysam.view("-@", str(snakemake.threads), "-c", "-F", "260", prededup[i])
    
    # count reads post deduplication
    post_count = pysam.view("-@", str(snakemake.threads), "-c", "-F", "260", dedup[i])
    
    # add data to df
    df.append([sample, pre_count, "pre-deduplication"])
    df.append([sample, post_count, "post-deduplication"])

# convert df to dataframe
df = pd.DataFrame(df, columns=["sample","count","type"])

# save df to csv
df.to_csv(snakemake.output[0], index=False)


