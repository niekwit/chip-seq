import re

# Retrieve parameters from Snakemake
window_size = int(snakemake.params["size"]) * 1000 
bed_file = snakemake.output[0]
chrom_sizes_file = snakemake.input[0]

# Read chromosome sizes into a dictionary
chrom_sizes = {}
with open(chrom_sizes_file) as f:
    for line in f:
        chrom, size = line.strip().split()
        
        if re.match("^GL|^KI", chrom):
            continue
        
        chrom_sizes[chrom] = int(size)

# Create windows for each chromosome and write to a BED file
# If the last window goes beyond the chromosome size, it is truncated
counter = 1
with open(bed_file, "w") as f:
    for chrom, size in chrom_sizes.items():
        for start in range(0, size, window_size):
            name = f"window_{counter}"
            end = min(start + window_size, size)
            f.write(f"{chrom}\t{start + 1}\t{end}\t{name}\n")
            counter += 1
            
print(f"{counter -1} genomic windows ({window_size} each) created in {bed_file}...")