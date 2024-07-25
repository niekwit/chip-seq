import pybedtools
import pandas as pd


bed_files = snakemake.input["beds"]
conditions = snakemake.params["conditions"]
out_bed_files = snakemake.output
k = snakemake.params["keep"]
extend_by = snakemake.params["extend_by"]
max_size = snakemake.params["max_size"]
chrom_sizes = snakemake.input["chrom_sizes"]

# Load chrom sizes into dictionary
chrom_sizes = {}
with open(chrom_sizes) as f:
    for line in f:
        (chrom, size) = line.split()
        chrom_sizes[chrom] = int(size)

for condition in conditions:
    beds = [bed for bed in bed_files if condition in bed]
    
    # Load each individual bed file into a pandas data frame
    ind_peaks = []
    for bed in beds:
        df = pd.read_csv(bed, sep="\t", header=None, low_memory=False)
        ind_peaks.append(df)
    
    # Get overlapping regions between all bed files
    x = pybedtools.BedTool()
    input_list = [pybedtools.BedTool(bed).fn for bed in beds]
    consensus_peaks = x.multi_intersect(i=input_list).to_dataframe()
    
    # Get total number of overlapping regions before filtering
    total = len(consensus_peaks)
    
    # Filter out regions that are not overlapping in at least k bed files
    # name column contains the number of overlapping regions
    consensus_peaks = consensus_peaks[consensus_peaks["name"] >= k]
    
    # Number of peaks not in k bed files
    skipped_peaks = total - len(consensus_peaks)
    
    # Instead of only the overlapping region, get the region that contains 
    # the whole area of the overlapping regions
    regions = []
    extended_peaks = 0
    for row in consensus_peaks.itertuples():
        chrom = row.chrom
        start = row.start
        end = row.end
    
        starts = []
        ends = []
        for df in ind_peaks:
                df = df[(df[0] == chrom) & (df[1] <= int(end)) & (df[2] >= int(start))]
                if not df.empty:
                    starts.append(df[1].min())
                    ends.append(df[2].max())
        start = min(starts)
        end = max(ends)
    
        # Check if regions are short enough to extend
        if end - start < max_size:
            start = start - extend_by
            if start < 1:
                start = 1
            # Do not extend end beyond chromosome boundary
            end = end + extend_by
            if end > chrom_sizes[chrom]:
                end = chrom_sizes[chrom]
            extended_peaks += 1
            
        regions.append([chrom, start, end])
    
    log = [x for x in snakemake.log if condition in x][0]
    with open(log, "a") as log:
        log.write(f"Total number of peaks analysed: {total}\n")
        log.write(f"Skipped peaks: {skipped_peaks}\n")
        log.write(f"Overlapping peaks: {len(consensus_peaks)}\n")
        log.write(f"Extended peaks: {extended_peaks}\n")