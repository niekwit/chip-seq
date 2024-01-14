import pandas as pd
import sys
import os
import re

def import_samples():
    """Imports all samples from config/samples.csv
    """
    
    try:
        csv = pd.read_csv("config/samples.csv")
        SAMPLES = csv["Sample"].tolist()
        SAMPLES.extend(csv["Input"].tolist())
        
        # check if samples from samples.csv match fastq files (both R1 and R2) in reads folder
        for sample in SAMPLES:
            r1 = f"reads/{sample}_R1.fastq.gz"
            r2 = f"reads/{sample}_R2.fastq.gz"
            
            assert any(os.path.isfile(x) for x in [r1,r2]), f"ERROR: one or more fastq files (R1/R2) from sample {sample} not found in reads directory"

        return SAMPLES
        
    except FileNotFoundError:
        print("ERROR: config/samples.csv not found")
        sys.exit(1)
        

def import_ip_samples():
    """Imports all IP samples from config/samples.csv
    """
    csv = pd.read_csv("config/samples.csv")
    IP_SAMPLES = csv["Sample"].tolist()
    
    return IP_SAMPLES


def import_input_samples():
    """Imports all input samples from config/samples.csv
    """
    csv = pd.read_csv("config/samples.csv")
    INPUT_SAMPLES = csv["Input"].tolist()
    
    return INPUT_SAMPLES


def import_unique_ip_names(IP_SAMPLES):
    """Imports all unique IP sample names from config/samples.csv
    """
    # remove underscore and trailing numbers from IP sample names
    UNIQUE_IP_NAMES = list(set([re.sub(r"_*\d+$","",x) for x in IP_SAMPLES]))
    
    return UNIQUE_IP_NAMES


def macs2_genome(genome):
    
    if "hg" in genome:
        genome = "hs"
    