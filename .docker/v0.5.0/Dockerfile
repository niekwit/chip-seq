FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="7354850eac91724b4c9106327961656c4d38ad1fa6fdd1e73c7b482e23d005be"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/bedtools/intersect/environment.yaml
#   prefix: /conda-envs/4bc1c7d4971260f38b3717e4c2f2ddf4
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bedtools =2.31.1
RUN mkdir -p /conda-envs/4bc1c7d4971260f38b3717e4c2f2ddf4
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/bedtools/intersect/environment.yaml /conda-envs/4bc1c7d4971260f38b3717e4c2f2ddf4/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/bowtie2/align/environment.yaml
#   prefix: /conda-envs/79f24222e503aa3928bb131a6d96239b
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bowtie2 =2.5.4
#     - samtools =1.20
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/79f24222e503aa3928bb131a6d96239b
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/bowtie2/align/environment.yaml /conda-envs/79f24222e503aa3928bb131a6d96239b/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/bowtie2/build/environment.yaml
#   prefix: /conda-envs/84012557432533881e14141abce4094f
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bowtie2 =2.5.4
RUN mkdir -p /conda-envs/84012557432533881e14141abce4094f
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/bowtie2/build/environment.yaml /conda-envs/84012557432533881e14141abce4094f/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/deeptools/bamcoverage/environment.yaml
#   prefix: /conda-envs/4a098defd6d43db31846a0f3e27479bb
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - deeptools =3.5.5
RUN mkdir -p /conda-envs/4a098defd6d43db31846a0f3e27479bb
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/deeptools/bamcoverage/environment.yaml /conda-envs/4a098defd6d43db31846a0f3e27479bb/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/fastqc/environment.yaml
#   prefix: /conda-envs/90ea9a0a5f4de478eb6aa43932774175
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.12.1
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/90ea9a0a5f4de478eb6aa43932774175
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/fastqc/environment.yaml /conda-envs/90ea9a0a5f4de478eb6aa43932774175/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/macs2/callpeak/environment.yaml
#   prefix: /conda-envs/ff66d0870b4e454f6c8e413316e47391
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - macs2 =2.2.9.1
RUN mkdir -p /conda-envs/ff66d0870b4e454f6c8e413316e47391
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/macs2/callpeak/environment.yaml /conda-envs/ff66d0870b4e454f6c8e413316e47391/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/multiqc/environment.yaml
#   prefix: /conda-envs/5f5723fa7ebbd67a3632b0259e715b9f
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - multiqc =1.24.1
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/5f5723fa7ebbd67a3632b0259e715b9f
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/multiqc/environment.yaml /conda-envs/5f5723fa7ebbd67a3632b0259e715b9f/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/picard/markduplicates/environment.yaml
#   prefix: /conda-envs/8d28ed6a140e6a5ec677be879fa1616c
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - picard =3.2.0
#     - samtools =1.20
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/8d28ed6a140e6a5ec677be879fa1616c
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/picard/markduplicates/environment.yaml /conda-envs/8d28ed6a140e6a5ec677be879fa1616c/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/samtools/faidx/environment.yaml
#   prefix: /conda-envs/b29db9b4576aa24f6ca5cdabda9911b4
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.20
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/b29db9b4576aa24f6ca5cdabda9911b4
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/samtools/faidx/environment.yaml /conda-envs/b29db9b4576aa24f6ca5cdabda9911b4/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/samtools/index/environment.yaml
#   prefix: /conda-envs/adec7d399fa9492178d9caa96f93852d
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.20
RUN mkdir -p /conda-envs/adec7d399fa9492178d9caa96f93852d
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/samtools/index/environment.yaml /conda-envs/adec7d399fa9492178d9caa96f93852d/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/trim_galore/pe/environment.yaml
#   prefix: /conda-envs/7c8f9db6001702ef40e5b666df2ea0f5
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - trim-galore =0.6.10
RUN mkdir -p /conda-envs/7c8f9db6001702ef40e5b666df2ea0f5
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.3.0/bio/trim_galore/pe/environment.yaml /conda-envs/7c8f9db6001702ef40e5b666df2ea0f5/environment.yaml

# Conda environment:
#   source: workflow/envs/R.yaml
#   prefix: /conda-envs/77532f99f5fa631ed3f198f2efa2f451
#   name: r
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - r-base=4.3
#     - r-tidyverse=2.0.0
#     - r-cowplot=1.1.2
#     - r-reshape2=1.4.4
#     - r-rcolorbrewer=1.1_3
#     - r-ggrepel=0.9.4
#     - r-scales=1.3.0
#     - bioconductor-chipseeker=1.38.0
#     - bioconductor-genomicfeatures=1.54.1
#     - bioconductor-rtracklayer=1.62.0
#     - r-viridislite=0.4.2
RUN mkdir -p /conda-envs/77532f99f5fa631ed3f198f2efa2f451
COPY workflow/envs/R.yaml /conda-envs/77532f99f5fa631ed3f198f2efa2f451/environment.yaml

# Conda environment:
#   source: workflow/envs/deeptools.yaml
#   prefix: /conda-envs/d47898d8476a6b0c4616e73266df7377
#   name: deeptools
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - deeptools=3.5.2
#     - python=3.10
#     - pysam=0.21.0
#     - matplotlib=3.7.2
#     - pandas=2.0.0
#     - seaborn=0.12.2
#     - numpy=1.23.5
#     - wiggletools=1.2.11
#     - ucsc-wigtobigwig=447
#     - bedtools=2.31.1
#     - pybedtools 0.10.0
RUN mkdir -p /conda-envs/d47898d8476a6b0c4616e73266df7377
COPY workflow/envs/deeptools.yaml /conda-envs/d47898d8476a6b0c4616e73266df7377/environment.yaml

# Conda environment:
#   source: workflow/envs/plotting.yaml
#   prefix: /conda-envs/7a8bab1f3dfe1f9a74db4ef9b7fd29d1
#   name: plotting
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - r-tidyverse=2.0.0
#     - r-cowplot=1.1.1
#     - pysam=0.22.0
#     - pandas=2.1.1
RUN mkdir -p /conda-envs/7a8bab1f3dfe1f9a74db4ef9b7fd29d1
COPY workflow/envs/plotting.yaml /conda-envs/7a8bab1f3dfe1f9a74db4ef9b7fd29d1/environment.yaml

# Conda environment:
#   source: workflow/envs/trim_galore.yml
#   prefix: /conda-envs/b3c6ffed432a95e3b1381ffeb757799f
#   name: trim_galore
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - trim-galore=0.6.10
#     - conda-forge::pigz=2.8
RUN mkdir -p /conda-envs/b3c6ffed432a95e3b1381ffeb757799f
COPY workflow/envs/trim_galore.yml /conda-envs/b3c6ffed432a95e3b1381ffeb757799f/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/4bc1c7d4971260f38b3717e4c2f2ddf4 --file /conda-envs/4bc1c7d4971260f38b3717e4c2f2ddf4/environment.yaml && \
    mamba env create --prefix /conda-envs/79f24222e503aa3928bb131a6d96239b --file /conda-envs/79f24222e503aa3928bb131a6d96239b/environment.yaml && \
    mamba env create --prefix /conda-envs/84012557432533881e14141abce4094f --file /conda-envs/84012557432533881e14141abce4094f/environment.yaml && \
    mamba env create --prefix /conda-envs/4a098defd6d43db31846a0f3e27479bb --file /conda-envs/4a098defd6d43db31846a0f3e27479bb/environment.yaml && \
    mamba env create --prefix /conda-envs/90ea9a0a5f4de478eb6aa43932774175 --file /conda-envs/90ea9a0a5f4de478eb6aa43932774175/environment.yaml && \
    mamba env create --prefix /conda-envs/ff66d0870b4e454f6c8e413316e47391 --file /conda-envs/ff66d0870b4e454f6c8e413316e47391/environment.yaml && \
    mamba env create --prefix /conda-envs/5f5723fa7ebbd67a3632b0259e715b9f --file /conda-envs/5f5723fa7ebbd67a3632b0259e715b9f/environment.yaml && \
    mamba env create --prefix /conda-envs/8d28ed6a140e6a5ec677be879fa1616c --file /conda-envs/8d28ed6a140e6a5ec677be879fa1616c/environment.yaml && \
    mamba env create --prefix /conda-envs/b29db9b4576aa24f6ca5cdabda9911b4 --file /conda-envs/b29db9b4576aa24f6ca5cdabda9911b4/environment.yaml && \
    mamba env create --prefix /conda-envs/adec7d399fa9492178d9caa96f93852d --file /conda-envs/adec7d399fa9492178d9caa96f93852d/environment.yaml && \
    mamba env create --prefix /conda-envs/7c8f9db6001702ef40e5b666df2ea0f5 --file /conda-envs/7c8f9db6001702ef40e5b666df2ea0f5/environment.yaml && \
    mamba env create --prefix /conda-envs/77532f99f5fa631ed3f198f2efa2f451 --file /conda-envs/77532f99f5fa631ed3f198f2efa2f451/environment.yaml && \
    mamba env create --prefix /conda-envs/d47898d8476a6b0c4616e73266df7377 --file /conda-envs/d47898d8476a6b0c4616e73266df7377/environment.yaml && \
    mamba env create --prefix /conda-envs/7a8bab1f3dfe1f9a74db4ef9b7fd29d1 --file /conda-envs/7a8bab1f3dfe1f9a74db4ef9b7fd29d1/environment.yaml && \
    mamba env create --prefix /conda-envs/b3c6ffed432a95e3b1381ffeb757799f --file /conda-envs/b3c6ffed432a95e3b1381ffeb757799f/environment.yaml && \
    mamba clean --all -y