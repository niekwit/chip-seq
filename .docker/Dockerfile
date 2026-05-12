FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="cbe1306d767275c7c93fbf400a1d65b6a956d1aa506be91a6ccbd53ee0a4743f"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/bedtools/intersect/environment.yaml
#   prefix: /conda-envs/8aac1e4743179f2ba704ad36d26b54d1
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bedtools =2.31.1
RUN mkdir -p /conda-envs/8aac1e4743179f2ba704ad36d26b54d1
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/bedtools/intersect/environment.yaml /conda-envs/8aac1e4743179f2ba704ad36d26b54d1/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/bowtie2/align/environment.yaml
#   prefix: /conda-envs/da7819ccf0f2b653dc605d765769a2e6
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bowtie2 =2.5.3
#     - samtools =1.20
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/da7819ccf0f2b653dc605d765769a2e6
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/bowtie2/align/environment.yaml /conda-envs/da7819ccf0f2b653dc605d765769a2e6/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/bowtie2/build/environment.yaml
#   prefix: /conda-envs/398c3dba37080bfd985cf97615815622
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bowtie2 =2.5.3
RUN mkdir -p /conda-envs/398c3dba37080bfd985cf97615815622
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/bowtie2/build/environment.yaml /conda-envs/398c3dba37080bfd985cf97615815622/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/deeptools/bamcoverage/environment.yaml
#   prefix: /conda-envs/e84a94c164a7278b8266143f989704c2
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - deeptools =3.5.5
RUN mkdir -p /conda-envs/e84a94c164a7278b8266143f989704c2
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/deeptools/bamcoverage/environment.yaml /conda-envs/e84a94c164a7278b8266143f989704c2/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/deeptools/computematrix/environment.yaml
#   prefix: /conda-envs/ac43e5777290902d31614a40c46d3b27
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - deeptools =3.5.5
RUN mkdir -p /conda-envs/ac43e5777290902d31614a40c46d3b27
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/deeptools/computematrix/environment.yaml /conda-envs/ac43e5777290902d31614a40c46d3b27/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/deeptools/plotfingerprint/environment.yaml
#   prefix: /conda-envs/7cb095aaa511a7da68018746ff70c626
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - deeptools =3.5.5
RUN mkdir -p /conda-envs/7cb095aaa511a7da68018746ff70c626
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/deeptools/plotfingerprint/environment.yaml /conda-envs/7cb095aaa511a7da68018746ff70c626/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/deeptools/plotheatmap/environment.yaml
#   prefix: /conda-envs/0d2d19e93d8855f411360cb686b17a7d
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - deeptools =3.5.5
RUN mkdir -p /conda-envs/0d2d19e93d8855f411360cb686b17a7d
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/deeptools/plotheatmap/environment.yaml /conda-envs/0d2d19e93d8855f411360cb686b17a7d/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/macs2/callpeak/environment.yaml
#   prefix: /conda-envs/0dcab85fa385f9120d35e4d38b7b2b13
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - macs2 =2.2.9.1
RUN mkdir -p /conda-envs/0dcab85fa385f9120d35e4d38b7b2b13
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/macs2/callpeak/environment.yaml /conda-envs/0dcab85fa385f9120d35e4d38b7b2b13/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/picard/markduplicates/environment.yaml
#   prefix: /conda-envs/bcfbc3642d4958486dba442bafb24402
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - picard =3.1.1
#     - samtools =1.20
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/bcfbc3642d4958486dba442bafb24402
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/picard/markduplicates/environment.yaml /conda-envs/bcfbc3642d4958486dba442bafb24402/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/samtools/faidx/environment.yaml
#   prefix: /conda-envs/a302e5bfdf888bec35ee16e4d2d11f55
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.20
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/a302e5bfdf888bec35ee16e4d2d11f55
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/samtools/faidx/environment.yaml /conda-envs/a302e5bfdf888bec35ee16e4d2d11f55/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/samtools/index/environment.yaml
#   prefix: /conda-envs/42a2caff52fc46f2f25454995edb10df
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.20
RUN mkdir -p /conda-envs/42a2caff52fc46f2f25454995edb10df
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/samtools/index/environment.yaml /conda-envs/42a2caff52fc46f2f25454995edb10df/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/trim_galore/pe/environment.yaml
#   prefix: /conda-envs/975a4a78cc02737b050b54b2645b9ecc
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - trim-galore =0.6.10
RUN mkdir -p /conda-envs/975a4a78cc02737b050b54b2645b9ecc
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/trim_galore/pe/environment.yaml /conda-envs/975a4a78cc02737b050b54b2645b9ecc/environment.yaml

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

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/8aac1e4743179f2ba704ad36d26b54d1 --file /conda-envs/8aac1e4743179f2ba704ad36d26b54d1/environment.yaml && \
    conda env create --prefix /conda-envs/da7819ccf0f2b653dc605d765769a2e6 --file /conda-envs/da7819ccf0f2b653dc605d765769a2e6/environment.yaml && \
    conda env create --prefix /conda-envs/398c3dba37080bfd985cf97615815622 --file /conda-envs/398c3dba37080bfd985cf97615815622/environment.yaml && \
    conda env create --prefix /conda-envs/e84a94c164a7278b8266143f989704c2 --file /conda-envs/e84a94c164a7278b8266143f989704c2/environment.yaml && \
    conda env create --prefix /conda-envs/ac43e5777290902d31614a40c46d3b27 --file /conda-envs/ac43e5777290902d31614a40c46d3b27/environment.yaml && \
    conda env create --prefix /conda-envs/7cb095aaa511a7da68018746ff70c626 --file /conda-envs/7cb095aaa511a7da68018746ff70c626/environment.yaml && \
    conda env create --prefix /conda-envs/0d2d19e93d8855f411360cb686b17a7d --file /conda-envs/0d2d19e93d8855f411360cb686b17a7d/environment.yaml && \
    conda env create --prefix /conda-envs/0dcab85fa385f9120d35e4d38b7b2b13 --file /conda-envs/0dcab85fa385f9120d35e4d38b7b2b13/environment.yaml && \
    conda env create --prefix /conda-envs/bcfbc3642d4958486dba442bafb24402 --file /conda-envs/bcfbc3642d4958486dba442bafb24402/environment.yaml && \
    conda env create --prefix /conda-envs/a302e5bfdf888bec35ee16e4d2d11f55 --file /conda-envs/a302e5bfdf888bec35ee16e4d2d11f55/environment.yaml && \
    conda env create --prefix /conda-envs/42a2caff52fc46f2f25454995edb10df --file /conda-envs/42a2caff52fc46f2f25454995edb10df/environment.yaml && \
    conda env create --prefix /conda-envs/975a4a78cc02737b050b54b2645b9ecc --file /conda-envs/975a4a78cc02737b050b54b2645b9ecc/environment.yaml && \
    conda env create --prefix /conda-envs/77532f99f5fa631ed3f198f2efa2f451 --file /conda-envs/77532f99f5fa631ed3f198f2efa2f451/environment.yaml && \
    conda env create --prefix /conda-envs/d47898d8476a6b0c4616e73266df7377 --file /conda-envs/d47898d8476a6b0c4616e73266df7377/environment.yaml && \
    conda env create --prefix /conda-envs/7a8bab1f3dfe1f9a74db4ef9b7fd29d1 --file /conda-envs/7a8bab1f3dfe1f9a74db4ef9b7fd29d1/environment.yaml && \
    conda env create --prefix /conda-envs/b3c6ffed432a95e3b1381ffeb757799f --file /conda-envs/b3c6ffed432a95e3b1381ffeb757799f/environment.yaml && \
    conda clean --all -y