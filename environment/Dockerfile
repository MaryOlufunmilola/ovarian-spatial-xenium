# hash:sha256:3251fffbcd17acacbf9e6710999b0c65376424b3aa1e1b097792ebde5a7bbe3c
FROM registry.codeocean.com/codeocean/py-r:python3.12.8-R4.4.2-JupyterLab4.3.5-RStudiorstudio-server-2024.12.1-563-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libbz2-dev=1.0.8-5build1 \
        libcairo2-dev=1.16.0-5ubuntu2.1 \
        libcurl4-openssl-dev \
        libfontconfig1-dev=2.13.1-4.2ubuntu5 \
        libfreetype6-dev=2.11.1+dfsg-1ubuntu0.3 \
        libfribidi-dev=1.0.8-2ubuntu3.1 \
        libgeos-dev=3.10.2-1 \
        libglpk-dev=5.0-1 \
        libharfbuzz-dev=2.7.4-1ubuntu3.2 \
        libhdf5-dev=1.10.7+repack-4ubuntu2 \
        libjpeg-dev=8c-2ubuntu10 \
        liblzma-dev=5.2.5-2ubuntu1.1 \
        libnlopt-dev=2.7.1-3build1 \
        libpng-dev=1.6.37-3ubuntu0.5 \
        libssl-dev=3.0.2-0ubuntu1.25 \
        libtiff5-dev=4.3.0-6ubuntu0.13 \
        libxml2-dev=2.9.13+dfsg-1ubuntu0.12 \
        libxt-dev=1:1.2.1-1 \
        zlib1g-dev=1:1.2.11.dfsg-2ubuntu9.2 \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("RColorBrewer", "1.1-3")' \
    && Rscript -e 'remotes::install_version("Seurat", "5.5.1")' \
    && Rscript -e 'remotes::install_version("SeuratObject")' \
    && Rscript -e 'remotes::install_version("arrow", "24.0.0")' \
    && Rscript -e 'remotes::install_version("devtools", "2.5.2")' \
    && Rscript -e 'remotes::install_version("dplyr", "1.2.1")' \
    && Rscript -e 'remotes::install_version("future", "1.70.0")' \
    && Rscript -e 'remotes::install_version("ggplot2", "4.0.3")' \
    && Rscript -e 'remotes::install_version("ggrepel", "0.9.8")' \
    && Rscript -e 'remotes::install_version("harmony")' \
    && Rscript -e 'remotes::install_version("pacman", "0.5.1")' \
    && Rscript -e 'remotes::install_version("parallelly", "1.47.0")' \
    && Rscript -e 'remotes::install_version("pheatmap", "1.0.13")' \
    && Rscript -e 'remotes::install_version("purrr", "1.2.2")' \
    && Rscript -e 'remotes::install_version("scCustomize", "3.3.0")' \
    && Rscript -e 'remotes::install_version("ggtext", "0.1.2")' \
    && Rscript -e 'remotes::install_version("svglite", "2.2.2")' \
    && Rscript -e 'remotes::install_version("hdf5r", "1.3.12")' \
    && Rscript -e 'remotes::install_version("officer")' \
    && Rscript -e 'remotes::install_version("rvg")' \
    && Rscript -e 'remotes::install_version("tidyverse")' \
    && Rscript -e 'remotes::install_version("BiocManager", "1.30.27")'

RUN Rscript -e 'remotes::install_github( \
        "immunogenomics/harmony", \
        ref = "df19af23ae0639bd6ea2da63898f973f08c85862")' \
    && Rscript -e 'remotes::install_github( \
        "immunogenomics/presto", \
        ref = "a24772a135c7895a8183b007376050556c60a05b")' \
    && Rscript -e 'remotes::install_github( \
        "jokergoo/ComplexHeatmap", \
        ref = "ad11b26afe8c63dc05b20035503352ca740f602e")' \
    && Rscript -e 'remotes::install_github( \
        "jokergoo/circlize", \
        ref = "9b2157843eba242f09adc66048e3470a1561f759")'

RUN conda install -y -c conda-forge -c bioconda openssl zlib \
        bioconductor-rhdf5lib \
        bioconductor-rhdf5 \
        bioconductor-hdf5array \
        bioconductor-glmgampoi \
    && Rscript -e "BiocManager::install(c( \
        'Biobase', \
        'BiocGenerics', \
        'BiocParallel', \
        'BiocNeighbors', \
        'SingleR', \
        'UCell' \      
    ), update=FALSE, ask=FALSE)" \
    && Rscript -e "install.packages('NMF', repos='https://cloud.r-project.org')" \
    && Rscript -e "devtools::install_github('jinworks/CellChat', upgrade='never')"