FROM ubuntu:22.04

ENV R_VERSION=4.5.0
RUN apt-get update && apt-get install -y curl
RUN curl -O https://cdn.rstudio.com/r/ubuntu-2204/pkgs/r-${R_VERSION}_1_amd64.deb
RUN apt-get update && apt-get install -y ./r-${R_VERSION}_1_amd64.deb
ENV TZ=America/Los_Angeles
ENV DEBIAN_FRONTEND=noninteractive
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R
RUN ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript

# Install Seurat's system dependencies
RUN apt-get update
RUN apt-get install -y \
    build-essential \
    git \
    libabsl-dev\ 
    libboost-all-dev \
    libcurl4-gnutls-dev \
    libfftw3-dev \
    libfontconfig1-dev \
    libfreetype-dev \
    libfribidi-dev \
    libgit2-dev \
    libgsl-dev \
    libharfbuzz-dev \
    libhdf5-dev \
    libjpeg-dev \
    libpng-dev \
    libssl-dev \
    libtiff5-dev \
    libxml2-dev \
    libcairo-dev \
    libgmp-dev \
    openjdk-8-jdk \
    pkg-config \
    python3-dev \
    python3-pip \
    wget

RUN R -e 'install.packages("devtools", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("reticulate", version = "1.34.0", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("pals", version = "1.8", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("clustree", version = "0.5.0", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("ggplot2", version = "3.4.1", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("Matrix", version = "1.6-4", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("matrixStats", version = "0.63.0", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("SeuratObject", version = "5.1.0", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("Seurat", version = "5.3.0", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_version("argparse", version = "2.2.1", repos = "http://cran.us.r-project.org")'
# Install remotes to install seurat-disk
RUN R -e 'install.packages("remotes", repos = "http://cran.us.r-project.org")'
RUN R -e 'remotes::install_github("mojaveazure/seurat-disk")'

# Install BiocManager, then SeuratDisk and its HDF5 dependencies
RUN R -e 'install.packages("BiocManager", repos = "http://cran.us.r-project.org")'
RUN R -e 'BiocManager::install("Signac", update=FALSE, ask=FALSE, dependencies=TRUE, Ncpus=8, verbose=TRUE)'
RUN R -e "remotes::install_cran('anndata')"
RUN R -e "remotes::install_bioc('DropletUtils')"
RUN R -e "remotes::install_bioc('scDblFinder')"
#RUN pip install anndata==0.9.2
#RUN R -e "reticulate::import('anndata')"
# Install Python dependencies
RUN pip3 install --no-cache-dir --upgrade pip setuptools wheel 
RUN pip3 install --no-cache-dir numpy==1.24.3 anndata==0.10.3 numba==0.58.1 pandas==2.1.1 rpy2 pyroe scanpy==1.9.6 plotly==5.22.0 kaleido==0.2.1
# Copy the R script into the Docker image
COPY filterEmptyCells.R /usr/local/bin/filterEmptyCells.R
COPY findValidCells.R /usr/local/bin/findValidCells.R
COPY runFindValidCells.sh /usr/local/bin/runFindValidCells.sh
COPY read_counts.py /usr/local/bin/read_counts.py
COPY combineFilters.py /usr/local/bin/combineFilters.py
COPY filter.sh /usr/local/bin/filter.sh
