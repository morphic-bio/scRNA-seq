ARG KALLISTO_VERSION=0.50.1
ARG BUS_VERSION=0.43.2
ARG KB_PYTHON_VERSION=0.28.0
ARG STAR_VERSION=2.7.11a
ARG PYROE_VERSION=0.9.3
ARG ALEVINFRY_VERSION=0.8.2
ARG SIMPLEAF_VERSION=0.16.2
ARG PISCEM_VERSION=0.7.3
ARG SALMON_VERSION=1.10.1

#builders

#star builder from biodepot - has multiple versions of star
FROM biodepot/starbins:latest AS star_builder

#alevin-fry
FROM rust:1.76.0-slim-bookworm AS alevin_builder
ARG DEBIAN_FRONTEND=noninteractive

ARG ALEVINFRY_VERSION
ARG SIMPLEAF_VERSION
ARG PISCEM_VERSION
ARG SALMON_VERSION
ARG PYROE_VERSION

RUN apt-get update && apt-get install -y curl build-essential g++ libbz2-dev libz-dev curl libboost-all-dev liblzma-dev ca-certificates zlib1g-dev libcurl4-openssl-dev curl unzip autoconf apt-transport-https gnupg software-properties-common wget cmake git

RUN git clone --recursive https://github.com/COMBINE-lab/piscem.git
RUN cd piscem* && git checkout tags/v${PISCEM_VERSION}
RUN  cd piscem* && cargo build --release
RUN git clone https://github.com/COMBINE-lab/alevin-fry.git
RUN cd alevin-fry && git checkout tags/v${ALEVINFRY_VERSION}
RUN cd alevin-fry && cargo build --release
RUN git clone https://github.com/COMBINE-lab/simpleaf.git
RUN cd simpleaf && git checkout tags/v${SIMPLEAF_VERSION}
RUN cd simpleaf && cargo build --release
#RUN python3 -m pip install pyroe==${PYROE_VERSION} --break-system-packages
RUN curl -L https://github.com/COMBINE-lab/salmon/archive/refs/tags/v${SALMON_VERSION}.tar.gz |  tar -zvxf -
RUN cd salmon-${SALMON_VERSION} && \
    mkdir build && \
    cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local/salmon && make && make install

#build kallisto and bustools
FROM debian:bookworm-slim AS kallisto_builder
ARG KALLISTO_VERSION
ARG BUS_VERSION

RUN apt-get update && apt-get install -y build-essential cmake autoconf2.69 automake curl \
 zlib1g-dev libhdf5-dev libbz2-dev liblzma-dev

#install kallisto
RUN curl -L https://github.com/pachterlab/kallisto/archive/refs/tags/v${KALLISTO_VERSION}.tar.gz | tar -zvxf - 
RUN cd kallisto*/ext/htslib && autoheader && autoconf
RUN cd kallisto* && mkdir -p build && cd build && cmake .. -DUSE_HDF5=ON && make && make install

#install bustools
RUN curl -L https://github.com/BUStools/bustools/archive/refs/tags/v${BUS_VERSION}.tar.gz | tar -zvxf - 
RUN cd bustools* && mkdir -p build && cd build && cmake .. && make && make install


#put them all together

FROM python:3.11.13-slim-bookworm
ARG DEBIAN_FRONTEND=noninteractive
ARG KB_PYTHON_VERSION
ARG PYROE_VERSION
ARG STAR_VERSION
ARG LD_LIBRARY_PATH
RUN apt-get update && \
    apt-get install -y wget time libhwloc15 libhdf5-dev zlib1g-dev libbz2-dev liblzma-dev curl && \
    apt-get remove -y curl && apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*
RUN pip install --upgrade pip setuptools wheel
RUN pip install numba==0.57.1 numcodecs==0.11.0 numpy==1.24.3 pandas==1.5.3 anndata==0.11.4 kb-python==${KB_PYTHON_VERSION} pyroe==${PYROE_VERSION} --break-system-packages

#Note that bookworm uses autoconf 2.71 which is buggy and causes problems with the htslib compile

#put all the copy commands at the end to avoid rebuilding everything
#still need a wrapper because of comma delimited input

#kallisto and bustools
#the copy is necessary in order to get arm64 linux binaries when necessary
COPY --from=kallisto_builder /usr/local/bin/kallisto /usr/local/lib/python*/site-packages/kb_python/bins/linux/kallisto/kallisto
COPY --from=kallisto_builder /usr/local/bin/bustools /usr/local/lib/python*/site-packages/kb_python/bins/linux/bustools/bustools
RUN ln -s /usr/local/lib/python3*/site-packages/kb_python/bins/linux/kallisto/kallisto /usr/local/bin/kallisto
RUN ln -s /usr/local/lib/python3*/site-packages/kb_python/bins/linux/bustools/bustools /usr/local/bin/bustools
COPY runKallisto.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/runKallisto.sh
RUN pip install --upgrade pandas

#star_solo
COPY --from=star_builder /bin/$STAR_VERSION /usr/local/bin
COPY runstar.sh /usr/local/bin/runstar.sh
COPY version.sh /usr/local/bin/version.sh
ENV PATH=$PATH:/usr/local/bin/Linux_x86_64_static

#alevin-fry
COPY --from=alevin_builder /alevin-fry /usr/local/alevin-fry
COPY --from=alevin_builder /simpleaf /usr/local/simpleaf
COPY --from=alevin_builder /usr/local/salmon/ /usr/local/salmon
COPY --from=alevin_builder /piscem/target/release/piscem /usr/local/bin/piscem
ENV LD_LIBRARY_PATH=/usr/local/salmon/lib:$LD_LIBRARY_PATH 
ENV PATH=/usr/local/salmon/bin:$PATH
ENV PATH=/usr/local/alevin-fry/target/release:$PATH
ENV PATH=/usr/local/simpleaf/target/release:$PATH
COPY runMultiAlign.sh /usr/local/bin/runMultiAlign.sh
COPY runMultiIndex.sh /usr/local/bin/runMultiIndex.sh
