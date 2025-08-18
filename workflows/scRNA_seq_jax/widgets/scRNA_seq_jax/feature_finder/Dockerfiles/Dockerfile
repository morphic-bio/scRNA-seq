FROM debian:bookworm-slim AS builder
ARG VERSION=1.0.0
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -y update && apt-get -y install  build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libgomp1 libglib2.0-dev curl libcairo2-dev libx11-dev libxext-dev libxrender-dev
RUN curl -L https://github.com/morphic-bio/process_features/archive/refs/tags/v.${VERSION}.tar.gz | tar -xz
WORKDIR /process_features-v.${VERSION}
RUN make && mv assignBarcodes /usr/local/bin/assignBarcodes