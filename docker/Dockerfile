FROM ubuntu:20.04
LABEL maintainer="shnegi@ucsc.edu"

# Prevent dpkg from trying to ask any questions, ever
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    gcc \
    git \
    make \
    bzip2 \
    tabix \
    python3 \
    python3-pip \
    libncurses5-dev \
    libncursesw5-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    autoconf \
    build-essential \
    pkg-config \
    apt-transport-https software-properties-common dirmngr gpg-agent \
    && rm -rf /var/lib/apt/lists/*

## bcftools
RUN wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 && \
        tar -xjf bcftools-1.17.tar.bz2 && \
        cd bcftools-1.17 && \
        ./configure && make && make install && \
        cd .. && rm -rf bcftools-1.17.tar.bz2

# python packages
RUN pip3 install --no-cache-dir numpy pandas seaborn matplotlib scipy argparse

# add scripts
WORKDIR /opt/scripts/SegMeth-v1.0
ADD scripts/* /opt/scripts/SegMeth-v1.0

WORKDIR /home