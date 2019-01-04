FROM ubuntu:18.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -q -y \
    build-essential \
    libpq-dev \
    python-dev \
    python-pip \
    bedtools \
    liblzma-dev \
    python-biopython \
    python-matplotlib \
    python-reportlab \
    python-tk \
    r-base-core \
    zlib1g-dev \
    curl \
    git \
  && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite(c('DNAcopy', 'cghFLasso'))"

RUN pip install python-dateutil --upgrade \
  && pip install \
    numpy \
    scipy \
    pandas \
    pybedtools \
    future \
    futures \
    pyfaidx \
    pysam \
    hmmlearn \
    cnvkit==0.9.5

# Let matplotlib build its font cache
RUN cnvkit.py version

SHELL ["/bin/bash", "-c"]

RUN curl https://wiki.dnanexus.com/images/files/dx-toolkit-v0.253.0-ubuntu-14.04-amd64.tar.gz | tar xzf - -C /usr/local/ --no-same-owner --no-same-permissions \
  && source /usr/local/dx-toolkit/environment

ARG CACHEBUST=1

#RUN git clone https://github.com/commandlinegirl/sashimi.git
RUN git clone -b master --single-branch https://github.com/commandlinegirl/sashimi.git
RUN mv sashimi/src/* /opt

WORKDIR /opt

#ENTRYPOINT cnvkit.py
