FROM ubuntu:16.04

RUN apt-get update && apt-get install -q -y \
    curl \
    git \
    build-essential \
    libpq-dev \
    python-dev \
    python-pip \
    bedtools \
  && rm -rf /var/lib/apt/lists/*

RUN pip install python-dateutil --upgrade \
  && pip install \
    numpy \
    scipy \
    pandas \
    pybedtools

SHELL ["/bin/bash", "-c"]

RUN curl https://wiki.dnanexus.com/images/files/dx-toolkit-v0.253.0-ubuntu-14.04-amd64.tar.gz | tar xzf - -C /usr/local/ --no-same-owner --no-same-permissions \
  && source /usr/local/dx-toolkit/environment

ARG CACHEBUST=1

RUN git clone -b add_smooth_column --single-branch https://github.com/commandlinegirl/sashimi.git
RUN mv sashimi/src/* /opt

WORKDIR /opt


