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

ARG CACHEBUST=1

RUN git clone https://github.com/commandlinegirl/sashimi.git
RUN mv sashimi/src/* /opt

WORKDIR /opt


