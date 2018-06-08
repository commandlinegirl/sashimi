FROM ubuntu:16.04


RUN apt-get update && apt-get install -q -y \
    curl \
    git \
    build-essential \
    libpq-dev \
    python-dev \
    libxml2-dev \
    libxslt1-dev \
    libldap2-dev \
    libsasl2-dev \
    liblzma-dev \
    libffi-dev \
    python-pip \
    bedtools

RUN pip install python-dateutil --upgrade

RUN pip install \
    numpy \
    scipy \
    pandas \
    pybedtools \
    cython


#RUN git clone https://github.com/daler/pybedtools.git
#RUN cd pybedtools && git pull && python setup.py develop && cd

# Download dx-toolkit and source dx at login
RUN curl https://wiki.dnanexus.com/images/files/dx-toolkit-v0.253.0-ubuntu-14.04-amd64.tar.gz | tar xzf - -C /usr/local/ --no-same-owner --no-same-permissions && echo "source /usr/local/dx-toolkit/environment" >> ~/.bashrc

ADD lib /opt

SHELL ["/bin/bash", "-c"]

#ENTRYPOINT ["python"]

WORKDIR /opt
