FROM python:2.7.15
LABEL maintainer="Harry Jubb<hj4@sanger.ac.uk>"

# INSTALL OPENBABEL WITH PYTHON BINDINGS
RUN apt-get update && \
  apt-get install -y --no-install-recommends \
  cmake \
  libxml2-dev \
  zlib1g-dev \
  libeigen3-dev \
  libcairo2-dev \
  swig \
  && rm -rf /var/lib/apt/lists/*
RUN mkdir /openbabel
WORKDIR /openbabel
RUN wget http://sourceforge.net/projects/openbabel/files/openbabel/2.4.0/openbabel-openbabel-2-4-0.tar.gz
RUN tar -xzf openbabel-openbabel-2-4-0.tar.gz
RUN mkdir build
WORKDIR build
RUN cmake ../openbabel-openbabel-2-4-0 -DRUN_SWIG=ON -DPYTHON_BINDINGS=ON
RUN make
RUN make install
ENV PYTHONPATH="/usr/local/lib"
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

RUN mkdir /arpeggio
WORKDIR /arpeggio

COPY requirements.txt /arpeggio
RUN pip install --no-cache-dir -r requirements.txt

COPY config.py /arpeggio/
COPY arpeggio.py /arpeggio/
