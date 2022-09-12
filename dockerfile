FROM continuumio/miniconda3

RUN apt-get update -y && apt install build-essential -y && rm -rf /var/lib/apt/lists/*
RUN pip install numpy cython
RUN pip install scgenome
