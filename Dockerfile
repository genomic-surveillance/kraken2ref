FROM python:3.12

RUN apt-get update -y && apt-get -y install wget curl vim git libncurses-dev bwa samtools libgsl-dev
ADD kraken2ref/ /home/kraken2ref/kraken2ref
ADD .git /home/kraken2ref/.git
COPY pyproject.toml /home/kraken2ref/
COPY setup.cfg  /home/kraken2ref/
COPY requirements.txt  /home/kraken2ref/

RUN cd  /home/kraken2ref/ && pip install .
WORKDIR /home