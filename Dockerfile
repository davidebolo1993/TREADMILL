#Build with:
#sudo docker build -t davidebolo1993/treadmill .

FROM ubuntu:18.04

# File author/maintainer info
MAINTAINER Davide Bolognini <davidebolognini7@gmail.com>

# Install dependencies
RUN apt-get update && apt-get install -y nano curl git wget build-essential g++ cmake zlib1g-dev && apt-get clean
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
RUN bash Miniconda-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda
RUN conda create -y -n treadmillenv python=3.6
RUN echo "source activate treadmillenv" > ~/.bashrc
ENV PATH /miniconda/envs/treadmillenv/bin:$PATH
RUN git clone --recursive https://github.com/jts/nanopolish.git && cd nanopolish && make all
ENV PATH nanopolish:$PATH
ENV PATH nanopolish/scripts:$PATH
RUN conda install -y -n treadmillenv -c bioconda samtools bcftools bedtools bedops minimap2 ngmlr last pysam pyfaidx cyvcf2 pybedtools numpy longshot
RUN conda install -y -n treadmillenv -c r r
RUN pip install --user whatshap
ENV PATH /root/.local/bin:$PATH
#more to go

#Pull with:
#sudo docker pull davidebolo1993/treadmill
