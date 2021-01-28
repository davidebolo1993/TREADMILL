#Build with:
#sudo docker build -t davidebolo1993/treadmill .

FROM ubuntu:18.04

# File author/maintainer info
MAINTAINER Davide Bolognini <davidebolognini7@gmail.com>
# Install dependencies
RUN apt-get update && apt-get install -y nano curl git wget build-essential g++ cmake zlib1g-dev gfortran && apt-get clean
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
RUN bash Miniconda-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda
RUN conda create -y -n treadmillenv python=3.8
RUN echo "source activate treadmillenv" > ~/.bashrc
ENV PATH /miniconda/envs/treadmillenv/bin:$PATH
RUN git clone --recursive https://github.com/jts/nanopolish.git && cd nanopolish && make all
ENV PATH nanopolish:$PATH
ENV PATH nanopolish/scripts:$PATH
RUN conda install -y -n treadmillenv -c bioconda samtools bcftools bedtools minimap2 ngmlr pysam pyfaidx pybedtools mappy numpy
RUN conda install -y -n treadmillenv -c r r
RUN conda install -y -n treadmillenv -c bioconda bioconductor-karyoploter r-optparse r-ggrepel
RUN conda install -y -n treadmillenv -c r r-rjson r-ggplot2 r-plyr r-data.table 
RUN conda install -y -n treadmillenv -c conda-forge r-ggforce 
RUN git clone --recursive https://github.com/davidebolo1993/TREADMILL && cd TREADMILL && ./configure && python setup.py install #this takes care of all the other python modules
#Pull with:
#sudo docker pull davidebolo1993/treadmill
