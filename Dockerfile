#Build with:
#sudo docker build -t davidebolo1993/treadmill .

FROM ubuntu:18.04

# File author/maintainer info
MAINTAINER Davide Bolognini <davidebolognini7@gmail.com>

# Install dependencies
ENV DEBIAN_FRONTEND=noninteractive 
RUN apt-get update && apt-get install -y nano curl git build-essential g++ cmake libz-dev libcurl4-openssl-dev libssl-dev libbz2-dev liblzma-dev libncurses5-dev && apt-get clean
RUN curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda
RUN conda create -y -n treadmillenv python=3.8
RUN echo "source activate treadmillenv" > ~/.bashrc
ENV PATH /miniconda/envs/treadmillenv/bin:$PATH

#get htslib
RUN curl -LO https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 && tar -vxjf htslib-1.11.tar.bz2 && cd htslib-1.11 && make && make install
#get samtools
RUN curl -LO https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && tar -vxjf samtools-1.11.tar.bz2 && cd samtools-1.11 && make && make install
#get bcftools
RUN curl -LO https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2 && tar -vxjf bcftools-1.11.tar.bz2 && cd bcftools-1.11 && make && make install
#install bedtools through miniconda
RUN conda install -y -n treadmillenv -c bioconda bedtools
#install also r-base and its dependencies
RUN conda install -y -n treadmillenv -c r r-base r-data.table r-scales r-tseries r-gtools r-rcolorbrewer r-ggplot2 r-rjson r-plyr
RUN conda install -y -n treadmillenv -c conda-forge r-ggforce
RUN conda install -y -n treadmillenv -c bioconda r-optparse r-changepoint
RUN R -e "install.packages('ggrepel',dependencies=TRUE, repos='http://cran.rstudio.com/')"
#get TREADMILL
RUN git clone --recursive https://github.com/davidebolo1993/TREADMILL.git && cd TREADMILL && ./configure && python setup.py install

#also install deepsignal in order to calculate methylation frequencies
RUN conda create -y -n deepsignalenv python=3.6
ENV PATH /miniconda/envs/deepsignalenv/bin:$PATH
RUN conda install -y -n deepsignalenv -c bioconda ont-tombo ont-fast5-api
RUN conda install -y -n deepsignalenv -c anaconda tensorflow==1.13.1
RUN pip install deepsignal

#Pull with:
#sudo docker pull davidebolo1993/treadmill

#Then run:
#sudo docker run davidebolo1993/treadmill TREADMILL --help

#Or load the environment
#sudo docker run -ti davidebolo1993/treadmill
#$(treadmillenv) TREADMILL --help