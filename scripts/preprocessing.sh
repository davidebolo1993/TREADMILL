#!/bin/bash

#parse command line options

usage() { echo "Usage: $0 [-d <fast5.dir>] [-r <fastq.dir>] [-g <genome.fa>] [-t <threads.int>] [-o <output.dir>] [-l <label>]" 1>&2; exit 1;}

while getopts "hd:r:g:t:o:l:" opt; do
    case "${opt}" in

        d)
            d=$(readlink -f ${OPTARG})
            ;;
        
        r)
            r=$(readlink -f ${OPTARG})
            ;;

        g)
            g=$(readlink -f ${OPTARG})
            ;;

        t)
            t=${OPTARG}
            ;;

        o)
            o=$(readlink -f ${OPTARG})
            ;;

        l)
            l=${OPTARG}
            ;;

        \?)
            echo "Invalid Option: -$OPTARG" 1>&2
            usage
            exit 1
            ;;

        h )
            usage
            exit 0
            ;;
    esac
done

shift $((OPTIND-1))
mkdir ${o}

### Nanopolish ###

echo "Running nanopolish pipeline to call base modifications"

#cat all the fastq in the ${r} directory
now=$(date +"%T")

echo "[${now}] [cat]"

set -x

cat ${r}/*.fastq > ${r}/${l}.fastq

set +x

now=$(date +"%T")

echo "[${now}] [cat] Done"

#create nanopolish index. The default sequencing summary does not seem to work to speed up indexing, the one from Albacore is probably different

echo "[${now}] [nanopolish index]"

set -x

docker run -v ${r}/:/$(basename ${r})/ -v ${d}/:/$(basename ${d})/ -ti davidebolo1993/treadmill nanopolish index -v -d /$(basename ${d})/ /$(basename ${r})/${l}.fastq

set +x

now=$(date +"%T")

echo "[${now}] [nanopolish index] Done"

#align with minimap2. Running with the default parameter settings and --MD tag enabled for further variant calling??

mkdir ${o}/minimap2

now=$(date +"%T")

echo "[${now}] [minimap2]"

set -x

docker run -v $(dirname ${g})/:/genome/ -v ${r}/:/$(basename ${r})/ -v ${o}/minimap2/:/$(basename ${o})/minimap2/ -ti davidebolo1993/treadmill minimap2 -ax map-ont -o /$(basename ${o})/minimap2/${l}.sam --MD -t {t} /genome/$(basename ${g}) /$(basename ${r})/${l}.fastq

set +x

now=$(date +"%T")

echo "[${now}] [minimap2] Done"
