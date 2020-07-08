#!/bin/bash

#parse command line options

usage() { echo "Usage: $0 [-d <fast5.dir>] [-r <fastq.dir>] [-g <genome.fa>] [-b <regions.optional.bed>] [-t <threads.int>] [-o <output.dir>] [-l <label>]" 1>&2; exit 1;}

while getopts "hd:r:g:b:t:o:l:" opt; do
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

	b)
	    b=$(readlink -f ${OPTARG})
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

### NANOPOLISH ###

echo "Running nanopolish pipeline to call base modifications"

mkdir ${o}/minimap2
mkdir ${o}/nanopolish

#cat all the fastq in the ${r} directory
now=$(date +"%T")

echo "[${now}] [cat]"

cat ${r}/*.fastq > ${r}/${l}.fastq

now=$(date +"%T")

echo "[${now}] [cat] Done"

#create nanopolish index. The default sequencing summary does not seem to work to speed up indexing, the one from albacore is probably different

echo "[${now}] [nanopolish index]"

docker run -v ${r}/:/$(basename ${r})/ -v ${d}/:/$(basename ${d})/ davidebolo1993/treadmill nanopolish index -v -d /$(basename ${d})/ /$(basename ${r})/${l}.fastq

now=$(date +"%T")

echo "[${now}] [nanopolish index] Done"

#align with minimap2. Running with the default parameter settings and --MD tag enabled for users who want to perform further variant calling (?)

echo "[${now}] [minimap2 -x map-ont]"

docker run -v $(dirname ${g})/:/genome/ -v ${r}/:/$(basename ${r})/ -v ${o}/minimap2/:/$(basename ${o})/minimap2/ davidebolo1993/treadmill minimap2 -ax map-ont -o /$(basename ${o})/minimap2/${l}.sam --MD -t ${t} /genome/$(basename ${g}) /$(basename ${r})/${l}.fastq

now=$(date +"%T")

echo "[${now}] [minimap2 -x map-ont] Done"

#convert sam to bam

echo "[${now}] [samtools view]"

docker run -v ${o}/minimap2/:/$(basename ${o})/minimap2/ davidebolo1993/treadmill samtools view -o /$(basename ${o})/minimap2/${l}.bam -@ ${t} /$(basename ${o})/minimap2/${l}.sam
rm ${o}/minimap2/${l}.sam

now=$(date +"%T")

echo "[${now}] [samtools view] Done"

#sort bam

echo "[${now}] [samtools sort]"

docker run -v ${o}/minimap2/:/$(basename ${o})/minimap2/ davidebolo1993/treadmill samtools sort -o /$(basename ${o})/minimap2/${l}.srt.bam -@ ${t} /$(basename ${o})/minimap2/${l}.bam
rm ${o}/minimap2/${l}.bam

now=$(date +"%T")

echo "[${now}] [samtools sort] Done"

#index

echo "[${now}] [samtools index]"

docker run -v ${o}/minimap2/:/$(basename ${o})/minimap2/ davidebolo1993/treadmill samtools index -@ ${t} /$(basename ${o})/minimap2/${l}.srt.bam

now=$(date +"%T")

echo "[${now}] [samtools index] Done"

#call methylation genome-wide or in regions

echo "[${now}] [nanopolish call-methylation]"

if [ -f "${b}" ]

then

	while IFS=$'\t' read -r chrom start end; do

        	docker run -v $(dirname ${g})/:/genome/ -v ${o}/minimap2/:/$(basename ${o})/minimap2/ -v ${o}/nanopolish/:/$(basename ${o})/nanopolish/ -v ${d}/:/$(basename ${d})/ -v ${r}/:/$(basename ${r})/ davidebolo1993/treadmill nanopolish call-methylation -v --progress -r /$(basename ${r})/${l}.fastq -b /$(basename ${o})/minimap2/${l}.srt.bam -g /genome/$(basename ${g}) -q cpg -t ${t} -w "${chrom}":"$(printf "%'d" ${start})"-"$(printf "%'d" ${end})" > ${o}/nanopolish/${l}.${chrom}":"${start}"-"${end}.methylation_calls.tsv 2> /dev/null

	done < ${b}

else

	docker run -v $(dirname ${g})/:/genome/ -v ${o}/minimap2/:/$(basename ${o})/minimap2/ -v ${o}/nanopolish/:/$(basename ${o})/nanopolish/ -v ${d}/:/$(basename ${d})/ -v ${r}/:/$(basename ${r})/ davidebolo1993/treadmill nanopolish call-methylation -v --progress -r /$(basename ${r})/${l}.fastq -b /$(basename ${o})/minimap2/${l}.srt.bam -g /genome/$(basename ${g}) -q cpg -t ${t} > ${o}/nanopolish/${l}.methylation_calls.tsv 2> /dev/null

fi

now=$(date +"%T")

echo "[${now}] [nanopolish call-methylation] Done"
