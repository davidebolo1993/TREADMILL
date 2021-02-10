#!/bin/bash

while getopts ":d:r:t:T:o:h" opt; do
  case ${opt} in
    d )
      fast5=$OPTARG
      ;;
    r ) # process option t
      fastq=$OPTARG
      ;;
    t )
      threads=$OPTARG
      ;;
    T )
      tres=$OPTARG
      ;;
    o )
      out=$OPTARG
      ;;
    h )
      echo "Usage: cmd [-h] [-d <fast5.dir> -r <fastq.dir> -T <treadmill.dir> -t <threads.int> -o <output.dir>]"
      exit 1
      ;;
    \? )
      echo "Invalid option: -$OPTARG" 1>&2
      exit 1
      ;;
    : )
      echo "Invalid Option: -$OPTARG requires an argument" 1>&2
      exit 1
      ;;
      esac
done

if [ $OPTIND -eq 1 ]; then echo "Usage: cmd [-h] [-d <fast5.dir> -r <fastq.dir> -T <treadmill.dir> -t <threads.int> -o <output.dir>" && exit 1; fi

shift $((OPTIND -1))

fast5dir=$(readlink -f ${fast5})
fastqdir=$(readlink -f ${fastq})
treadmilldir=$(readlink -f ${tres})
outdir=$(readlink -f ${out})

mkdir -p ${outdir}

echo ${fast5dir}, ${fastqdir}, ${treadmilldir} ${threads}


echo "Converting multi-fast5 to single-fast5"

singlef5dir="${outdir}/single_fast5"
multi_to_single_fast5 -i ${fast5dir} -s ${singlef5dir} -t ${threads} --recursive

echo "Basecalling using Guppy (high-accuracy for modified bases)"

singlefqdir="${outdir}/single_fastq"
/opt/ont/minknow/guppy/bin/guppy_basecaller -i ${singlef5dir} -r -s ${singlefqdir} --config dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg
cat ${singlefqdir}/*.fastq > ${singlefqdir}/collapsed.fastq

echo "Running Tombo prior to modification calls"

tombo preprocess annotate_raw_with_fastqs --fast5-basedir ${singlef5dir} --fastq-filenames ${singlefqdir}/collapsed.fastq --basecall-group Basecall_1D_000 --basecall-subgroup BaseCalled_template --overwrite --processes ${threads}
tombo resquiggle ${singlef5dir} ${singlefqdir}/collapsed.fastq --processes ${threads} --corrected-group RawGenomeCorrected_001 --basecall-group Basecall_1D_000 --overwrite

echo "Calling CpG methylations using DeepSignal"

