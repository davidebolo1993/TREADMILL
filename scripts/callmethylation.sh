#!/bin/bash

is_single="false" #default is multi-fast5

while getopts ":d:sr:t:F:T:m:p:o:h" opt; do
  case ${opt} in
    d )
      fast5=$OPTARG
      ;;
    s )
      is_single='true'
      ;;
    r ) 
      fastq=$OPTARG
      ;;
    t )
      threads=$OPTARG
      ;;
    F )
      tresfa=$OPTARG
      ;;
    T )
      trestsv=$OPTARG
      ;;
    m )
      model=$OPTARG
      ;;
    p )
      script=$OPTARG
      ;;
    o )
      out=$OPTARG
      ;;
    h )
      echo "Usage: cmd [-h] -d <fast5.in.dir> -r <fastq.in.dir> -F <treadmill.in.decoy.fa> -T <treadmill.in.reads.tsv> -o <output.dir> -t <threads.int> -m <deepsignal.model.cpkt> -p <deepsignal.callfrequency.py> -s"
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

if [ $OPTIND -eq 1 ]; then echo "Usage: cmd [-h] -d <fast5.in.dir> -r <fastq.out.dir/fastq.in.file> -F <treadmill.in.decoy.fa> -T <treadmill.in.reads.tsv> -o <output.dir> -t <threads.int> -m <deepsignal.model.cpkt> -p <deepsignal.callfrequency.py> -s" && exit 1; fi

shift $((OPTIND -1))

### CHECK THAT ALL THE COMMANDS CAN BE RUN PROPERLY BEFORE STARTING

#CHECK DEEPSIGNAL INSTALLATION

deepsignal --help > /dev/null 2>&1

if [ $? -eq 0 ]; then
   
   echo "Found existing DeepSignal installation"

else
   
   echo "Can't execute DeepSignal" && exit 1

fi


#CHECK TOMBO INSTALLATION

tombo --help > /dev/null 2>&1

if [ $? -eq 0 ]; then
   
   echo "Found existing Tombo installation"

else
   
   echo "Can't execute Tombo" && exit 1

fi


#CHECK ONT_FAST5_API INSTALLATION IF IS MULTI-READ FAST5


if [ ${is_single} == "false" ]; then

  multi_to_single_fast5 --help > /dev/null 2>&1

  if [ $? -eq 0 ]; then

    echo "Found existing multi_to_single_fast5 installation"

  else

    echo "Can't execute multi_to_single_fast5" && exit 1

  fi

fi

### CHECK IF ALL INPUT FILES EXIST

#CHECK IF FAST5 DIR EXISTS

fast5dir=$(readlink -f ${fast5})

if [ ! -d ${fast5dir} ]; then #if this does not exist, skip
  
  echo "Specified FAST5 directory not found"
  exit 1

fi

#CHECK IF FASTQ DIR EXISTS

fastqdir=$(readlink -f ${fastq})

if [ ${is_single} == "false" ]; then #create fastq directory

  mkdir -p ${fastqdir}

fi

#CHECK IF DECOY FASTA EXIST

decoyfa=$(readlink -f ${tresfa})

if [ ! -f ${decoyfa} ]; then 

  echo "Specified decoy FASTA not found"
  exit 1

fi

#CHECK IF TSV FILE EXIST

alleletsv=$(readlink -f ${trestsv})

if [ ! -f ${alleletsv} ]; then 

  echo "Specified read names in TSV.GZ not found"
  exit 1

fi

dmodel=$(readlink -f ${model})

if [ ! -f ${dmodel}".index" ]; then 

  echo "Specified model file for DeepSignal not found"
  exit 1

fi

scriptmeth=$(readlink -f ${script})

if [ ! -f ${scriptmeth} ]; then 

  echo "Specified call_modification_frequency.py script from DeepSignal not found"
  exit 1

fi

outdir=$(readlink -f ${out})
mkdir -p ${outdir}


### IF MULTI-READ FAST5 FILES, CONVERT TO SINGLE

if [ ${is_single} == "false" ]; then 

  echo "Converting multi-fast5 to single-fast5"
  singlef5dir="${outdir}/single_fast5"
  multi_to_single_fast5 -i ${fast5dir} -s ${singlef5dir} -t ${threads} --recursive
  echo "Guppy basecalling"
  singlefqdir="${outdir}/single_fastq"
  /opt/ont/minknow/guppy/bin/guppy_basecaller -i ${singlef5dir} -r -s ${singlefqdir} --config dna_r9.4.1_450bps_hac.cfg --cpu_threads_per_caller ${threads}

else

  echo "Skipping conversion of multi-fast5 to single-fast5"
  singlef5dir=${fast5dir}
  singlefqdir=${fastqdir}

fi

cat ${singlefqdir}/*.fastq > ${singlefqdir}/collapsed.fastq

echo "Filtering out undesired FASTQ"

zcat ${alleletsv} | cut -f 2 | sort | uniq > wanted.txt
grep -f wanted.txt -A 3 ${singlefqdir}/collapsed.fastq | grep -v "^--" > ${singlefqdir}/filtered.fastq
rm wanted.txt && rm collapsed.fastq


echo "Pre-processing with Tombo"

tombo preprocess annotate_raw_with_fastqs --fast5-basedir ${singlef5dir} --fastq-filenames ${singlefqdir}/filtered.fastq --basecall-group Basecall_1D_000 --basecall-subgroup BaseCalled_template --overwrite --processes ${threads}
tombo resquiggle ${singlef5dir} ${decoyfa} --processes ${threads} --corrected-group RawGenomeCorrected_001 --basecall-group Basecall_1D_000 --overwrite


echo "Calling allele-specific methylations with DeepSignal"

deepsignal call_mods --input_path ${singlef5dir} --model_path ${dmodel} --result_file ${outdir}/deepsignal.modcalls.tsv --reference_path ${decoyfa} --corrected_group RawGenomeCorrected_001 --nproc ${threads} --is_gpu no
echo -e "chrom\tpos\tstrand\tpos_in_strand\tprob_0_sum\tprob_1_sum\tcount_modified\tcount_unmodified\tcoverage\tmodification_frequency\tk_mer\tallele_name" > ${outdir}/deepsignal.modfreqs.tsv

zcat ${alleletsv} | awk '{print >> $1".treadmill.txt"; close($1".treadmill.txt")}'
names=$(ls *.treadmill.txt)

for txt in ${names}; do

  awk '{print >> "reads"$3".txt"; close("reads"$3".txt")}' ${txt}

  groups=$(ls reads*.txt)

  for reads in ${groups}; do

    cut -f 2 ${reads} | sort | uniq > names.txt
    grep -f names.txt ${outdir}/deepsignal.modcalls.tsv > regioncalls.tmp.tsv && rm names.txt
    python ${scriptmeth} --input_path regioncalls.tmp.tsv --result_file regionfreqs.tmp.tsv --prob_cf 0 && rm regioncalls.tmp.tsv
    cat regionfreqs.tmp.tsv | awk -v var=${reads} 'FS=OFS="\t"''{print $0, var}' >> ${outdir}/deepsignal.modfreqs.tsv
    rm ${reads}

  done

  rm ${txt}

done