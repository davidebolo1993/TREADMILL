#!/bin/bash

is_single="false" #default is multi-fast5

while getopts ":d:st:F:T:o:h" opt; do
  case ${opt} in
    d )
      fast5=$OPTARG
      ;;
    s )
      is_single='true'
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
    o )
      out=$OPTARG
      ;;
    h )
      echo "Usage: callmethylation.sh [-h] -d <fast5.in.dir> -F <treadmill.in.decoy.fa> -T <treadmill.in.reads.tsv.gz> -o <output.dir> -t <threads.int> -s <is.single.fast5.bool>"
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

if [ $OPTIND -eq 1 ]; then echo "Usage: callmethylation.sh [-h] -d <fast5.in.dir> -F <treadmill.in.decoy.fa> -T <treadmill.in.reads.tsv.gz> -o <output.dir> -t <threads.int> -s <is.single.fast5.bool>" && exit 1; fi

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


#CHECK ONT_FAST5_API INSTALLATION (MULTI-READ/SUBSET)


multi_to_single_fast5 --help > /dev/null 2>&1

if [ $? -eq 0 ]; then

  echo "Found existing ont_fast5_api installation"

else

  echo "Can't execute ont_fast5_api" && exit 1

fi


#CHECK IF GDOWN EXISTS

gdown --help > /dev/null 2>&1

if [ $? -eq 0 ]; then

  echo "Found existing gdown installation"

else

  echo "Can't execute gdown" && exit 1

fi


#CHECK IF DEEPSIGNAL SCRIPT IS IN PATH

call_modification_frequency.py -h > /dev/null 2>&1

if [ $? -eq 0 ]; then

  echo "Found existing call_modification_frequency.py script"

else

  echo "Can't execute call_modification_frequency.py" && exit 1

fi

### CHECK IF ALL INPUT FILES EXIST

#CHECK IF FAST5 DIR EXISTS

fast5dir=$(readlink -f ${fast5})

if [ ! -d ${fast5dir} ]; then #if this does not exist, skip
  
  echo "Specified FAST5 directory not found"
  exit 1

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


outdir=$(readlink -f ${out})
mkdir -p ${outdir}


#download Guppy basecaller if this does not exist

if [ ! -f ont-guppy/bin/guppy_basecaller ]; then

  wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_4.4.2_linux64.tar.gz
  tar -xzf ont-guppy_4.4.2_linux64.tar.gz
  rm ont-guppy_4.4.2_linux64.tar.gz

fi

basecaller=$(readlink -f ont-guppy/bin/guppy_basecaller)

#download model for DeepSignal if this does not exist

if [ ! -f model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+/bn_17.sn_360.epoch_9.ckpt.index ]; then

  gdown https://drive.google.com/uc?id=1meh07c9TsdIWelVTNW2M7ZKU6F-C4Muw
  tar -xzf model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+.tar.gz
  rm model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+.tar.gz

fi

dmodel=$(readlink -f model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+/bn_17.sn_360.epoch_9.ckpt)

### IF MULTI-READ FAST5s, CONVERT TO SINGLE

if [ ${is_single} == "false" ]; then 

  echo "Converting multi-fast5 to single-fast5"
  singlef5dir="${outdir}/single_fast5"
  multi_to_single_fast5 -i ${fast5dir} -s ${singlef5dir} -t ${threads} --recursive

else

  echo "Skipping conversion of multi-fast5 to single-fast5"
  singlef5dir=${fast5dir}

fi

echo "Guppy basecalling with high-accuracy"
singlefqdir="${outdir}/single_fastq"
${basecaller} -i ${singlef5dir} -r -s ${singlefqdir} --config dna_r9.4.1_450bps_hac.cfg --cpu_threads_per_caller ${threads}
#this works using GPU if correctly configured (tested)
#${basecaller} -i ${singlef5dir} -r -s ${singlefqdir} --config dna_r9.4.1_450bps_hac.cfg --gpu_runners_per_device ${threads} -x cuda:all


cat ${singlefqdir}/*.fastq > ${singlefqdir}/collapsed.fastq

echo "Filtering out undesired FASTQ"

zcat ${alleletsv} | cut -f 2 | sort | uniq > wanted.txt
grep -f wanted.txt -A 3 ${singlefqdir}/collapsed.fastq | grep -v "^--" > ${singlefqdir}/filtered.fastq
rm wanted.txt && rm ${singlefqdir}/collapsed.fastq

echo "Pre-processing with Tombo"

tombo preprocess annotate_raw_with_fastqs --fast5-basedir ${singlef5dir} --fastq-filenames ${singlefqdir}/filtered.fastq --basecall-group Basecall_1D_000 --basecall-subgroup BaseCalled_template --overwrite --processes ${threads}
tombo resquiggle ${singlef5dir} ${decoyfa} --processes ${threads} --corrected-group RawGenomeCorrected_001 --basecall-group Basecall_1D_000 --overwrite

echo "Calling CpG methylation using DeepSignal"

#cpu
deepsignal call_mods --input_path ${singlef5dir} --model_path ${dmodel} --result_file ${outdir}/deepsignal.modcalls.tsv --reference_path ${decoyfa} --corrected_group RawGenomeCorrected_001 --nproc ${threads} --is_gpu no
#modify the above command coherently to get calls using GPU instead of CPU. This requires tensorflow-gpu installation as from DeepSignal instructions
#deepsignal call_mods --input_path ${singlef5dir} --model_path ${dmodel} --result_file ${outdir}/deepsignal.modcalls.tsv --reference_path ${decoyfa} --corrected_group RawGenomeCorrected_001 --nproc ${threads} --is_gpu yes

echo -e "chrom\tpos\tstrand\tpos_in_strand\tprob_0_sum\tprob_1_sum\tcount_modified\tcount_unmodified\tcoverage\tmodification_frequency\tk_mer\tallele_name" > ${outdir}/deepsignal.modfreqs.tsv

echo "Extracting allele-specific methylation frequencies using TREADMILL output files"
  
zcat ${alleletsv} | awk '{print >> $1".treadmill.txt"; close($1".treadmill.txt")}'
names=$(ls *.treadmill.txt)

for txt in ${names}; do

  awk '{print >> "reads"$3".txt"; close("reads"$3".txt")}' ${txt}

  groups=$(ls reads*.txt)

  for reads in ${groups}; do

    cut -f 2 ${reads} | sort | uniq > names.txt
    grep -f names.txt ${outdir}/deepsignal.modcalls.tsv > regioncalls.tmp.tsv && rm names.txt
    call_modification_frequency.py --input_path regioncalls.tmp.tsv --result_file regionfreqs.tmp.tsv --prob_cf 0 && rm regioncalls.tmp.tsv
    cat regionfreqs.tmp.tsv | awk -v var=${reads} 'FS=OFS="\t"''{print $0, var}' >> ${outdir}/deepsignal.modfreqs.tsv && rm regionfreqs.tmp.tsv
    rm ${reads}

  done

  rm ${txt}

done
