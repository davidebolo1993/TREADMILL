#!/bin/bash

while getopts ":d:r:t:T:o:m:s:h" opt; do
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
    m )
      model=$OPTARG
      ;;
    s )
      script=$OPTARG
      ;;
    h )
      echo "Usage: cmd [-h] [-d <fast5.in.dir> -r <fastq.out.dir> -T <treadmill.in.dir> -t <threads.int> -o <output.dir> -m <deepsignal.model.cpkt> -s <deepsignal.frequency.py>]"
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

if [ $OPTIND -eq 1 ]; then echo "Usage: cmd [-h] [-d <fast5.in.dir> -r <fastq.out.dir> -T <treadmill.in.dir> -t <threads.int> -o <output.dir> -m <deepsignal.model.cpkt> -s <deepsignal.frequency.py>" && exit 1; fi

shift $((OPTIND -1))

fast5dir=$(readlink -f ${fast5})
fastqdir=$(readlink -f ${fastq})
treadmilldir=$(readlink -f ${tres})
outdir=$(readlink -f ${out})
dmodel=$(readlink -f ${model})
scriptmeth=$(readlink -f ${script})

mkdir -p ${outdir}

echo "Converting multi-fast5 to single-fast5"

singlef5dir="${outdir}/single_fast5"
multi_to_single_fast5 -i ${fast5dir} -s ${singlef5dir} -t ${threads} --recursive

echo "Guppy basecalling"

singlefqdir="${outdir}/single_fastq"
/opt/ont/minknow/guppy/bin/guppy_basecaller -i ${singlef5dir} -r -s ${singlefqdir} --config dna_r9.4.1_450bps_hac.cfg --cpu_threads_per_caller ${threads} 
cat ${singlefqdir}/*.fastq > ${singlefqdir}/collapsed.fastq

echo "Filtering out undesired FASTQ"

zcat ${treadmilldir}/TREADMILL.tsv.gz | cut -f 2 | sort | uniq > wanted.txt
grep -f wanted.txt -A 3 ${singlefqdir}/collapsed.fastq | grep -v "^--" > ${singlefqdir}/filtered.fastq
rm wanted.txt

echo "Pre-processing with Tombo"

tombo preprocess annotate_raw_with_fastqs --fast5-basedir ${singlef5dir} --fastq-filenames ${singlefqdir}/filtered.fastq --basecall-group Basecall_1D_000 --basecall-subgroup BaseCalled_template --overwrite --processes ${threads}
tombo resquiggle ${singlef5dir} ${treadmilldir}/fake.fa --processes ${threads} --corrected-group RawGenomeCorrected_001 --basecall-group Basecall_1D_000 --overwrite

echo "Calling allele-specific methylations with DeepSignal"

deepsignal call_mods --input_path ${singlef5dir} --model_path ${dmodel} --result_file ${outdir}/deepsignal.modcalls.tsv --reference_path ${treadmilldir}/fake.fa --corrected_group RawGenomeCorrected_001 --nproc ${threads} --is_gpu no

echo -e "chrom\tpos\tstrand\tpos_in_strand\tprob_0_sum\tprob_1_sum\tcount_modified\tcount_unmodified\tcoverage\tmodification_frequency\tk_mer" > ${outdir}/a1.deepsignal.modfreqs.tsv
echo -e "chrom\tpos\tstrand\tpos_in_strand\tprob_0_sum\tprob_1_sum\tcount_modified\tcount_unmodified\tcoverage\tmodification_frequency\tk_mer" > ${outdir}/a2.deepsignal.modfreqs.tsv

zcat ${treadmilldir}/TREADMILL.tsv.gz | awk '{print >> $1".treadmill.txt"; close($1".treadmill.txt")}'
names=$(ls *.treadmill.txt)

for txt in ${names}; do

  awk '{print >> "reads"$3".txt"; close("reads"$3".txt")}' ${txt}
  cut -f 2 reads1.txt | sort | uniq > names1.txt && rm reads1.txt
  cut -f 2 reads2.txt | sort | uniq > names2.txt && rm reads2.txt

  grep -f names1.txt ${outdir}/deepsignal.modcalls.tsv > regioncalls.tmp.tsv
  python ${scriptmeth} --input_path regioncalls.tmp.tsv --result_file regionfreqs.tmp.tsv --prob_cf 0 && rm regioncalls.tmp.tsv
  cat regionfreqs.tmp.tsv >> ${outdir}/a1.deepsignal.modfreqs.tsv

  isdiff=$(cmp names1.txt names2.txt)

  rm names1.txt

  if [ -z ${isdiff} ]; then #if is empty, then a single allele and we do not have to split

    cat regionfreqs.tmp.tsv >> ${outdir}/a2.deepsignal.modfreqs.tsv && rm regionfreqs.tmp.tsv && rm names2.txt

  else 

    rm regionfreqs.tmp.tsv # the one from the other allele has already been stored properly
    grep -f names2.txt ${outdir}/deepsignal.modcalls.tsv > regioncalls.tmp.tsv && rm names2.txt
    python ${scriptmeth} --input_path regioncalls.tmp.tsv --result_file regionfreqs.tmp.tsv --prob_cf 0 && rm regioncalls.tmp.tsv
    cat regionfreqs.tmp.tsv >> ${outdir}/a2.deepsignal.modfreqs.tsv && rm regionfreqs.tmp.tsv

  fi

  rm ${txt}

done
