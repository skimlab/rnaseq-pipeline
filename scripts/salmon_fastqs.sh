#!/bin/bash

# pathname should be in "absolute" path
SALMON_INDEX="/Users/skim/opt/Salmon/annots/Homo_sapiens/gencode38/gencode.v38.pc_transcripts_quasi_index_salmon_1_5_2/"
SALMON_OUTPUT="outs"
FASTQ="fastq"
NUM_THREADS=8
SALMON_ARGS=""

current_env()
{
cat << EOF

[ Current Env ]
  SALMON_INDEX  : ${SALMON_INDEX}
  SALMON_OUTPUT : ${SALMON_OUTPUT}
  FASTQ         : ${FASTQ}
  NUM_THREADS   : ${NUM_THREADS}
  SALMON_ARGS   : ${SALMON_ARGS}

EOF
}

help_msg() 
{
  cmdline=$0

cat << EOF
salmon_fastq.sh [options]

[options]
  -i, --salmon_index    : a folder where salmon index files are stored
  -o, --output          : a folder where salmon outputs to be saved
  -q, --fastq           : a folder for sample folders each of which contains
                            fastq files (R1 and R2)
  -p, --num_threads     : the number of threads to use for parallel processing
                            this is passed to 'salmon'.
  --seqBias, --gcBias   : passed to 'salmon'
  -h, --help            : print this message.

EOF
}

POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -i|--salmon_index)
      SALMON_INDEX="$2"
      shift # past argument
      shift # past value
      ;;
    -o|--output)
      SALMON_OUTPUT="$2"
      shift # past argument
      shift # past value
      ;;
    -q|--fastq)
      FASTQ="$2"
      shift # past argument
      shift # past value
      ;;
    -p|--num_threads)
      NUM_THREADS="$2"
      shift # past argument
      shift # past value
      ;;
    --seqBias)
      SALMON_ARGS="${SALMON_ARGS} --seqBias"
      shift # past argument with no value
      ;;
    --gcBias)
      SALMON_ARGS="${SALMON_ARGS} --gcBias"
      shift # past argument with no value
      ;;
    --default)
      DEFAULT=YES
      shift # past argument with no value
      ;;
    -h|--help)
      echo "$(help_msg)"
      echo "$(current_env)"
      exit 1
      ;;
    -*) # unknown option
      echo "Error: Unsupported flag $1" >& 2
      echo "$(help_msg)"
      echo "$(current_env)"
      exit 1
      ;;
    *)    # unknown option
      POSITIONAL+=("$1") # save it in an array for later
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL[@]}" # restore positional parameters

if [[ -n $1 ]]; then
	echo "$(help_msg)"
	exit 1
fi

echo "$(current_env)"

for f in ${FASTQ}/*; do
  f_out=${SALMON_OUTPUT}/${f##*/}
  mkdir ${f_out}

  FASTQ_R1=${f}/*R1*.fastq.gz
  FASTQ_R2=${f}/*R2*.fastq.gz

  echo ${FASTQ_R1}
  echo ${FASTQ_R2}

  salmon quant -i ${SALMON_INDEX} \
         -l A -1 ${FASTQ_R1} -2 ${FASTQ_R2} \
         -p ${NUM_THREADS} --validateMappings -o ${f_out} ${SALMON_ARGS}
done
