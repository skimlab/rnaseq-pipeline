#!/bin/bash

# pathname should be in "absolute" path
FASTQC_OUTPUT="fastqc"
FASTQ="fastq"
NUM_THREADS=2

current_env()
{
  echo "[Current Env]"

  echo "  FASTQC_OUTPUT : ${FASTQC_OUTPUT}"
  echo "  FASTQ         : ${FASTQ}"
  echo "  NUM_THREADS   : ${NUM_THREADS}"
}

help_msg() 
{
  cmdline=$0

cat << EOF
fastqc_fastq.sh [options]

[options]
  -o, --output          : a folder where fastqc outputs to be saved
  -q, --fastq           : a folder for sample folders each of which contains
                            fastq files (R1 and R2)
  -p, --num_threads     : the number of threads to use for parallel processing
                            this is passed to 'fastqc'.
  -h, --help            : print this message.

EOF
}

POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -o|--output)
      FASTQC_OUTPUT="$2"
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
  f_out=${FASTQC_OUTPUT}/${f##*/}

  echo "${f_out}"

  mkdir ${f_out}

  fastqc -t ${NUM_THREADS} -o ${f_out} ${f}/*.fastq.gz 
done
