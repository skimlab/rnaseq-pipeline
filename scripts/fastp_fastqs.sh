#!/bin/bash

# pathname should be in "absolute" path
FASTP_OUTPUT="fastp"
FASTQ="fastq"
WRITE_FASTQ="NO"

current_env()
{
  echo "[Current Env]"

  echo "  FASTP_OUTPUT : ${FASTP_OUTPUT}"
  echo "  FASTQ        : ${FASTQ}"
  echo "  WRITE_FASTQ  : ${WRITE_FASTQ} "
}

help_msg() 
{
  cmdline=$0

cat << EOF
fastp_fastq.sh [options]

[options]
  -o, --output          : a folder where fastp outputs to be saved
  -q, --fastq           : a folder for sample folders each of which contains
                            fastq files (R1 and R2)
  -w, --write           : write filtered fastq files
  -h, --help            : print this message.

EOF
}

POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -o|--output)
      FASTP_OUTPUT="$2"
      shift # past argument
      shift # past value
      ;;
    -q|--fastq)
      FASTQ="$2"
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
    -w|--write)
      WRITE_FASTQ="YES"
      shift # past argument with no value
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
  sample_id=${f##*/}
  f_out=${FASTP_OUTPUT}/${sample_id}

  echo "${f_out}"
  mkdir ${f_out}

  FASTQ_R1=${f}/${sample_id}_R1.fastq.gz
  FASTQ_R2=${f}/${sample_id}_R2.fastq.gz

  FASTQ_R1_TRIMMED=${f}/${sample_id}_R1_trimmed.fastq.gz
  FASTQ_R2_TRIMMED=${f}/${sample_id}_R2_trimmed.fastq.gz

  JSON=${f_out}/fastp.json
  HTML=${f_out}/fastp.html

  # echo ${FASTQ_R1}
  # echo ${FASTQ_R2}

  if [[ ${WRITE_FASTQ} == "YES" ]]; then
    fastp -j ${JSON} -h ${HTML} -i ${FASTQ_R1} -I ${FASTQ_R2} -o ${FASTQ_R1_TRIMMED} -O ${FASTQ_R2_TRIMMED}
  else
    fastp -j ${JSON} -h ${HTML} -i ${FASTQ_R1} -I ${FASTQ_R2}
  fi

done
