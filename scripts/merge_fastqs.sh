#!/bin/bash

OUTPUT_FASTQ="fastq"
INPUT_FASTQ="fastq_input"

help_msg() 
{
  cmdline=$0

cat << EOF
merge_fastq.sh [options] sample_id

[options]
  -o, --output_fastq   : a folder to save merged fastq files.
  -i, --input_fastq    : a folder to save merged fastq files.

EOF
}

POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -o|--output_fastq)
      OUTPUT_FASTQ="$2"
      shift # past argument
      shift # past value
      ;;
    -i|--input_fastq)
      INPUT_FASTQ="$2"
      shift # past argument
      shift # past value
      ;;
    --default)
      DEFAULT=YES
      shift # past argument with no value
      ;;
    -*|--*)
      # unknown option
      echo "Error: Unsupported flag $1" >& 2
      exit 1
      ;;
    *)    # unknown option
      POSITIONAL+=("$1") # save it in an array for later
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL[@]}" # restore positional parameters

if [[ -z $1 ]]; then
	echo "$(help_msg)"
	exit 1
else
	SAMPLE_ID=$1
fi

echo "OUTPUT_FASTQ  : $OUTPUT_FASTQ"
echo "INPUT_FASTQ   : $INPUT_FASTQ"
echo "SAMPLE_ID     : $SAMPLE_ID"

mkdir ${OUTPUT_FASTQ}/${SAMPLE_ID}
cat ${INPUT_FASTQ}/${SAMPLE_ID}*/*R1*.fastq.gz > ${OUTPUT_FASTQ}/${SAMPLE_ID}/${SAMPLE_ID}_R1.fastq.gz
cat ${INPUT_FASTQ}/${SAMPLE_ID}*/*R2*.fastq.gz > ${OUTPUT_FASTQ}/${SAMPLE_ID}/${SAMPLE_ID}_R2.fastq.gz
