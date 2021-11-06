#!/bin/bash

OUTPUT_FASTQ="fastq"
INPUT_FASTQ="fastq_input"
SAMPLE_ID="sample_id"

current_env()
{
cat << EOF

[Current Env]

  OUTPUT_FASTQ  : ${OUTPUT_FASTQ}
  INPUT_FASTQ   : ${INPUT_FASTQ}
  SAMPLE_ID     : ${SAMPLE_ID}

  fastq files from :

$(ls -d ./${INPUT_FASTQ}/${SAMPLE_ID}_L*) 

  will be merged and saved into:

./${OUTPUT_FASTQ}/${SAMPLE_ID}/.

EOF
}

help_msg() 
{
  cmdline=$0

cat << EOF
merge_fastq.sh [options] sample_id

[options]
  -o, --output_fastq   : a folder to save merged fastq files.
  -i, --input_fastq    : a folder from input fastq files (fastq file folders) are saved.

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
    -h|--help)
      echo "$(help_msg)"
      echo "$(current_env)"
      exit 1
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

echo "$(current_env)"

mkdir ${OUTPUT_FASTQ}/${SAMPLE_ID}
cat ${INPUT_FASTQ}/${SAMPLE_ID}_L*/*R1*.fastq.gz > ${OUTPUT_FASTQ}/${SAMPLE_ID}/${SAMPLE_ID}_R1.fastq.gz
cat ${INPUT_FASTQ}/${SAMPLE_ID}_L*/*R2*.fastq.gz > ${OUTPUT_FASTQ}/${SAMPLE_ID}/${SAMPLE_ID}_R2.fastq.gz

ls -l ${OUTPUT_FASTQ}/${SAMPLE_ID}
