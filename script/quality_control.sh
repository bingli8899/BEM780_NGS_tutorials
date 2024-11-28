#!/bin/bash

INPUT_DIR=$1
THREADS=$2
OUTPUT_DIR=$3

CLEANED_DIR="${OUTPUT_DIR}/cleaned"
FASTQC_DIR="${OUTPUT_DIR}/fastqc"
MULTIQC_DIR="${OUTPUT_DIR}/multiqc"

mkdir -p "$CLEANED_DIR" "$FASTQC_DIR" "$MULTIQC_DIR"

echo "Running fastp..."
for FILE in "${INPUT_DIR}"/*.{fastq,fastq.gz}; do # Run fastp for each FASTQ or FASTQ.GZ file 
    if [[ -f "$FILE" ]]; then
        BASENAME=$(basename "$FILE")
        BASENAME_NO_EXT="${BASENAME%.fastq}" # Remove .fastq or .fastq.gz
        BASENAME_NO_EXT="${BASENAME_NO_EXT%.gz}" # Handle .gz if present
        fastp -i "$FILE" -o "${CLEANED_DIR}/${BASENAME_NO_EXT}_cleaned.fastq" --thread "$THREADS"
    fi
done

if [[ -f fastp.html ]]; then
    mv fastp.html "${OUTPUT_DIR}/fastp.html"
fi

if [[ -f fastp.json ]]; then
    mv fastp.json "${OUTPUT_DIR}/fastp.json"
fi


fastqc "${CLEANED_DIR}"/*.fastq -t "$THREADS" -o "$FASTQC_DIR" # Run FastQC on cleaned files 
multiqc "$FASTQC_DIR" -o "$MULTIQC_DIR" # Run MultiQC to aggregate FastQC reports 

echo "Oyay! Outputs are in $OUTPUT_DIR."
