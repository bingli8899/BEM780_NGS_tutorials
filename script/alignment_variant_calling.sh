#!/bin/bash
# The below code runs the alignment and variant calling pipeline. 

INPUT_FASTA="$1" # Input FASTA directory -- Need the actual input directory 
REF_GENOME="$2"  # Reference genome file -- Need the file name 
THREADS="$3"  # Number of threads
OUTPUT_DIR="$4"  # Output directory

mkdir -p "$OUTPUT_DIR" 

bowtie2-build "$REF_GENOME" "${OUTPUT_DIR}/reference_index" # Index the refernce 

for INPUT_FILE in "${INPUT_DIR}"/*.{fasta,fa,fastq,fastq.gz}; do

    if [[ ! -e "$INPUT_FILE" ]]; then
        continue
    fi

    FILE_TYPE=$(echo "$INPUT_FILE" | awk -F . '{if (NF>1) print $NF}')
    echo "Processing file: $INPUT_FILE"

    if [[ "$FILE_TYPE" == "fasta" || "$FILE_TYPE" == "fa" ]]; then # Handle fasta, fastq, and fastq.gz files 
        echo "Detected FASTA file."
        TEMP_FASTQ="${INPUT_FILE%.*}.fastq"
        seqtk seq -A "$INPUT_FILE" > "$TEMP_FASTQ"
        INPUT_FILE="$TEMP_FASTQ"
    elif [[ "$FILE_TYPE" == "gz" ]]; then
        echo "Detected FASTQ.GZ file."
    elif [[ "$FILE_TYPE" == "fastq" ]]; then
        echo "Detected FASTQ file."
    else
        echo "Unsupported file type: $FILE_TYPE. Skipping."
        continue
    fi

    BASENAME=$(basename "$INPUT_FILE" | sed 's/\..*//')

    bowtie2 -x "${OUTPUT_DIR}/reference_index" \
        -U "$INPUT_FILE" \
        -S "${OUTPUT_DIR}/${BASENAME}.sam" \
        -p "$THREADS"

    if [[ -f "$TEMP_FASTQ" ]]; then
        echo "Cleaning up temporary files..."
        rm "$TEMP_FASTQ"
    fi
    
done

samtools view -Sb "${OUTPUT_DIR}/${BASENAME}.sam" > "${OUTPUT_DIR}/${BASENAME}.bam"  # convert sam file to bam files 
samtools sort "${OUTPUT_DIR}/${BASENAME}.bam" -o "${OUTPUT_DIR}/${BASENAME}_sorted.bam" # sort bam files 
samtools index "${OUTPUT_DIR}/${BASENAME}_sorted.bam" # Index bam files 

# Variant calling with bcftools 
bcftools mpileup -f "$REF_GENOME" "${OUTPUT_DIR}/${BASENAME}_sorted.bam" | \
    bcftools call -mv -Oz -o "${OUTPUT_DIR}/${BASENAME}_variants.vcf.gz"
bcftools index "${OUTPUT_DIR}/${BASENAME}_variants.vcf.gz"

echo "Nice! All output stored to $OUTPUT_DIR."
