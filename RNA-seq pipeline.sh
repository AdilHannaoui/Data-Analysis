#!/bin/bash

# This pipeline contains a Bash script designed to streamline the RNA-seq data analysis workflow.
# The script automates several key steps in the RNA-seq analysis pipeline, including quality control,
# trimming of raw sequencing reads, alignment to a reference genome, and counting gene expression levels.


SECONDS=0

# Define paths
WORKDIR="/path/to/your/working/directory"               # Working directory
DATA_DIR="$WORKDIR/data/experiment_data"                 # Directory containing the input fastq files
TRIMMO_JAR="$WORKDIR/Trimmomatic/trimmomatic.jar"       # Path to Trimmomatic jar file
HISAT2_INDEX="$WORKDIR/HISAT2/index/genome"             # Path to HISAT2 reference genome index
GTF_FILE="$WORKDIR/HISAT2/annotation.gtf"                # Path to GTF file for gene annotations
OUTPUT_DIR="$WORKDIR/output"                             # Directory for all output files
THREADS=6

# Create output directories
mkdir -p "$OUTPUT_DIR/fastqc_pre" "$OUTPUT_DIR/fastqc_post" "$OUTPUT_DIR/hisat2" "$OUTPUT_DIR/featurecounts" "$OUTPUT_DIR/logs" "$OUTPUT_DIR/fastq_trimmed"

# Change to the working directory
cd "$WORKDIR" || { echo "Failed to change directory to $WORKDIR"; exit 1; }

# Loop through all fastq files in the data directory
for FASTQ_FILE in "$DATA_DIR"/*.fastq; do
	SAMPLE_NAME=$(basename "$FASTQ_FILE" .fastq)

	echo "Processing sample: $SAMPLE_NAME"

	# Quality control before trimming
	fastqc "$FASTQ_FILE" -o "$OUTPUT_DIR/fastqc_pre/" > "$OUTPUT_DIR/logs/fastqc_pre_$SAMPLE_NAME.log" 2>&1

	# Trimmomatic
	TRIMMED_FASTQ_FILE="$OUTPUT_DIR/fastq_trimmed/${SAMPLE_NAME}_trimmed.fastq"
	java -jar "$TRIMMO_JAR" SE -threads "$THREADS" "$FASTQ_FILE" "$TRIMMED_FASTQ_FILE" \
	ILLUMINACLIP:"$WORKDIR/Trimmomatic/adapters/TruSeq3-SE.fa":2:30:10 SLIDINGWINDOW:4:20 MINLEN:20 TRAILING:10 -phred33

	# Error check for Trimmomatic
	if [ ! -s "$TRIMMED_FASTQ_FILE" ]; then
	    echo "Error: Trimmomatic did not generate the output file correctly for $SAMPLE_NAME."
	    exit 1
	fi

	# Quality control after trimming
	fastqc "$TRIMMED_FASTQ_FILE" -o "$OUTPUT_DIR/fastqc_post/" > "$OUTPUT_DIR/logs/fastqc_post_$SAMPLE_NAME.log" 2>&1

	# Alignment with HISAT2
	BAM_FILE="$OUTPUT_DIR/hisat2/${SAMPLE_NAME}_trimmed.bam"
	hisat2 -q --rna-strandness R -x "$HISAT2_INDEX" -U "$TRIMMED_FASTQ_FILE" | samtools sort -o "$BAM_FILE"

	# Error check for HISAT2
	if [ $? -ne 0 ]; then
	    echo "Error: HISAT2 or Samtools did not run correctly for $SAMPLE_NAME."
	    exit 1
	fi
	echo "HISAT2 finished for $SAMPLE_NAME!"

	# Count gene occurrences
	FEATURECOUNTS_OUTPUT="$OUTPUT_DIR/featurecounts/${SAMPLE_NAME}_featurecounts.txt"
	featureCounts -T "$THREADS" -S 2 -a "$GTF_FILE" -o "$FEATURECOUNTS_OUTPUT" "$BAM_FILE"

	# Error check for featureCounts
	if [ $? -ne 0 ]; then
	    echo "Error: featureCounts did not run correctly for $SAMPLE_NAME."
	    exit 1
	fi

	# Extract counts and save to a separate file
	cut -f1,7 "$FEATURECOUNTS_OUTPUT" | tail -n +2 > "$OUTPUT_DIR/featurecounts/Counts_${SAMPLE_NAME}.txt"
	echo "featureCounts finished for $SAMPLE_NAME!"

done

# Total execution time
duration=$SECONDS
echo "Total execution time: $(($duration / 60)) minutes and $(($duration % 60)) seconds."
