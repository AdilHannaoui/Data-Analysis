#!/bin/bash

# This workflow can be used for both RIP-seq and ChIP-seq analyses, as both experiments aim to identify enriched regions through peak calling.
#The script follows essential data processing steps, including quality assessment, sequence trimming, alignment, and peak detection,
#thus enabling reproducible identification of enriched regions. By maintaining a standardized workflow, the analysis ensures comparability
#and efficiency in studies of protein-DNA or protein-RNA interactions.


# Initialize timer
SECONDS=0

# Working directories and configurations
WORKDIR="/path/to/working_directory"
FASTQ_DIR="$WORKDIR/data/samples"
TRIMMO_JAR="$WORKDIR/Trimmomatic-0.39/"
BOWTIE2_INDEX="$WORKDIR/Bowtie2/genome_index/genome"
OUTPUT_DIR="$WORKDIR/output/study"
THREADS=6
MAX_REPLICATES=3  # Set the number of replicates to process; here we assume 3 replicates per sample.

# Header file to add to the count matrix after creation
HEADER_FILE="$WORKDIR/header.txt"

# Create necessary output directories
mkdir -p "$OUTPUT_DIR/fastqc_pre" "$OUTPUT_DIR/fastqc_post" "$OUTPUT_DIR/bowtie2" \
         "$OUTPUT_DIR/logs" "$OUTPUT_DIR/macs2" "$OUTPUT_DIR/fastq_trimmed"

# Change to the working directory
cd "$WORKDIR" || exit

# Step 1: Quality assessment, trimming, and alignment
for FASTQ_FILE in "$FASTQ_DIR"/*.fastq; do
    SAMPLE_NAME=$(basename "$FASTQ_FILE" .fastq)
    
    echo "Processing sample: $SAMPLE_NAME"

    # Initial quality check
    fastqc "$FASTQ_FILE" -o "$OUTPUT_DIR/fastqc_pre/" > "$OUTPUT_DIR/logs/fastqc_pre_$SAMPLE_NAME.log" 2>&1

    # Trimming with Trimmomatic (adjust parameters as needed)
    TRIMMED_FASTQ_FILE="$OUTPUT_DIR/fastq_trimmed/${SAMPLE_NAME}.fastq"
    if [ ! -s "$TRIMMED_FASTQ_FILE" ]; then
        java -jar "$TRIMMO_JAR/trimmomatic-0.39.jar" SE -threads "$THREADS" -phred33 "$FASTQ_FILE" \
            "$TRIMMED_FASTQ_FILE" ILLUMINACLIP:$TRIMMO_JAR/adapters/TruSeq3-SE.fa:2:30:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:20
        echo "Trimmomatic completed for sample $SAMPLE_NAME!"

        # Error check for Trimmomatic
        if [ ! -s "$TRIMMED_FASTQ_FILE" ]; then
            echo "Error: Trimmomatic did not generate the output file correctly for $SAMPLE_NAME."
            exit 1
        fi
    fi

    # Post-trimming quality check
    fastqc "$TRIMMED_FASTQ_FILE" -o "$OUTPUT_DIR/fastqc_post/" > "$OUTPUT_DIR/logs/fastqc_post_$SAMPLE_NAME.log" 2>&1

    # Alignment with Bowtie2 and BAM file generation
    BAM_FILE="$OUTPUT_DIR/bowtie2/${SAMPLE_NAME}.bam"
    if [ ! -s "$BAM_FILE" ]; then
        bowtie2 -p "$THREADS" -x "$BOWTIE2_INDEX" -U "$TRIMMED_FASTQ_FILE" | samtools sort -o "$BAM_FILE"
        samtools index "$BAM_FILE"
    fi
done

# Step 2: Peak calling with MACS2
for REP in $(seq 1 $MAX_REPLICATES); do
    for IP_FILE in "$OUTPUT_DIR/bowtie2"/*_IP${REP}.bam; do
        BASE_NAME=$(basename "$IP_FILE" | sed "s/_IP${REP}.*//")
        IN_FILE="${OUTPUT_DIR}/bowtie2/${BASE_NAME}_IN${REP}.bam"

        # Check for the existence of IN and IP files
        if [[ -f "$IN_FILE" && -f "$IP_FILE" ]]; then
            echo "Processing sample: $BASE_NAME, replicate: $REP"
            echo "Using IP file: $IP_FILE"
            echo "Using IN file: $IN_FILE"

            # Run MACS2 for peak calling
            macs2 callpeak -t "$IP_FILE" -c "$IN_FILE" --format BAM \
                --name "$OUTPUT_DIR/macs2/${BASE_NAME}_rep${REP}" --gsize 1.2e7 --pvalue 0.01 --nomodel > \
                "$OUTPUT_DIR/logs/narrowpeaks_${BASE_NAME}.log" 2>&1
        else
            echo "IN file for $IP_FILE for replicate $REP not found"
        fi
    done
done

# Step 3: Peak intersection and count generation
for FILE in "$OUTPUT_DIR/macs2"/*_rep1_peaks.narrowPeak; do
    BASE_NAME=$(basename "$FILE" | sed 's/_rep1_peaks.narrowPeak//')
    NARROWPEAK_FILES=()

    # Collect narrowPeak files for replicates
    for REP in $(seq 1 $MAX_REPLICATES); do
        NARROWPEAK_FILE="$OUTPUT_DIR/macs2/${BASE_NAME}_rep${REP}_peaks.narrowPeak"
        if [[ -f "$NARROWPEAK_FILE" ]]; then
            NARROWPEAK_FILES+=("$NARROWPEAK_FILE")
        fi
    done

    # Check if there are enough narrowPeak files for intersection
    if [[ ${#NARROWPEAK_FILES[@]} -ge 2 ]]; then
        echo "Intersecting peaks for sample $BASE_NAME"

        # Intersect and remove duplicates
        bedtools intersect -a "${NARROWPEAK_FILES[0]}" -b "${NARROWPEAK_FILES[@]:1}" -wa | sort -u > \
            "$OUTPUT_DIR/macs2/${BASE_NAME}_peaks.bed"

        BAM_FILES=()
        for REP in $(seq 1 $MAX_REPLICATES); do
            BAM_FILES+=("$OUTPUT_DIR/bowtie2/${BASE_NAME}_IP${REP}.bam")
            BAM_FILES+=("$OUTPUT_DIR/bowtie2/${BASE_NAME}_IN${REP}.bam")
        done

        # Generate count matrix
        bedtools multicov -bams "${BAM_FILES[@]}" -bed "$OUTPUT_DIR/macs2/${BASE_NAME}_peaks.bed" > \
            "$OUTPUT_DIR/macs2/${BASE_NAME}_counts.txt"

        # Add header to count matrix
        if [[ -f "$HEADER_FILE" ]]; then
            cat "$HEADER_FILE" "$OUTPUT_DIR/macs2/${BASE_NAME}_counts.txt" > \
                "$OUTPUT_DIR/macs2/${BASE_NAME}_counts_with_header.txt"
            rm -rf "$OUTPUT_DIR/macs2/${BASE_NAME}_counts.txt"
        else
            echo "Header file $HEADER_FILE not found. Cannot add header."
        fi

        # Display BAM files used
        echo "BAM files used for bedtools multicov: ${BAM_FILES[*]}"
    else
        echo "Not enough narrowPeak files to intersect for sample $BASE_NAME"
    fi
done

# Total execution time
duration=$SECONDS
echo "Total execution time: $(($duration / 60)) minutes and $(($duration % 60)) seconds."
