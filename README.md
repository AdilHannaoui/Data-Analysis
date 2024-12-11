---
title: "ChIP-seq workflow"
author: "Adil Hannaoui Anaaoui"
date: "December 11th, 2024"
---

## ChIP-seq
<div align="justify">
This repository compiles a workflow designed for ChIP-seq data analysis. It includes a detailed description of the techniques and tools used, explores available alternatives, and provides examples of expected results.<br><br>

Chromatin immunoprecipitation (ChIP) experiments are performed to identify DNA bound to specific (chromatin) proteins of interest. The first step involves isolating the chromatin and immunoprecipitating (IP) fragements with an antibody against the protein of interest. In ChIP-seq, the immunoprecipitated DNA fragments are then sequenced, followed by identification of enriched regions of DNA or peaks. These peak calls can then be used to make biological inferences by determining the associated genomic features and/or over-represented sequence motifs.<br><br>
</div>

![ChIP Workflow](./img/ChIP_technique.png)<br><br>


```bash
$ fastqc reads.fastq
```

```bash
$ bowtie2 -p 2
          -q
          --local 
		      -x ~/direcory/to/genome/index 
		      -U ~/reads.fastq 
		      -S ~/reads_mapped.sam
```

```bash
$ samtools view -h \
	              -S \
	              -b \
	              -o /reads_mapped.bam /reads_mapped.sam
```

```bash
$ sambamba sort -t 2 \ 
	              -o /reads_mapped.bam /reads_mapped_sorted.bam 
```

```bash
$ sambamba view -h \
	              -t 2 \
	              -f bam \ 
                -F "[XS] == null and not unmapped  and not duplicate" \ 
	              reads_mapped_sorted.bam  > reads_mapped_sorted_aligned.bam 
```

```bash
$ macs2 callpeak -t bowtie2/H1hesc_Nanog_Rep1_aln.bam \
                 -c bowtie2/H1hesc_Input_Rep1_aln.bam \
 	               -f BAM \
                 -g 12000000 \
                 -p 0.01
	               -n sample_name 
```
