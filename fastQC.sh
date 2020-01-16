#!/bin/bash

#author: Michael Skaro
# The purpose of this shell script is prepare an RNA seq raw data output file(.fq) for data analysis in DESeq2(R package). 
# The GDC pipeline will be used to analyze this data, so tha itmay be compared to other cancer
# data downloaded from TCGA
# Required packages: FastQC, Trimmomatic, Star, htseq
# data required: gdc.v22.gtf, star 2 pass indicies(downloadable for humans)

# This will be piece mealed into 56 commands as I have no access to sapelo2 right now. Once I have been migrated properly
# i will iterate instead of the individual commands. 


#date: 1/14/2020

fastqc -f fastq -o QC_out ERR031017_1.fastq;
fastqc -f fastq -o QC_out ERR031017_2.fastq;
fastqc -f fastq -o QC_out ERR031018_1.fastq;
fastqc -f fastq -o QC_out ERR031018_2.fastq;
fastqc -f fastq -o QC_out ERR031019_1.fastq;
fastqc -f fastq -o QC_out ERR031019_2.fastq;
fastqc -f fastq -o QC_out ERR031022_1.fastq;
fastqc -f fastq -o QC_out ERR031022_2.fastq;
fastqc -f fastq -o QC_out ERR031023_1.fastq;
fastqc -f fastq -o QC_out ERR031023_2.fastq;
fastqc -f fastq -o QC_out ERR031024_1.fastq;
fastqc -f fastq -o QC_out ERR031024_2.fastq;
fastqc -f fastq -o QC_out ERR031025_1.fastq;
fastqc -f fastq -o QC_out ERR031025_2.fastq;
fastqc -f fastq -o QC_out ERR031026_1.fastq;
fastqc -f fastq -o QC_out ERR031026_2.fastq;
fastqc -f fastq -o QC_out ERR031027_1.fastq;
fastqc -f fastq -o QC_out ERR031027_2.fastq;
fastqc -f fastq -o QC_out ERR031028_1.fastq;
fastqc -f fastq -o QC_out ERR031028_2.fastq;
fastqc -f fastq -o QC_out ERR031029_1.fastq;
fastqc -f fastq -o QC_out ERR031029_2.fastq;
fastqc -f fastq -o QC_out ERR031031_1.fastq;
fastqc -f fastq -o QC_out ERR031031_2.fastq;
fastqc -f fastq -o QC_out ERR031032_1.fastq;
fastqc -f fastq -o QC_out ERR031032_2.fastq;
fastqc -f fastq -o QC_out ERR031033_1.fastq;
fastqc -f fastq -o QC_out ERR031033_2.fastq;
fastqc -f fastq -o QC_out ERR031035_1.fastq;
fastqc -f fastq -o QC_out ERR031035_2.fastq;
fastqc -f fastq -o QC_out ERR031038_1.fastq;
fastqc -f fastq -o QC_out ERR031038_2.fastq;
fastqc -f fastq -o QC_out ERR031039_1.fastq;
fastqc -f fastq -o QC_out ERR031039_2.fastq;
fastqc -f fastq -o QC_out ERR031040_1.fastq;
fastqc -f fastq -o QC_out ERR031040_2.fastq;
fastqc -f fastq -o QC_out ERR031041_1.fastq;
fastqc -f fastq -o QC_out ERR031041_2.fastq;
fastqc -f fastq -o QC_out ERR031042_1.fastq;
fastqc -f fastq -o QC_out ERR031042_2.fastq;
fastqc -f fastq -o QC_out ERR031043_1.fastq;
fastqc -f fastq -o QC_out ERR031043_2.fastq;
fastqc -f fastq -o QC_out ERR031044_1.fastq;
fastqc -f fastq -o QC_out ERR031044_2.fastq;
fastqc -f fastq -o QC_out ERR299295_1.fastq;
fastqc -f fastq -o QC_out ERR299295_2.fastq;
fastqc -f fastq -o QC_out ERR299296_1.fastq;
fastqc -f fastq -o QC_out ERR299296_2.fastq;
fastqc -f fastq -o QC_out ERR299297_1.fastq;
fastqc -f fastq -o QC_out ERR299297_2.fastq;
fastqc -f fastq -o QC_out ERR299298_1.fastq;
fastqc -f fastq -o QC_out ERR299298_2.fastq;
fastqc -f fastq -o QC_out ERR299299_1.fastq;
fastqc -f fastq -o QC_out ERR299299_2.fastq;



#END
