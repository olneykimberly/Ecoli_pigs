#!/bin/bash

# change directory
cd ../fastq

# create file with list of R1 samples
ls | grep 'L003_R1_001.fastq.gz' > R1Samples_snRNA.txt

# loops through list and print first line
touch sampleReadInfo_snRNA.txt
for sample in `cat R1Samples_snRNA.txt`; do
    #printf "${sample}\t"
    zcat ${sample} | head -1 >> sampleReadInfo_snRNA.txt	
done;

paste -d "\t" R1Samples_snRNA.txt sampleReadInfo_snRNA.txt > sampleReadGroupInfo_snRNA.txt
rm R1Samples_snRNA.txt
rm sampleReadInfo_snRNA.txt
mv sampleReadGroupInfo_snRNA.txt ../scripts/sampleReadGroupInfo_snRNA.txt
