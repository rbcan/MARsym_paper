#!/bin/bash

source ~/.bashrc

DATE=`date +%d%m%Y`
logfile=$JOB_ID.haplotyping.$DATE.log
exec > $logfile 2>&1

# specify all of the following "variables"
ref="" # reference assembly: *.fa
gene="" # region for gene to be extracted from reference e.g. NODE_1:6932-7552
bam="" # sample name in *.bam
length="" # gene length in number of nucleotides, e.g. 1500

mkdir ViQuaS
cd ViQuaS # copy reference assembly to working directory
mkdir fasta
mkdir bam

samtools faidx ./$ref.fa $gene > ./fasta/$gene.fa
samtools view -b ./bam/$bam.bam $gene > ./bam/$gene.bam
samtools index ./bam/$gene.bam

samtools view -h -F 4 -b ./bam/$gene.bam > ./bam/$gene.onlymapped.bam

timeout 5m Rscript /home/ansorge/programs/ViQuaS1.3/ViQuaS.R ./fasta/$gene.fa ./bam/$gene.onlymapped.bam 5 0.7 1 $length
mv ViQuaS-Richness.txt Richness.$gene
mv ViQuaS-Spectrum.fa Spectrum.$gene
rm -r SIM
rm -r Samples
rm -r viquas
