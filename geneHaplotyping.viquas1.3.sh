#!/bin/bash

source ~/.bashrc

DATE=`date +%d%m%Y`
logfile=geneHaplotyping.$DATE.log
exec > $logfile 2>&1

echo =================================================================
echo "Script started $(date)"
DATE=`date +%d%m%Y`
echo "Your logfile is geneHaplotyping."$DATE".log" 
echo =================================================================

cd /PATH/ViQuaS1.3/ # this folder contains the ViQuaS program parts - important to run script in that folder

mkdir ./fasta.files/
mkdir ./bam.files/

# preparation file sample.viquasIn for all samples s1, s2, s3, ... sn: Contig <tab> geneStart <tab> geneEnd; e.g. Node_15  500 1500
# this example runs with sample s1

gene_pos=$(awk -F "\t" '{print $2":"$3"-"$4}' s1.viquasIn | tr '\n' '\t')

for gene in $gene_pos
do
samtools faidx $ref.fasta $gene > ./fasta.files/$ref.$gene.fa
samtools view -b /PATH/100x.s1.bam $gene > ./bam.files/s1.$gene.bam # input from mapping required either use downsampled or realigned bam-file
samtools index ./bam.files/s1.$gene.bam # creates a bam-file for each gene
samtools view -h -F 4 -b ./bam.files/s1.$gene.bam > ./bam.files/s1.$gene.onlymapped.bam # keeps only the gene region
samtools index ./bam.files/s1.$gene.onlymapped.bam
rm ./bam.files/s1.$gene.bam
done

for gene in $gene_pos
do
calc=$(echo $gene | sed 's/.*://g')
a=$(echo $((-1*($calc)))) # calculates gene length
timeout 5m Rscript ./ViQuaS.R ./fasta.files/$ref.$gene.fa ./bam.files/s1.$gene.onlymapped.bam 5 0.7 1 $a
mv ViQuaS-Richness.txt Richness.s1.$gene.txt
mv ViQuaS-Spectrum.fa Spectrum.s1.$gene.fa
rm -r SIM
rm -r Samples
rm -r viquas
done


echo =================================================================
echo "Script ended $(date)"
DATE=`date +%d%m%Y`
echo "Your logfile is geneHaplotyping."$DATE".log" 
echo =================================================================



