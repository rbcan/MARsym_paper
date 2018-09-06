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

ref="<ref>" # Replace <ref> with reference name for which you have ref.fasta
sample="s1" # this example runs with sample s1, replace with your corresponding sample

cd /PATH/ViQuaS1.3/ # this folder contains the ViQuaS program parts - important to run script in that folder

mkdir ./fasta.files/
mkdir ./bam.files/

# preparation file $sample.viquasIn for all samples s1, s2, s3, ... sn: Contig <tab> geneStart <tab> geneEnd; e.g. Node_15  500 1500


gene_pos=$(awk -F "\t" '{print $2":"$3"-"$4}' $sample.viquasIn | tr '\n' '\t')

for gene in $gene_pos
do
samtools faidx $ref.fasta $gene > ./fasta.files/$ref.$gene.fa
samtools view -b /PATH/100x.$sample.bam $gene > ./bam.files/$sample.$gene.bam # input from mapping required either use downsampled or realigned bam-file
samtools index ./bam.files/$sample.$gene.bam # creates a bam-file for each gene
samtools view -h -F 4 -b ./bam.files/$sample.$gene.bam > ./bam.files/$sample.$gene.onlymapped.bam # keeps only the gene region
samtools index ./bam.files/$sample.$gene.onlymapped.bam
rm ./bam.files/$sample.$gene.bam
done

for gene in $gene_pos
do
calc=$(echo $gene | sed 's/.*://g')
a=$(echo $((-1*($calc)))) # calculates gene length
timeout 5m Rscript ./ViQuaS.R ./fasta.files/$ref.$gene.fa ./bam.files/$sample.$gene.onlymapped.bam 5 0.7 1 $a
mv ViQuaS-Richness.txt Richness.$sample.$gene.txt
mv ViQuaS-Spectrum.fa Spectrum.$sample.$gene.fa
rm -r SIM
rm -r Samples
rm -r viquas
done


echo =================================================================
echo "Script ended $(date)"
DATE=`date +%d%m%Y`
echo "Your logfile is geneHaplotyping."$DATE".log" 
echo =================================================================



