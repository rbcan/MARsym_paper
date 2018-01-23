#!/bin/bash



mkdir ./fasta.files/
mkdir ./bam.files/

# sample.viquasIn: Contig <tab> geneStart <tab> geneEnd; e.g. Node_15  500 1500
gene_pos=$(awk -F "\t" '{print $2":"$3"-"$4}' s1.viquasIn | tr '\n' '\t')

for gene in $gene_pos
do
samtools faidx $ref.fasta $gene > ./fasta.files/$ref.$gene.fa
samtools view -b ./100x.s1.map.real.bam $gene > ./bam.files/s1.$gene.bam
samtools index ./bam.files/s1.$gene.bam
samtools view -h -F 4 -b ./bam.files/s1.$gene.bam > ./bam.files/s1.$gene.onlymapped.bam
samtools index ./bam.files/s1.$gene.onlymapped.bam
rm ./bam.files/s1.$gene.bam
done

for gene in $gene_pos
do
calc=$(echo $gene | sed 's/.*://g')
a=$(echo $((-1*($calc))))
timeout 5m Rscript ./ViQuaS.R ./fasta.files/$ref.$gene.fa ./bam.files/s1.$gene.onlymapped.bam 5 0.7 1 $a
mv ViQuaS-Richness.txt Richness.s1.$genetxt
mv ViQuaS-Spectrum.fa Spectrum.s1.$gene.fa
rm -r SIM
rm -r Samples
rm -r viquas
done



