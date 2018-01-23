#!/bin/bash

#log into 15 threads with qrsh or try qsub

echo "Please enter your input REF reference name (of which you have REF.fasta), include PATH if not in your working directory: "
read REF
echo "your sample read files should be named like this: SAMPLE_NAME.f.fq.gz (forward) and SAMPLE_NAME.r.fq.gz (reverse)"
echo "Please enter your input SAMPLE name (s) (e.g. "BazF BazG BazH"): "
read SAMPLES

echo =================================================================
echo "Script started $(date)"
DATE=`date +%d%m%Y`
echo "Your logfile is "$REF".mapping.$DATE.log" 
echo =================================================================


logfile=$REF.mapping.$DATE.log
exec > $logfile 2>&1

echo =================================================================
echo "Script started $(date)"
DATE=`date +%d%m%Y`
echo "Your logfile is "$REF".mapping.$DATE.log" 
echo =================================================================

echo " "
echo "your reference is "$REF
echo "your read inputs are: "$sample".f.fq.gz, "$sample".r.fq.gz"
echo " "

### Adapt your the files according to your files:
### replace READs_f and READS_r with your clean forward and reverse reads (I have reads with quality 20 usually)

############ Create coverage depth files
echo "coverage depth for *.q20.bam:" > coverage_depth_bam.txt
echo "coverage depth for *.q20.bam:" > coverage_depth_rmdup.txt

############ Modify read headers
#echo "modify read headers"
#sed 's/ 1:N:0:[0-9]*/\/1/g' input.reads.R1.fastq > output.reads.R1.fastq
#sed 's/ 2:N:0:[0-9]*/\/2/g' input.reads.R2.fastq > output.reads.R2.fastq

############ Quality 20
echo "$(date)"
echo "set read quality 20"

for sample in $SAMPLES
do
bbduk.sh threads=10 in1=$sample.f.fq.gz in2=$sample.r.fq.gz out=$sample.q20.fq.gz qtrim=rl trimq=20
#bbduk.sh -Xmx8g ref=/home/rebecca/work/tools/bbmap/resources/phix174_ill.ref.fa.gz,/home/rebecca/work/tools/bbmap/resources/truseq.fa.gz in1=$sample.f.fq.gz in2=$sample.r.fq.gz out=$sample.q20.fq.gz qtrim=rl trimq=20
#reformat.sh -Xmx8g in=$sample.q20.fq.gz out1=$sample.fq.gz ftl=10
done

  
############ mapping with bbmap
echo ====================================================================================================
echo "$(date)"
echo "mapping with bbmap in progress"
echo " " 

for sample in $SAMPLES
  do
    echo "*****************************************************************"
    echo "mapping and sorting $sample"
    echo "*****************************************************************"
    bbmap.sh threads=15 ref=$REF.fa in=$sample.q20.fq.gz out=$sample.q20.bam
    samtools sort $sample.q20.bam $sample.q20S
    mv $sample.q20S.bam $sample.q20.bam
    ########## index *.bam file
    samtools index $sample.q20.bam
    ########## remove PCR duplicates
    echo "*****************************************************************"
    echo "remove PCR duplicates $sample"
    echo "*****************************************************************"
    samtools rmdup -S $sample.q20.bam $sample.rmdup.bam
    ########## check the coverage depth and the change after the duplicate removal
    echo "*****************************************************************"
    echo "check coverage depth $sample"
    echo "*****************************************************************"
    samtools depth $sample.q20.bam | awk -v var="$sample" '{sum+=$NF;cnt++}END{print $sample" "sum/cnt}' >> coverage_depth_bam.txt
    samtools depth $sample.rmdup.bam | awk -v var="$sample" '{sum+=$NF;cnt++}END{print $sample".rmdup "sum/cnt}' >> coverage_depth_rmdup.txt
  done

echo " "
echo ====================================================================================================

echo mapping finished!
echo =================================================================
echo "Script finished $(date)"
DATE=`date +%d%m%Y`
echo "Your logfile is "$REF".mapping.$DATE.log" 
echo =================================================================

