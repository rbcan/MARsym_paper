#!/bin/bash

#log into 15 threads with qrsh or try qsub

DATE=`date +%d%m%Y`
logfile=mapping."$DATE".log
exec > $logfile 2>&1

echo =================================================================
echo "Script started $(date)"
DATE=`date +%d%m%Y`
echo "Your logfile is mapping."$DATE".log" 
echo =================================================================

ref="<ref>" # Replace <ref> with reference name for which you have <ref>.fasta
### info: this script uses a consensus reference for all samples. If you have one reference per sample, you need to replave all $ref in the code below with $sample 
SAMPLES="<s1> <s2> <s3> ... <sn>" # Replace <s1> <s2> <s3> ... <sn> with your sample IDs separated by a space each. You have to provide sn.f.fq.gz and sn.r.fq.gz for each of your samples n.
out="mapping" # Specify name variable for outputfolder

echo " "
echo "your reference is "$ref
echo " "
for sample in $SAMPLES
do
echo " "
echo "your read inputs are: "$sample".f.fq.gz, "$sample".r.fq.gz"
echo " "
done

############ Create coverage depth files
echo "coverage depth for *.q20.bam:" > coverage_depth_bam.txt
echo "coverage depth for *.q20.bam:" > coverage_depth_rmdup.txt

############ If needed, set read quality to q20. Not always advisable - depends on sample type.
echo "$(date)"
echo "set read quality to 20"

for sample in $SAMPLES
do
/PATH/bbduk.sh threads=15 in1=$sample.f.fq.gz in2=$sample.r.fq.gz out=$sample.q20.fq.gz qtrim=rl trimq=20
done

############ mapping with bbmap - adjust read name if not using q20 reads
echo ====================================================================================================
echo "$(date)"
echo "start mapping with bbmap"
echo " " 

for sample in $SAMPLES
  do
    echo "*****************************************************************"
    echo "mapping and sorting $sample"
    echo "*****************************************************************"
    /PATH/bbmap.sh threads=15 ref=$ref.fa in=$sample.q20.fq.gz out=$sample.q20.bam
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
echo "Your logfile is mapping."$DATE".log" 
echo =================================================================

