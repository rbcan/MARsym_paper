#!/bin/bash

source ~/.bashrc

DATE=`date +%d%m%Y`
logfile=SNPcalling.$DATE.log
exec > $logfile 2>&1

echo =================================================================
echo "Script started $(date)"
DATE=`date +%d%m%Y`
echo "Your logfile is SNPcalling."$DATE".log" 
echo =================================================================


ref="<ref>" # Replace <ref> with reference name for which you have <ref>.fasta
### info: this script uses a consensus reference for all samples. If you have one reference per sample, you need to replave all $ref in the code below with $sample 
SAMPLES="<s1> <s2> <s3> ... <sn>" # Replace <s1> <s2> <s3> ... <sn> with your sample IDs separated by a space each. If you ran MARsym_mapping beforehand your outputs should be <s1>.id95.rmdup.bam <s2>.id95.rmdup.bam <s3>.id95.rmdup.bam ... <sn>.id95.rmdup.bam. Rename in code if your mapped bam-files are named differently
out="SNPs" # Specify name variable for outputfolder
down="100" # Read coverage for which you want to downsample. Ideally this is the maximum coverage of your sample with the lowest coverage.


# create lists for the SNP count output
echo "SNP counts" > SNPcounts.txt
echo -e "\tsample\tSNPs\tSNPsPerkbp" >> SNPcounts.txt

# start SNP calling prep and procedure
for sample in $SAMPLES
        do
        ########## create reference dictionary, necessary for the use of GATK
        ########## tools needed: picard-tools
        ##### important: the <ref>.dict file has to be in the same folder as <ref>.fasta)

        echo "create reference dictionary"
	samtools faidx $ref.fasta
        java -jar /PATH/picard-tools-1.119/CreateSequenceDictionary.jar R=$ref.fasta O=$ref.dict

        echo ====================================================================================================

        ########## add headers to bamfile, required by GATK
        ########## tools needed: picard-tools
        ##### use indexed bamfile where PCR duplicates were removed

        echo $(date)" add readgroup headers to bamfile "$sample

        java -jar /PATH/picard-tools-1.119/AddOrReplaceReadGroups.jar I=$sample.id95.rmdup.bam O=$sample.RG.bam ID=$sample SM=$sample LB=$sample PL=illumina PU=barcode VALIDATION_STRINGENCY=LENIENT

        echo $(date)" finished"

        echo ====================================================================================================

        echo $(date)" index bamfile "$sample

        samtools index $sample.RG.bam
        samtools faidx $sample.fasta

        echo ====================================================================================================

        ## Step1:
        ########## Realign reads around INDELS (equivalent to the program coval)
        echo $(date)" Realign around indels"
        ########## create targets

        echo $(date)" create targets "$sample

        java -jar /PATH/GenomeAnalysisTK.jar \ # was used with GATK v3.3.0, McKenna et al. (2010)
        -T RealignerTargetCreator \
        -R ./$ref.fasta \
        -I ./$sample.RG.bam \
        --minReadsAtLocus 4 \
        -o ./$sample.intervals

        echo ====================================================================================================

        ########## realgining around indel intervals
        # LOD: threshold above lewhich cleaner should clean

        echo $(date)" realgining around indel intervals "$sample

        java -jar /PATH/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R ./$ref.fasta \
        -I ./$sample.RG.bam \
        -targetIntervals ./$sample.intervals \
        -LOD 3.0 \
        -o ./$sample.real.bam 

        echo ====================================================================================================

        ## Step 1.2:
        ########## Cutting the average coverage to the same value to ensure comparability between samples
        depth=$(samtools depth $sample.real.bam | awk '{sum+=$3;cnt++}END{print sum/cnt}')
        factor=$(python -c "print $down/$depth")
        samtools view -s $factor -u $sample.real.bam -o $down\x.$sample.bam
        samtools index $down\x.$sample.bam

        echo ====================================================================================================

        ## Step2:
        ########## SNP and INDEL calling
        ##### if needed, change settings for ploidy and alternate alleles
        ##### if the coverage step 1.2 was NOT done, please change $down\x.$sample into $sample.real
        # -stand_call_conf (min phred-scaled confidence at which variants should be called: 30)
        # -stand_emit_conf (min phred-sclaed confidence at which variants should be emitted --> if below stand_call: LowQual: 10)
        # -max_alternate_alleles (max number of alternate alleles to genotype)
        # -ploidy (ploidy per sample / number of strains: 10)
        # -gt_mode DISCOVERY (to discover new alternate alleles)
        echo $(date)" call variants with HaplotypeCaller "$sample
        echo "variants with confidence >10 and <30 --> LowQual"
        echo "variants with confidence <10 --> not called"
        echo "ploidy is 10"

        java -jar -Xmx30G /PATH/GenomeAnalysisTK.jar \ # adapt -Xmx30G as needed
        -T HaplotypeCaller \
        -R ./$ref.fasta \
        -I ../$down\x.$sample.id95.bam \
        -gt_mode DISCOVERY \
        -stand_call_conf 30 \
	-stand_emit_conf 10 \
        -ploidy 10 \
        --max_alternate_alleles 10 \
        -o $down\x.$sample.rawVar_q30_ploidy10.vcf

        echo ====================================================================================================

        ## Step3: extract only
        ########## extract SNPs/INDELs

        ########## SNPs
        echo "extract only SNPs"

        echo "ploidy 10"
        java -jar /PATH/GenomeAnalysisTK.jar \
        -T SelectVariants \
        -R ./$ref.fasta \
        -V $down\x.$sample.rawVar_q30_ploidy10.vcf \
        -selectType SNP \
        -o $down\x.$sample.rawSNPs_ploidy10.vcf

        echo ====================================================================================================

        ########## INDELs
        echo "extract only INDELs"

        echo "ploidy 10"
        java -jar /PATH/GenomeAnalysisTK.jar \
        -T SelectVariants \
        -R ./$ref.fasta \
        -V $down\x.$sample.rawVar_q30_ploidy10.vcf \
        -selectType INDEL \
        -o $down\x.$sample.rawINDELs_ploidy10.vcf

        echo ====================================================================================================

        ## Step4:
        ########## Filter variants with hardFilter approach

        ########## SNPs
        # note to not use haplotype filter (only for diploid organisms)
        # MQRankSum to 20 (suggested was 12.5 but for my SNPs it filtered out too many. After manual checks 20 seems about right for 16S, 23S and recA
        # note, the SNPs are still in the file but marked with the filter that excluded them. To esctract only SNPs that passed all filters do step 5
        echo "filtering SNPs"
        echo "filters: QD < 2.0, FS > 60.0, MQ < 40.0, MQRankSum < -20, ReadPosRankSum < -8"

        echo "ploidy 10"
        java -jar /PATH/GenomeAnalysisTK.jar \
        -T VariantFiltration \
        -R ./$ref.fasta \
        -V $down\x.$sample.rawSNPs_ploidy10.vcf \
        --filterExpression "QD < 2.0" \
        --filterName "QualByDepth" \
        --filterExpression "FS > 60.0" \
        --filterName "Fisher" \
        --filterExpression "MQ < 40.0" \
        --filterName "RMS_mapQual" \
        --filterExpression "MQRankSum < -20" \
        --filterName "MQRS" \
        --filterExpression "ReadPosRankSum < -8" \
        --filterName "RPRS" \
        -o $down\x.$sample.filtSNPs_ploidy10.vcf

       echo ====================================================================================================

        ########## INDELs
        echo "filtering INDELs"
        echo "filters QD < 2.0, FS > 200.0, ReadPosRankSum < -20.0"#
        echo "ploidy 10"
        java -jar /PATH/GenomeAnalysisTK.jar \
        -T VariantFiltration \
        -R ./$ref.fasta \
        -V $down\x.$sample.rawINDELs_ploidy10.vcf \
        --filterExpression "QD < 2.0" \
        --filterName "QualByDepth" \
        --filterExpression "FS > 200.0" \
        --filterName "Fisher" \
        --filterExpression "ReadPosRankSum < -20.0" \
        --filterName "RPRS" \
        -o $down\x.$sample.filtINDELs_ploidy20.vcf

        echo ====================================================================================================

        ## Step 5:
        ########## save only variants that "PASS" filters into separate file
        echo "save Vars that passed all filters in separate files"

        cat $down\x.$sample.filtSNPs_ploidy10.vcf | grep 'PASS\|^#' > ./$down\x.$sample.filtSNPs_PASS_ploidy10.vcf
        cat $down\x.$sample.filtINDELs_ploidy10.vcf | grep 'PASS\|^#' > ./$down\x.$sample.filtINDELs_PASS_ploidy10.vcf
        echo ====================================================================================================

done

for sample in $SAMPLES
        do
        ########## count variants
        echo "variant numbers saved to SNPcounts.txt"

	genomelength=$(/PATH/ucsc_tools-master/executables/faCount $ref.fasta | tail -1 | awk '{print $2-$(NF-1)}')
       	rawploidy10=$(grep -v "^#" ./$down\x.$sample.rawSNPs_ploidy10.vcf -c | awk '{print $0"\t"$1/"'$genomelength'"*1000}')
        echo -e "rawSNPs_pl10\t"$sample"\t" $rawploidy10 >> SNPcounts.txt
        rawINDploidy10=$(grep -v "^#" ./$down\x.$sample.rawINDELs_ploidy20.vcf -c | awk '{print $0"\t"$1/"'$genomelength'"*1000}')
        echo -e "rawINDELs_pl10\t"$sample"\t"$rawINDploidy10 >> SNPcounts.txt

        filtploidy10=$(grep -v "^#" ./$down\x.$sample.filtSNPs_PASS_ploidy10.vcf -c | awk '{print $0"\t"$1/"'$genomelength'"*1000}')
	echo -e "filtSNPs_pl10\t"$sample"\t"$filtploidy10 >> SNPcounts.txt
        filtINDploidy10=$(grep -v "^#" ./$down\x.$sample.filtINDELs_PASS_ploidy10.vcf -c | awk '{print $0"\t"$1/"'$genomelength'"*1000}')
	echo -e "filtINDELs_pl10\t"$sample"\t"$filtINDploidy10 >> SNPcounts.txt

done

head -50 SNPcounts.txt

echo =================================================================
echo "Script finished $(date)"
DATE=`date +%d%m%Y`
echo "Your logfile is SNPcalling."$DATE".log" 
echo =================================================================

