# MARsym_snpcalling

# Purpose of script:
Shows detailed commands that were used for SNP calling in the MARsym paper.
Script calls SNPs on mapped reads to a reference with a ploidy setting of 10. Input files can be created with MARsym_mapping.

# Programs that need to be installed to execute script:
- GenomeAnalysisToolKit (GATK) v3.3.0 https://software.broadinstitute.org/gatk/download/ (for later versions command details might need to be adatpted)
- picard-tools v1.119 https://github.com/broadinstitute/picard (for later versions command details might need to be adatpted)
- samtools https://github.com/samtools/samtools
- ucsc tools (executable faCount) https://github.com/adamlabadorf/ucsc_tools

# Required input files:
- Reference sequence in fasta format: ref.fasta
info: the script uses a consensus reference for all samples. If you have one reference per sample, you need to replave all $ref in the code with $sample
- indexed bamfile of reads mapping to <ref>.fasta where PCR duplicates were removed in sorted and indexed bam format: s1.id95.rmdup.bam s2.id95.rmdup.bam s3.id95.rmdup.bam ... sn.id95.rmdup.bam, where s1, s2, s3, sn are the individual $samples. Input bam-files can be created with MARsym_mapping

# Output files of sample s1
- ref.dict - reference dictionary (for all samples s1, s2, s3, ... sn) 
- SNPcounts.txt - final SNP counts: absolute and per kbp (for all samples s1, s2, s3, ... sn)
- SNPcalling.$DATE.log - logfile (for all samples s1, s2, s3, ... sn)
- s1.real.bam - bamfiles with readgroups and realigned reads around INDELs
- $downx.s1.bam - downsampled bam file (with readgroups and realigned reads around INDELs) to target read coverage $downx
- $downx.s1.rawVar_q30_ploidy10.vcf - raw variants called with ploidy 10
- $downx.s1.rawSNPs_ploidy10.vcf - raw SNPs called with ploidy 10
- $downx.s1.rawINDELs_ploidy10.vcf - raw INDELs called with ploidy 10
- $downx.s1.filtSNPs_ploidy10.vcf - filtered SNPs: SNPs that don't pass the filter are flagged with corresponding filter(s)
- $downx.s1.filtINDELs_ploidy20.vcf - filtered INDELs: SNPs that don't pass the filter are flagged with corresponding filter(s)
- $downx.s1.filtSNPs_PASS_ploidy10.vcf - file with only those SNPs that passed the all filters 
- $downx.s1.filtINDELs_PASS_ploidy10.vcf - file with only those INDELs that passed the all filters 


