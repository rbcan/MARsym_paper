# 01. Mapping: mapping.bbmapv36x.sh 

## [01] Purpose of script:
Shows detailed commands that were used for read mapping in the MARsym paper.
Reads of quality q20 were mapped to a consensus reference sequence with minimum nucleotide identity of 0.95 using BBMap.

## [01] Programs that need to be installed to execute script:
- BBmap v36.x https://sourceforge.net/projects/bbmap/ (for later versions command details might need to be adatpted)
- samtools https://github.com/samtools/samtools

## [01] Required input files:
- Illumina reads in gzipped fastq format for each sample s1, s2, s3, ... sn (e.g. s1.fq.gz)
- Reference sequence in fasta format: ref.fasta
info: the script uses a consensus reference for all samples. If you have one reference per sample, you need to replace all $ref in the code with $sample

## [01] Example of output files for sample s1
- mapping.$DATE.log - *logfile for mapping (all samples s1, s2, s3, ... sn)*
- coverage_depth_bam.txt - *coverage depths of the mapping file of each sample (s1, s2, s3, ... sn)*
- coverage_depth_rmdup.txt - *coverage depths of the mapping file of each sample (s1, s2, s3, ... sn) after PCR duplicate removal*
- s1.q20.bam - *sorted and indexed mapping file*
- s1.rmdup.bam - *sorted and indexed mapping file after PCR duplicate removal*


# 02. SNPcalling: SNPcalling.gatkv3.3.0.sh

## [02] Purpose of script:
Shows detailed commands that were used for SNP calling in the MARsym paper.
Script calls SNPs on mapped reads to a reference with GATK and a ploidy setting of 10. Input files can be created with MARsym_mapping.

## [02] Programs that need to be installed to execute script:
- GenomeAnalysisToolKit (GATK) v3.3.0 https://software.broadinstitute.org/gatk/download/ (for later versions command details might need to be adatpted)
- picard-tools v1.119 https://github.com/broadinstitute/picard (for later versions command details might need to be adatpted)
- samtools https://github.com/samtools/samtools
- ucsc tools (executable faCount) https://github.com/adamlabadorf/ucsc_tools

## [02] Required input files:
- Reference sequence in fasta format: ref.fasta
info: the script uses a consensus reference for all samples. If you have one reference per sample, you need to replace all $ref in the code with $sample
- indexed bamfile of reads mapping to <ref>.fasta where PCR duplicates were removed in sorted and indexed bam format: s1.id95.rmdup.bam s2.id95.rmdup.bam s3.id95.rmdup.bam ... sn.id95.rmdup.bam, where s1, s2, s3, sn are the individual $samples. Input bam-files can be created with MARsym_mapping

## [02] Example of output files for sample s1 with target read coverage of 100x
- ref.dict - *reference dictionary (for all samples s1, s2, s3, ... sn)*
- SNPcounts.txt - *final SNP counts: absolute and per kbp (for all samples s1, s2, s3, ... sn)*
- SNPcalling.$DATE.log - *logfile (for all samples s1, s2, s3, ... sn)*
- s1.real.bam - *bamfiles with readgroups and realigned reads around INDELs*
- 100x.s1.bam - *downsampled bam file (with readgroups and realigned reads around INDELs) to target read coverage 100x*
- 100x.s1.rawVar_q30_ploidy10.vcf - *raw variants called with ploidy 10*
- 100x.s1.rawSNPs_ploidy10.vcf - *raw SNPs called with ploidy 10*
- 100x.s1.rawINDELs_ploidy10.vcf - *raw INDELs called with ploidy 10*
- 100x.s1.filtSNPs_ploidy10.vcf - *filtered SNPs: SNPs that don't pass the filter are flagged with corresponding filter(s)*
- 100x.s1.filtINDELs_ploidy20.vcf - *filtered INDELs: SNPs that don't pass the filter are flagged with corresponding filter(s)*
- 100x.s1.filtSNPs_PASS_ploidy10.vcf - *file with only those SNPs that passed the all filters* 
- 100x.s1.filtINDELs_PASS_ploidy10.vcf - *file with only those INDELs that passed the all filters* 

# 03. Strain number estimation: geneHaplotyping.viquas1.3.sh

## [03] Purpose of script:
Shows detailed commands that were used for estimating the number of gene versions for the provided set of genes using the tool ViQuaS in the MARsym paper.

## [03] Programs that need to be installed to execute script:
- ViQuaS https://academic.oup.com/bioinformatics/article/31/6/886/215466
- samtools https://github.com/samtools/samtools

## [03] Required input files:
- Reference sequence in fasta format: ref.fasta
info: the script uses a consensus reference for all samples. If you have one reference per sample, you need to replace all $ref in the code with $sample
- s1.real.bam - bamfiles with readgroups and realigned reads around INDELs
  OR (depends on sample/question)
- 100x.s1.bam - downsampled bam file (with readgroups and realigned reads around INDELs) with target read coverage 100x

## [03] Example of output files for sample s1
- Spectrum-file contains reconstructed fasta sequences and abundance of sequence
- Richness-file contains f_min value. We discarded all reconstructed sequences below this frequency

# 04. Identification of low-coverage genes
Shows detailed commands that were used for identification of genes with coverage below the range of coverage from gammaproteobacterial marker genes. These genes were classified as strain-specific in the MARsym paper. 

## [04] Programs that need to be installed to execute script:
- samtools https://github.com/samtools/samtools
- ucsc tools (executable faCount) https://github.com/adamlabadorf/ucsc_tools
- bedtools https://github.com/arq5x/bedtools2
- R https://www.r-project.org/
- PhylaAmphora https://github.com/martinwu/Phyla_AMPHORA

## [04] Required input files:
- Reference sequence in fasta format: ref.fasta
- Annotations of reference in gff3 format: ref.gff
- Predicted amino acid sequences of proteins: ref.faa
- s1.real.bam - bamfiles with readgroups and realigned reads around INDELs
  OR (depends on sample/question)
- 100x.s1.bam - downsampled bam file (with readgroups and realigned reads around INDELs) with target read coverage 100x
