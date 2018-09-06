#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 5
#$ -l h_vmem=50G

source ~/.bashrc

echo =================================================================
echo "Script started $(date)"
DATE=`date +%d%m%Y`
echo =================================================================

logfile=$sample.lowcov.$DATE.log
exec > $logfile 2>&1

mkdir -p lowcov
cd lowcov
sample_names="<s1> <s2> <s3> ... <sn>" # Replace <s1> <s2> <s3> ... <sn> with your sample IDs separated by a space each. You need to provide sn.bam
ref="<ref>" # Replace <ref> with reference name for which you have ref.fasta, ref.gff (gff3 file containing annotions), ref.faa (amino acid seuqences for annotated proteins)
### info: this script uses a consensus reference for the samples. If you have one reference per sample, you need to replave all $ref in the code below with $sample 
genome_files="/path/to/genome_files" # folder to your genome files ref.fasta, ref.gff, ref.faa
bamfolder="/path/to/genome_files" # folder to your sorted read alignment files sn.bam and index files sn.bai
tools="/path/to/tools" # folder to tools PhylaAmphora, ucsc tools

# step 1
# calculate coverage per base for each sample using bedtools
# required input files: bam file with reads aligned to reference, gff file with annotation
 echo "1. coverage per base calculations for $sample_names"
 for sample in $sample_names
 do
 echo "=> processing sample $sample"
 bedtools coverage -abam $bamfolder/$sample.bam -b $genome_files/$ref.gff -d > CovPerBase_$sample.txt
 awk -F ";Name=" '{print $1"\t"$NF}' CovPerBase_$sample.txt | awk -F "\t" '{print $9"\t"$NF}' > xx && mv xx CovPerBase_$sample.txt
 echo -e "gene\tcoverage" | cat - CovPerBase_$sample.txt > xx && mv xx CovPerBase_$sample.txt
 done

# step 2
# calculate coverage per gene for each sample using R
# required input files: CovPerBase file from step 1
echo "2. coverage per gene calculations for $sample_names"
for sample in $sample_names
do
echo "=> processing sample $sample"

echo "data <-read.table("INFILE.txt", sep="\t", header=TRUE)
genes = unique(data$gene, incomparables=FALSE)
outfile <- rbind(c("genes", "coverage"))
for(i in 1:(length(genes))){
	positions <- which(data$gene%in%genes[i])
	mean <- mean(data$coverage[positions])
	outfile <- rbind(outfile, rbind(c((as.character(genes[i])), mean)))
}
write.table(outfile, sep="\t", file="OUTFILE.txt")
proc.time()" > ./AverageCovPerGene_$sample\temp.R

sed -i "s/INFILE/CovPerBase_$sample/g" AverageCovPerGene_$sample\temp.R
sed -i "s/OUTFILE/CovPerGene_$sample/g" AverageCovPerGene_$sample\temp.R
R CMD BATCH AverageCovPerGene_$sample\temp.R
rm AverageCovPerGene_$sample\temp.R
done

# step 3
# extract marker genes with phyla_AMPHORA 
# input files: $ref.faa $ref.fa
echo "3. extracting marker genes using Phyla_AMPHORA for $sample_names"
for sample in $sample_names
do
mkdir phylamphgenes_$sample
cd phylamphgenes_$sample
perl $tools/Phyla_AMPHORA-master/Scripts/MarkerScanner.pl -Phylum 3 -DNA $genome_files/$ref.fa # I used Phylum 3 (gammaproteobacteria) - needs to be adjusted according to the target group
perl $tools/Phyla_AMPHORA-master/Scripts/MarkerAlignTrim.pl -WithReference -OutputFormat phylip
perl $tools/Phyla_AMPHORA-master/Scripts/Phylotyping.pl -CPUs 3 > $sample.phylotype.result
cat *.pep > $sample.all_markers.fa
grep ">" $sample.all_markers.fa | sed 's/>//g' > all_marker_ids.txt
makeblastdb -dbtype prot -in $sample.all_markers.fa
blastp -query $genome_files/$ref.faa -db $sample.all_markers.fa -evalue 0.0001 -out blast_$sample\_all_markers -outfmt 6
awk '{if ($3~"100.00") print $0}' blast_$sample\_all_markers > $sample.markers.txt
sed 's/ID=//g' $genome_files/$ref.gff | awk -F ";Name=" '{print $1"\t"$2}' | awk -F "\t" '{print $9"\t"$NF"\t"$4"\t"$5"\t"$1}' > $sample\_genelist
sed -i -e "1d" $sample\_genelist
sort -k1,1 $sample.markers.txt > xx && mv xx $sample.markers.txt
sort -k1,1 $sample\_genelist > xx && mv xx  $sample\_genelist
awk '{print $1}' $sample.markers.txt | join - $sample\_genelist -t '	' | uniq > PhyAMPH_genes100ID_$sample # TAB separator!!
cd ..
done

# step 4
# produce files for plotting and extracting genes/singlecopy genes based on their coverage
# required input files: CovPerBase file from step 2, phylamphgenes/$sample.markers.txt from step 3
echo "4. file formatting for $sample_names"
for sample in $sample_names
do
echo "=> processing sample $sample"
sed -i "1d" CovPerGene_$sample.txt
sed 's/"//g' CovPerGene_$sample.txt > xx && mv xx CovPerGene_$sample.txt
awk '{print $2"\t"$3}' CovPerGene_$sample.txt > xx && mv xx CovPerGene_$sample.txt
sort -k1,1 CovPerGene_$sample.txt > CovPerGene_$sample\_sorted.txt
awk -F ";Name=" '{print $1"\t"$2}' $genome_files/$ref.gff | awk -F "\t" '{print $(NF-1)"\t"$1"\t"$4"\t"$5"\t"$NF}' - | sort -k1,1 - | join - CovPerGene_$sample\_sorted.txt -t '	' > CovPerGene_$sample\_annot.txt # '	' has to be a tab
$tools/ucsc_tools-master/executables/faCount $genome_files/$ref.fa | awk '{if ($1!~"#" && $1!~"total") print $1"\t"$2}' | sort -k1,1  > contig.length.$sample
sort -k2,2 CovPerGene_$sample\_annot.txt | join -1 2 -2 1 - contig.length.$sample -t '	' | awk -F "\t" '{print $2"\t"$1"\t"$3"\t"$4"\t"$7"\t"$5"\t"$6}' | sort -k1,1 > xx && mv xx CovPerGene_$sample\_annot.txt # '	' has to be a tab
awk '{print $NF"\t"$0}' CovPerGene_$sample\_annot.txt | sort -n | awk -F "\t" '{print NR"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' > CovPerGene_$sample.plot
echo -e "order\tID\tcontig\tstart\tend\tcontlength\tgene\tcoverage" | cat - CovPerGene_$sample.plot > xx && mv xx CovPerGene_$sample.plot # echo -e
sort -k2,2 CovPerGene_$sample.plot > CovPerGene_$sample.plot.sort
awk '{if ($4>=100) print $1}' phylamphgenes_$sample/$sample.markers.txt | sed 's/fig|/ID=fig|/g' | sort | uniq | sort -k1,1 | join -1 1 -2 2 - CovPerGene_$sample.plot.sort -t '	' > CovPerGene_$sample.singlecopy.plot # '	' has to be a tab
echo -e "ID\torder\tcontig\tstart\tend\tcontlength\tannotation\tcoverage" | cat - CovPerGene_$sample.singlecopy.plot  > xx && mv xx CovPerGene_$sample.singlecopy.plot # echo -e
awk '{if ($0!~"hypothetical protein") print $0}' CovPerGene_$sample.singlecopy.plot > xx && mv xx CovPerGene_$sample.singlecopy.plot
awk '{if ($4>=100 || $5<=($6-100)) print $0}' CovPerGene_$sample.singlecopy.plot > xx && mv xx CovPerGene_$sample.singlecopy.plot
#### duplicate annotations removal:
awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8"\t"$7}' CovPerGene_$sample.singlecopy.plot | sed 's/ /xxx/g' | sort -k8 | uniq -u -f7 | sed 's/xxx/ /g' | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8"\t"$7}' | sort -k1 -r > xx && mv xx CovPerGene_$sample.singlecopy.plot 
done

# step 5
# define low coverage genes using R script
# required input files: CovPerGene_$sample.plot CovPerGene_$sample.singlecopy.plot
echo "5. define low coverage genes using R script for $sample_names"

echo "creating RplottingLowCov.R"

echo "library(ggplot2, lib.loc="'"~/R/"'")" > RplottingLowCov.R
echo "library(data.table, lib.loc="'"~/R/"'")" >> RplottingLowCov.R
echo "library(labeling, lib.loc="'"~/R/"'")" >> RplottingLowCov.R
echo "library(digest, lib.loc="'"~/R/"'")" >> RplottingLowCov.R


for sample in $sample_names
do
echo "" >> RplottingLowCov.R
echo "AllGenes"$sample" <- fread('./CovPerGene_"$sample".plot')" >> RplottingLowCov.R
echo "SingleCopy"$sample" <- fread('./CovPerGene_"$sample".singlecopy.plot')" >> RplottingLowCov.R
echo "setnames(AllGenes"$sample", 'gene', 'ID')" >> RplottingLowCov.R
echo "AllGenes"$sample"[ , Type := 'AllGenes']" >> RplottingLowCov.R
echo "SingleCopy"$sample"[ , Type := 'SingleCopy']" >> RplottingLowCov.R
echo "str(AllGenes"$sample")" >> RplottingLowCov.R
echo "str(SingleCopy"$sample")" >> RplottingLowCov.R
echo "bind."$sample" <- rbind(AllGenes"$sample", SingleCopy"$sample", fill = TRUE)" >> RplottingLowCov.R
done

for sample in $sample_names
do
echo "
library(gridExtra, lib.loc="'"~/R/"'")
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "'"guide-box"'")
  legend <- tmp$grobs[[leg]]
  return(legend)
}" >> RplottingLowCov.R
done

for sample in $sample_names
do
echo "
box"$sample" <- ggplot(
  data = "bind.$sample",
  mapping = aes(
    x = Type,
    y = coverage,
    fill = Type
    )
) + geom_boxplot(
  outlier.size = 1,
  #outlier.colour = 'black'
) + theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="'"none"'"
) + scale_y_log10()

scat"$sample" <- ggplot(
  data = "bind.$sample",
  mapping = aes(
    x = order,
    y = coverage,
    colour = Type
    )
  ) + geom_point(
    size = 1
    ) + theme(
      legend.position="'"right"'",
      legend.title=element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
      ) + scale_y_log10()
">> RplottingLowCov.R
done

echo "" >> RplottingLowCov.R

insamp=$(for sample in $sample_names; do echo "scat"$sample", box"$sample","; done | tr '\n' ' ')

for sample in $sample_names
do
echo "scat"$sample" <- scat"$sample" + theme(legend.position="'"none"'")" >> RplottingLowCov.R
done

myvar=$(for sample in $sample_names ; do echo $sample; done | wc -l)
echo "myvar: "$myvar 
myvar2=$(echo "scale=1; $myvar/2+0.0005" | bc)
echo "myvar2: "$myvar2
nrow=$(printf "%.0f\n" $(echo $myvar2 | bc -l))
echo "nrow: "$nrow
widths=$(for ((i=1; i<=2; i++)); do printf "%0.s62"; done | sed 's/62/6, 2/g' | sed 's/26/2, 6/g')
echo "widths: "$widths
heights=$(for ((i=1; i<=$nrow; i++)); do printf "%0.s10"; done | sed 's/01/0, 1/g')
echo "heights: "$heights
col=$(if [[ $myvar == 1 ]] ; then echo "2"; else echo "4"; fi)
echo "col: "$col
#echo "grid.arrange("$insamp" legend, ncol = "$col", nrow = "$nrow", widths=c("$widths"), heights=c("$heights"))" >> RplottingLowCov.R
echo "grid.arrange("$insamp" ncol = "$col", nrow = "$nrow", widths=c("$widths"), heights=c("$heights"))" >> RplottingLowCov.R
echo -e "\n"

for sample in $sample_names
do
echo "stats"$sample" <- boxplot.stats(SingleCopy"$sample""'$'"coverage
              )"'$'"stats
write.table(stats"$sample",sep="'"\t"'",file="'"./stats'$sample'"'")

out"$sample" <- boxplot.stats(SingleCopy"$sample""'$'"coverage
              )"'$'"out
write.table(out"$sample",sep="'"\t"'",file="'"./out'$sample'"'")
"
done >> RplottingLowCov.R

######

echo "run RplottingLowCov.R"

R CMD BATCH RplottingLowCov.R
mv Rplots.pdf $sample.plots.pdf
#####

# step 6
# extract low coverage genes with annotation, count them, separate from hypotheticals
# required input files: stats$sample from step 5, CovPerGene_$sample\_annot.txt from step 4
echo "6. extract low coverage genes with annotation, count them, separate from hypotheticals for $sample_names"
for sample in $sample_names
do
echo "=> extracting low coverage genes from sample $sample"
lowcov=`awk '{if ($1=="\"1\"") print $2}' stats$sample`
awk -v var=$lowcov -F "\t" '{if($7<=var) print $0}' CovPerGene_$sample\_annot.txt | sed 's/ID=//g' > LowCovGenes_$sample
awk -v var=$lowcov -F "\t" '{if($7<=var) print $0}' CovPerGene_$sample\_annot.txt | sed 's/ID=//g' | grep -v "hypothetical protein" | sort | uniq > LowCovGenes_noHypo_$sample
done
## make amino acid fasta files for non-hypothetical low coverage outliers
for sample in $sample_names
do
echo "=> make amino acid fasta files for non-hypothetical low coverage outliers from sample $sample"
grep -v "hypothetical protein" LowCovGenes_$sample | awk '{print $1}' | grep ".peg." > LowCovGeneID_$sample
xargs samtools faidx $genome_files/$ref.faa < LowCovGeneID_$sample > LowCovGenes_$sample.fa
done
## write hypothetical and non-hypothetical low coverage outlier numbers to file
echo "Low coverage genes $sample" > LowCovCat.txt
for sample in $sample_names
do
echo "=> write hypothetical and non-hypothetical low coverage outlier numbers to file from sample $sample"
echo "$sample:" >> LowCovCat.txt
echo " " >> LowCovCat.txt 
echo "hypotheticals" >> LowCovCat.txt 
hyp=`grep "hypothetical protein" LowCovGenes_$sample | sort | uniq | wc -l`
echo $hyp >> LowCovCat.txt
echo "annotated" >> LowCovCat.txt 
hyp=`grep -v "hypothetical protein" LowCovGenes_$sample | sort | uniq | wc -l`
echo $hyp >> LowCovCat.txt
echo "------------" >> LowCovCat.txt
done

echo =================================================================
echo "Script finished $(date)"
DATE=`date +%d%m%Y`
echo ================================================================= 


