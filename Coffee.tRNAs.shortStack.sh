#tRNA study
#screen -S tRNAs
screen -r tRNAs

mkdir tRNA
cd tRNA

cp ../DE/smallReference/Carabica.tRNA.fasta ./


grep ">" Carabica.tRNA.fasta | sed 's/>//g' | sed 's/TRNA_//g' > locus


sed -i"" 's/NC_039898.1/1/g' locus

sed -i"" 's/NC_039899.1/2/g' locus

sed -i"" 's/NC_039900.1/3/g' locus

sed -i"" 's/NC_039901.1/4/g' locus

sed -i"" 's/NC_039902.1/5/g' locus

sed -i"" 's/NC_039903.1/6/g' locus

sed -i"" 's/NC_039904.1/7/g' locus

sed -i"" 's/NC_039905.1/8/g' locus

sed -i"" 's/NC_039906.1/9/g' locus

sed -i"" 's/NC_039907.1/10/g' locus

sed -i"" 's/NC_039908.1/11/g' locus

sed -i"" 's/NC_039909.1/12/g' locus

sed -i"" 's/NC_039910.1/13/g' locus

sed -i"" 's/NC_039911.1/14/g' locus

sed -i"" 's/NC_039912.1/15/g' locus

sed -i"" 's/NC_039913.1/16/g' locus

sed -i"" 's/NC_039914.1/17/g' locus

sed -i"" 's/NC_039915.1/18/g' locus

sed -i"" 's/NC_039916.1/19/g' locus

sed -i"" 's/NC_039917.1/20/g' locus

sed -i"" 's/NC_039918.1/21/g' locus

sed -i"" 's/NC_039919.1/22/g' locus

cp locus bkup.locus

awk -F':|-' '{print $1":"$2+1"-"$3}' bkup.locus > temp
mv temp bkup.locus
#set the padding to 100

 awk -F':|-' '{print $1":"$2-100"-"$3+100}' locus > t

mv t locus


#instituir um foor loop aqui
~/bin/ShortStack/ShortStack --g ../processedGenome/ToBlakeSystem/whole.C.a.formated.genome.fa --bowtie_cores 60 --readfile ../PHAS/all.libs.concat.fasta --locifile locus --mmap u --nohp --mismatches 0 --outdir onlytRNA.noMultinorMismatch.ShortStack


##################################################
#format the original GFF according our annotation#
##################################################
cd /home/tcherubino/Carabica.smallRNA/processedGenome/originalGFF.NCBI

sed -i"" '/#/d' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039898.1/1/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039899.1/2/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039900.1/3/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039901.1/4/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039902.1/5/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039903.1/6/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039904.1/7/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039905.1/8/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039906.1/9/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039907.1/10/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039908.1/11/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039909.1/12/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039910.1/13/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039911.1/14/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039912.1/15/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039913.1/16/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039914.1/17/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039915.1/18/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039916.1/19/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039917.1/20/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039918.1/21/g' proccessed1_Cara_1.0_genomic.gff

sed -i"" 's/NC_039919.1/22/g' proccessed1_Cara_1.0_genomic.gff

awk '{
  if ($3 == "tRNA")
    print $0;
  }' proccessed1_Cara_1.0_genomic.gff > tRNA.gff

#get the tRNA coordinates from gff

cp ../processedGenome/originalGFF.NCBI/tRNA.gff ./orignal_tRNA.gff

awk '{print $1":"$4"-"$5"\t"$9}' orignal_tRNA.gff > tRNA_annot.tab

grep -f bkup.locus tRNA_annot.tab > unique.tRNA_annot.tab


##########################
######Strategy############
##########################

#alging all reads to the collapsed tRNA fasta file
bowtie-build Carabica.tRNA.fasta Carabica.tRNA

cat ../DE/*fasta > all.reads.fasta


#To alginm to the tRNA I won't allow mismatches (-v 0) and I want only yo know reads aligned or not
bowtie -x Carabica.tRNA -f all.reads.fasta -y -k 1 -v 0 -p 60 --no-unal -S > bowtie.result.tRNA.sam



# reads processed: 1233211059
# reads with at least one alignment: 21,855,655 (1.77%)
# reads that failed to align: 1211355404 (98.23%)
#Reported 21855655 alignments

/home/tcherubino/bin/FastQC/fastqc bowtie.result.tRNA.sam

samtools sort -@ 60 bowtie.result.tRNA.sam > bowtie.result.tRNA.firstRound.sorted.bam
rm bowtie.result.tRNA.sam


samtools index bowtie.result.tRNA.firstRound.sorted.bam

#retrieve the reads mapped only to tRNAs
samtools fasta -@ 60 -F 4 -t bowtie.result.tRNA.firstRound.sorted.bam > mappedTotRNAS.fasta

sed '/^>/d' mappedTotRNAS.fasta | sort | uniq -c > sequence_count

awk '{print $2"\t"$1}' sequence_count > mappedTotRNAS.tag

###########################
R

data <- read.delim("mappedTotRNAS.tag",header=F)

pdf("readsMappedTotRNAs.pdf")
plot(hist(data$V2))
plot(density(data$V2))
dev.off()

#Count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor. Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)

scallinFactor <- sum(data$V2)/1000000

data$V3 <- data$V2/scallinFactor

pdf("readsMappedTotRNAsCPM.pdf")
plot(hist(data$V3))
plot(density(data$V3))
dev.off()

pdf("sequencesLengths.pdf")
plot(hist(nchar(data$V1)))
plot(density(nchar(data$V1)))
dev.off()


#filter sequences larger than 24nt
SizeTreshhold <- data[nchar(data$V1)<= 24,]


#lets filter reads with CPM less than 1

CPM.1.treshold <- SizeTreshhold[SizeTreshhold$V3>1,]

CPM.1.treshold$V4 <- paste0(nchar(CPM.1.treshold$V1),"_",CPM.1.treshold$V1,"_",CPM.1.treshold$V2)


CPM.1 <- cbind(CPM.1.treshold$V4,CPM.1.treshold$V1)


write.table(CPM.1, file="processedReadsMappedTotTNAs.tab",row.names=F,col.names=F,sep="\t",quote=F)


awk -F'\t' '{printf ">%s\n%s\n",$1,$2}' processedReadsMappedTotTNAs.tab > mappedTotRNAPrecursors.fasta
"


bowtie -x ~/Carabica.smallRNA/processedGenome/ToBlakeSystem/whole.C.a.formated.genome -f ./mappedTotRNAPrecursors.fasta -y -a -v 3 -p 60 -S > bowtie.result.tRNA.pre.sam



# reads processed: 7458
# reads with at least one alignment: 7458 (100.00%)
# reads that failed to align: 0 (0.00%)
Reported 1,062,296 alignments
"

samtools sort -@ 60 bowtie.result.tRNA.pre.sam > bowtie.result.t_RNA.3miss.sorted.bam
"
samtools index -@ 60 bowtie.result.t_RNA.3miss.sorted.bam

samtools coverage bowtie.result.t_RNA.3miss.sorted.bam -m -o bowtie.result.t_RNA.3miss.sorted.coverage



##########################
########sPARTA on tRNAs###
##########################
mkdir sPARTA

cd sPARTA
pwd
#/home/tcherubino/Carabica.smallRNA/tRNA/sPARTA

#GENIC
python3 ~/bin/sPARTA/sPARTA.py -genomeFile /home/tcherubino/Carabica.smallRNA/processedGenome/whole.C.a.formated.genome.fa -annoType GFF -annoFile final.gff -genomeFeature 0 -miRNAFile  mappedTotRNAPrecursors.fasta -libs *chopped.txt -tarPred H -tarScore N --tag2FASTA --map2DD --validate -accel 60 -minTagLen 18

mv output/ tRNAS.genic.output

#INTERGENIC
python3 ~/bin/sPARTA/sPARTA.py -genomeFile /home/tcherubino/Carabica.smallRNA/processedGenome/whole.C.a.formated.genome.fa -annoType GFF -annoFile final.gff -genomeFeature 1 -miRNAFile mappedTotRNAPrecursors.fasta -libs *chopped.txt -tarPred H -tarScore N --tag2FASTA --map2DD -accel 60 -minTagLen 20 --validate


while read -r line;
do
grep -w $line C.arabica.geneDescriptions.tab >> tRNAtargetsDescription.tab
done < targets
