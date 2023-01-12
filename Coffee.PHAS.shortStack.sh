#screen -S 24PHAS
screen -r 24PHAS

cd /home/tcherubino/Carabica.smallRNA/PHAS

#sed 's/"//g' phasingPlusStrand.tab | sed '1d' |  awk '{print $2}' | awk -F":" '{print $1":"$2}' > getThose.txt

sed '1d' phasing.tab  |  awk '{if($2==24) print $1}'   > getThose.txt

#make multi line PHAS file single line

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ../DE/smallReference/PHAS.fa > sLPHAS.fa

grep -f getThose.txt -A 1 sLPHAS.fa > searchForSimilarity.fa

sed -i"" '/--/d' searchForSimilarity.fa

#/home/tcherubino/bin/mafft/mafft --globalpair --op 1.53 --ep 0.123 --leavegappyregion --thread 50 --nuc  searchForSimilarity.fa  > searchForSimilarity.gappy.fasta

#/home/tcherubino/bin/mafft/mafft --auto --op 1.53 --ep 0.123 --leavegappyregion --thread 50 --nuc  searchForSimilarity.fa  > searchForSimilarity.auto.gappy.fasta


#############################
##########MEME SUITE##########
################################

##############################################################
#will need to run in LFMP0 due to sudo restrictions on dashi##
#############################################################

#zoops	Zero or One Occurrence Per Sequence	MEME assumes that each sequence may contain at most one occurrence of each motif. This option is useful when you suspect that some motifs may be missing from some of the primary sequences. In that case, the motifs found will be more accurate than using the first option. This option takes more computer time than the first option (about twice as much) and is slightly less sensitive to weak motifs present in all of the primary sequences.


#~/meme/bin/meme searchForSimilarity.fa  -nmotifs 50 -minw 6 -maxw 26 -oc all.putative.24PHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 50

grep -f getTose.txt -A 1 PHAS.fa > 24-PHAS.fa


#20nt
~/meme/bin/meme searchForSimilarity.fa  -nmotifs 5 -w 20 -oc 20nt.all.putative.24PHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 50

tar -zcvf 20nt.all.putative.24PHAS.meme.out.tar.gz 20nt.all.putative.24PHAS.meme.out/

#21nt
~/meme/bin/meme searchForSimilarity.fa  -nmotifs 5 -w 21 -oc 21nt.all.putative.24PHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 50

tar -zcvf 21nt.all.putative.24PHAS.meme.out.tar.gz 21nt.all.putative.24PHAS.meme.out/

#22nt
~/meme/bin/meme searchForSimilarity.fa  -nmotifs 5 -w 22 -oc 22nt.all.putative.24PHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 50

tar -zcvf 22nt.all.putative.24PHAS.meme.out.tar.gz 22nt.all.putative.24PHAS.meme.out/

#23nt
~/meme/bin/meme searchForSimilarity.fa  -nmotifs 5 -w 23 -oc 23nt.all.putative.24PHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 50

tar -zcvf 23nt.all.putative.24PHAS.meme.out.tar.gz 23nt.all.putative.24PHAS.meme.out/


#24nt
~/meme/bin/meme searchForSimilarity.fa  -nmotifs 5 -w 24 -oc 24nt.all.putative.24PHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 50

tar -zcvf 24nt.all.putative.24PHAS.meme.out.tar.gz 24nt.all.putative.24PHAS.meme.out/


#####################################
###########onlyTheBest##############
###################################

#20nt
~/meme/bin/meme searchForSimilarity.fa  -nmotifs 1 -w 20 -oc 20.onlyTheBest.nt.all.putative.24PHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 50

tar -zcvf 20.onlyTheBest.nt.all.putative.24PHAS.meme.out.tar.gz 20.onlyTheBest.nt.all.putative.24PHAS.meme.out/

#21nt
~/meme/bin/meme searchForSimilarity.fa  -nmotifs 1 -w 21 -oc 21.onlyTheBest.nt.all.putative.24PHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 50

tar -zcvf 21.onlyTheBest.nt.all.putative.24PHAS.meme.out.tar.gz 21.onlyTheBest.nt.all.putative.24PHAS.meme.out/

#22nt
~/meme/bin/meme searchForSimilarity.fa  -nmotifs 1 -w 22 -oc 22.onlyTheBest.nt.all.putative.24PHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 50

tar -zcvf 22.onlyTheBest.nt.all.putative.24PHAS.meme.out.tar.gz 22.onlyTheBest.nt.all.putative.24PHAS.meme.out/

#23nt
~/meme/bin/meme searchForSimilarity.fa  -nmotifs 1 -w 23 -oc 23.onlyTheBest.nt.all.putative.24PHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 50

tar -zcvf 23.onlyTheBest.nt.all.putative.24PHAS.meme.out.tar.gz 23.onlyTheBest.nt.all.putative.24PHAS.meme.out/


#24nt
~/meme/bin/meme searchForSimilarity.fa  -nmotifs 1 -w 24 -oc 24.onlyTheBest.nt.all.putative.24PHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 50

tar -zcvf 24.onlyTheBest.nt.all.putative.24PHAS.meme.out.tar.gz 24.onlyTheBest.nt.all.putative.24PHAS.meme.out/

###################################################
##and a crazzy one to find miRNAs complementarity##
###################################################

~/meme/bin/meme searchForSimilarity.fa -nmotifs 10000 -minw 6 -maxw 26 -oc all.putative.24PHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 300

#tar -zcvf 20nt.all.putative.24PHAS.meme.out.tar.gz 20nt.all.putative.24PHAS.meme.out/


~/meme/bin/meme 24-PHAS.fa -nmotifs 1 -minw 20 -maxw 24 -oc all.24PHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 100

tar -zcvf all.24PHAS.meme.out.tar.gz all.24PHAS.meme.out


grep -f getPlusStrand.txt -A 1 24-PHAS.fa > plus24-phas.fa

grep -f getMinusStrand.txt -A 1 24-PHAS.fa > minus24-phas.fa

seqtk seq -r minus24-phas.fa > rev.comp.minus24-phas.fa


cat plus24-phas.fa rev.comp.minus24-phas.fa > processed24PHAS.fa

sed -i"" '/--/d' processed24PHAS.fa


~/meme/bin/meme processed24PHAS.fa -nmotifs 1 -minw 22 -maxw 24 -oc processed24PHAS.fa.meme.out -p 8 -dna -mod zoops -minsites 100

tar -zcvf processed24PHAS.fa.meme.out.tar.gz processed24PHAS.fa.meme.out


mafft --localpair --op 1.53 --ep 0.123 --thread 8 --nuc  processed24PHAS.fa > 24-phas.Alingment.fasta



##########################################################
############Alingment only to PHAS loci in the genome###################
##########################################################

cp ../DE/smallReference/PHAS.fa ./

grep ">" PHAS.fa | sed 's/>//g' | sed 's/PHAS_//g' > locus


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

#concatenate all fasta files
cat ../DE/*fasta > ./all.libs.concat.fasta



for file in `ls ../DE/*fasta`;
do
  ~/bin/ShortStack/ShortStack --g ../processedGenome/ToBlakeSystem/whole.C.a.formated.genome.fa --bowtie_cores 2 --readfile $file --locifile locus --mmap u --nohp --mismatches 0 --outdir `basename -s fasta $file`onlyPHAS.noMultinorMismatch.ShortStack &

done

mkdir shortStack

mv *ShortStack shortStack/


mkdir processBam
cd processBam
for folder in $(ls -d ../*onlyPHAS.noMultinorMismatch.ShortStack/)
do
    echo $folder
    cd $folder
    pwd

    samtools sort -@ 55 $(basename -s onlyPHAS.noMultinorMismatch.ShortStack $folder)bam > $(basename -s onlyPHAS.noMultinorMismatch.ShortStack $folder)sorted.bam

#keep only mapped

    samtools view -b --threads 55 -F 4 $(basename -s onlyPHAS.noMultinorMismatch.ShortStack $folder)sorted.bam > $(basename -s onlyPHAS.noMultinorMismatch.ShortStack $folder)sorted.Mapped.bam

    rm $(basename -s onlyPHAS.noMultinorMismatch.ShortStack $folder)sorted.bam

    cp $(basename -s onlyPHAS.noMultinorMismatch.ShortStack $folder)sorted.Mapped.bam ../processBam


done

cd ../processBam

samtools merge -@ 55 wholeGenome.merged.bam *sorted.Mapped.bam


samtools index -@ 55 wholeGenome.merged.bam

~/bin/ShortStack/ShortStack --g ../../../processedGenome/ToBlakeSystem/whole.C.a.formated.genome.fa --bowtie_cores 60 --bamfile wholeGenome.merged.bam --locifile locus --mmap u --nohp --mismatches 0 --outdir all.bam.PHAS


#####################################
#meanwyle...#########################
#####################################

awk -F":|-" '{print $1"\t"$2"\t"$3}' locus > locus.bed

#repetitive regions
 bedtools intersect -a locus.bed -b ../processedGenome/RepeatMasker/whole.C.a.formated.genome.fa.RepeatMasker.out.gff  -loj > intersected.with.repRegion.tab

#protein coding regions

grep gene ../processedGenome/ToBlakeSystem/C.a.ProteinCoding.genes.gff > gene.gff

bedtools intersect -a locus.bed -b gene.gff -loj > intersected.with.ProteinCodign.tab


bedtools intersect -a phasis21.bed -b shortstack21.bed > intersected.21.phasis.shortStack.tab
#140 phasis results overlap with 185 shortStack

####################################
######### some R analysis###########
####################################

data <- read.delim("TE.andSK.Results.txt")

sig.data <- data[data$PhaseScore>15,]

sig.num.data <- sig.data[,c(1,3,4,5,6,7,11,12,14,16,17),]

sig.num.data$is.TE <- factor(sig.num.data$is.TE)

plot(PhaseScore ~ is.TE,data=sig.num.data)

boxplot(PhaseScore ~ is.TE,data=sig.num.data)

plot(density(sig.num.data$PhaseScore))

t.test(PhaseScore ~ is.TE,data=sig.num.data)

plot(Length ~ Reads*is.TE,data=sig.num.data)

plot(density(sig.num.data$Length))

plot(hist(sig.num.data$Length))

plot(Length ~is.TE,data=sig.num.data)

t.test(Length ~is.TE,data=sig.num.data, alternative="two.sided", var.equal=FALSE)

############################
plot(Reads ~is.TE,data=sig.num.data)

t.test(Reads ~is.TE,data=sig.num.data, alternative="two.sided", var.equal=FALSE)

#############################
plot(RPM ~is.TE,data=sig.num.data)

t.test(RPM ~is.TE,data=sig.num.data, alternative="two.sided", var.equal=FALSE)

#############################
plot(Complexity ~is.TE,data=sig.num.data)

t.test(Complexity ~is.TE,data=sig.num.data, alternative="two.sided", var.equal=FALSE)

#############################
plot(UniqueReads ~is.TE,data=sig.num.data)

t.test(UniqueReads ~is.TE,data=sig.num.data, alternative="two.sided", var.equal=FALSE)
#############################
plot(FracTop ~is.TE,data=sig.num.data)

t.test(FracTop ~is.TE,data=sig.num.data, alternative="two.sided", var.equal=FALSE)

#############################

TEand21 <- sig.num.data[sig.num.data$is.TE == 1 & sig.num.data$DicerCall ==21,]

nrow(TEand21)

TEand24 <- sig.num.data[sig.num.data$is.TE == 1 & sig.num.data$DicerCall ==24,]

nrow(TEand24)

#######################################
sig21 <- sig.data[sig.data$DicerCall==21,]


t.test(PhaseScore ~ is.TE,data=sig21)

boxplot(PhaseScore ~ is.TE,data=sig21)

plot(density(sig21$PhaseScore))

boxplot(Length ~ is.TE,data=sig21)

plot(density(sig21$Length))

plot(hist(sig21$Length))

boxplot(Length ~is.TE,data=sig21)

t.test(Length ~is.TE,data=sig21, alternative="two.sided", var.equal=FALSE)

############################
boxplot(Reads ~is.TE,data=sig21)

t.test(Reads ~is.TE,data=sig21, alternative="two.sided", var.equal=FALSE)

#############################
boxplot(RPM ~is.TE,data=sig21)

t.test(RPM ~is.TE,data=sig21, alternative="two.sided", var.equal=FALSE)

#############################
boxplot(Complexity ~is.TE,data=sig21)

t.test(Complexity ~is.TE,data=sig21, alternative="two.sided", var.equal=FALSE)

#############################
boxplot(UniqueReads ~is.TE,data=sig21)

t.test(UniqueReads ~is.TE,data=sig21, alternative="two.sided", var.equal=FALSE)
#############################
boxplot(FracTop ~is.TE,data=sig21)

t.test(FracTop ~is.TE,data=sig21, alternative="two.sided", var.equal=FALSE)

#######################################
sig24 <- sig.data[sig.data$DicerCall==24,]

t.test(PhaseScore ~ is.TE,data=sig24)

boxplot(PhaseScore ~ is.TE,data=sig24)

plot(density(sig24$PhaseScore))

boxplot(Length ~ is.TE,data=sig24)

plot(density(sig24$Length))

plot(hist(sig24$Length))

boxplot(Length ~is.TE,data=sig24)

t.test(Length ~is.TE,data=sig24, alternative="two.sided", var.equal=FALSE)

############################
boxplot(Reads ~is.TE,data=sig24)

t.test(Reads ~is.TE,data=sig24, alternative="two.sided", var.equal=FALSE)

#############################
boxplot(RPM ~is.TE,data=sig24)

t.test(RPM ~is.TE,data=sig24, alternative="two.sided", var.equal=FALSE)

#############################
boxplot(Complexity ~is.TE,data=sig24)

t.test(Complexity ~is.TE,data=sig24, alternative="two.sided", var.equal=FALSE)

#############################
boxplot(UniqueReads ~is.TE,data=sig24)

t.test(UniqueReads ~is.TE,data=sig24, alternative="two.sided", var.equal=FALSE)
#############################
boxplot(FracTop ~is.TE,data=sig24)

t.test(FracTop ~is.TE,data=sig24, alternative="two.sided", var.equal=FALSE)


##########################
##########PCD PHAS########
##########################

data.PCD <- read.delim("PCD.andSK.Results.txt")

sig.data.PCD <- data.PCD[data.PCD$PhaseScore>15,]

sig.num.data.PCD <- sig.data.PCD[,c(1,3,4,5,6,7,11,12,14,16,17),]

sig.num.data.PCD$is.PCD <- factor(sig.num.data.PCD$is.PCD)

plot(PhaseScore ~ is.PCD,data=sig.num.data.PCD)

boxplot(PhaseScore ~ is.PCD,data=sig.num.data.PCD)

plot(density(sig.num.data.PCD$PhaseScore))

t.test(PhaseScore ~ is.PCD,data=sig.num.data.PCD)

boxplot(Length ~ Reads*is.PCD,data=sig.num.data.PCD)

plot(density(sig.num.data.PCD$Length))

plot(hist(sig.num.data.PCD$Length))

plot(Length ~is.PCD,data=sig.num.data.PCD)

t.test(Length ~is.PCD,data=sig.num.data.PCD, alternative="two.sided", var.equal=FALSE)

############################
plot(Reads ~is.PCD,data=sig.num.data.PCD)

t.test(Reads ~is.PCD,data=sig.num.data.PCD, alternative="two.sided", var.equal=FALSE)

#############################
plot(RPM ~is.PCD,data=sig.num.data.PCD)

t.test(RPM ~is.PCD,data=sig.num.data.PCD, alternative="two.sided", var.equal=FALSE)

#############################
plot(Complexity ~is.PCD,data=sig.num.data.PCD)

t.test(Complexity ~is.PCD,data=sig.num.data.PCD, alternative="two.sided", var.equal=FALSE)

#############################
plot(UniqueReads ~is.PCD,data=sig.num.data.PCD)

t.test(UniqueReads ~is.PCD,data=sig.num.data.PCD, alternative="two.sided", var.equal=FALSE)
#############################
plot(FracTop ~is.PCD,data=sig.num.data.PCD)

t.test(FracTop ~is.PCD,data=sig.num.data.PCD, alternative="two.sided", var.equal=FALSE)

#############################


PCDand21 <- sig.num.data.PCD[sig.num.data.PCD$is.PCD == 1 & sig.num.data.PCD$DicerCall ==21,]

nrow(PCDand21)

PCDand24 <- sig.num.data.PCD[sig.num.data.PCD$is.PCD == 1 & sig.num.data.PCD$DicerCall ==24,]

nrow(PCDand24)

#######################################
sig21 <- sig.data.PCD[sig.data.PCD$DicerCall==21,]


t.test(PhaseScore ~ is.PCD,data=sig21)

plot(PhaseScore ~ is.PCD,data=sig21)

boxplot(PhaseScore ~ is.PCD,data=sig21)

plot(density(sig21$PhaseScore))

boxplot(Length ~ is.PCD,data=sig21)

plot(density(sig21$Length))

plot(hist(sig21$Length))

boxplot(Length ~is.PCD,data=sig21)

t.test(Length ~is.PCD,data=sig21, alternative="two.sided", var.equal=FALSE)

############################
boxplot(Reads ~is.PCD,data=sig21)

t.test(Reads ~is.PCD,data=sig21, alternative="two.sided", var.equal=FALSE)

#############################
boxplot(RPM ~is.PCD,data=sig21)

t.test(RPM ~is.PCD,data=sig21, alternative="two.sided", var.equal=FALSE)

#############################
boxplot(Complexity ~is.PCD,data=sig21)

t.test(Complexity ~is.PCD,data=sig21, alternative="two.sided", var.equal=FALSE)

#############################
boxplot(UniqueReads ~is.PCD,data=sig21)

t.test(UniqueReads ~is.PCD,data=sig21, alternative="two.sided", var.equal=FALSE)
#############################
boxplot(FracTop ~is.PCD,data=sig21)

t.test(FracTop ~is.PCD,data=sig21, alternative="two.sided", var.equal=FALSE)

#######################################
sig24 <- sig.data.PCD[sig.data.PCD$DicerCall==24,]

t.test(PhaseScore ~ is.PCD,data=sig24)

boxplot(PhaseScore ~ is.PCD,data=sig24)

plot(density(sig24$PhaseScore))

boxplot(Length ~ is.PCD,data=sig24)

plot(density(sig24$Length))

plot(hist(sig24$Length))

boxplot(Length ~is.PCD,data=sig24)

t.test(Length ~is.PCD,data=sig24, alternative="two.sided", var.equal=FALSE)

############################
boxplot(Reads ~is.PCD,data=sig24)

t.test(Reads ~is.PCD,data=sig24, alternative="two.sided", var.equal=FALSE)

#############################
boxplot(RPM ~is.PCD,data=sig24)

t.test(RPM ~is.PCD,data=sig24, alternative="two.sided", var.equal=FALSE)

#############################
boxplot(Complexity ~is.PCD,data=sig24)

t.test(Complexity ~is.PCD,data=sig24, alternative="two.sided", var.equal=FALSE)

#############################
boxplot(UniqueReads ~is.PCD,data=sig24)

t.test(UniqueReads ~is.PCD,data=sig24, alternative="two.sided", var.equal=FALSE)
#############################
boxplot(FracTop ~is.PCD,data=sig24)

t.test(FracTop ~is.PCD,data=sig24, alternative="two.sided", var.equal=FALSE)


#############################################

TEloci <-    as.character(sig.num.data$X.Locus[sig.num.data$is.TE ==1 ])

notTEloci <- as.character(sig.num.data$X.Locus[sig.num.data$is.TE !=1 ])

PCDloci <- as.character(sig.num.data.PCD$X.Locus[sig.num.data.PCD$is.PCD ==1 ])

notPCDloci <- as.character(sig.num.data.PCD$X.Locus[sig.num.data.PCD$is.PCD !=1 ])

overlapTE <-  sig.num.data[ as.character(sig.num.data$X.Locus) %in% intersect(TEloci, PCDloci),]

overlapPCD <- sig.num.data.PCD[ as.character(sig.num.data.PCD$X.Locus) %in% intersect(TEloci, PCDloci),]


#not Genic and TE
length(intersect(notTEloci,notPCDloci))
################################################
#For 21PHAS

TE21loci <-    as.character(sig.num.data$X.Locus[sig.num.data$is.TE ==1 &  sig.num.data$DicerCall == 21])

notTE21loci <- as.character(sig.num.data$X.Locus[sig.num.data$is.TE !=1 &  sig.num.data$DicerCall == 21])

PCD21loci <- as.character(sig.num.data.PCD$X.Locus[sig.num.data.PCD$is.PCD ==1 &  sig.num.data.PCD$DicerCall == 21])

notPCD21loci <- as.character(sig.num.data.PCD$X.Locus[sig.num.data.PCD$is.PCD !=1 & sig.num.data.PCD$DicerCall == 21])

overlapTE21 <-  sig.num.data[ as.character(sig.num.data$X.Locus) %in% intersect(TE21loci, PCD21loci),]

overlapPCD21 <- sig.num.data.PCD[ as.character(sig.num.data.PCD$X.Locus) %in% intersect(TE21loci, PCD21loci),]


#not Genic and TE
length(intersect(notTE21loci,notPCD21loci))

################################################
#For 24PHAS

TE24loci <-    as.character(sig.num.data$X.Locus[sig.num.data$is.TE ==1 &  sig.num.data$DicerCall == 24])

notTE24loci <- as.character(sig.num.data$X.Locus[sig.num.data$is.TE !=1 &  sig.num.data$DicerCall == 24])

PCD24loci <- as.character(sig.num.data.PCD$X.Locus[sig.num.data.PCD$is.PCD ==1 &  sig.num.data.PCD$DicerCall == 24])

notPCD24loci <- as.character(sig.num.data.PCD$X.Locus[sig.num.data.PCD$is.PCD !=1 & sig.num.data.PCD$DicerCall == 24])

overlapTE24 <-  sig.num.data[ as.character(sig.num.data$X.Locus) %in% intersect(TE24loci, PCD24loci),]

overlapPCD24 <- sig.num.data.PCD[ as.character(sig.num.data.PCD$X.Locus) %in% intersect(TE24loci, PCD24loci),]


#not Genic and TE
length(intersect(notTE24loci,notPCD24loci))


################################################
geneDescriptions <- read.delim("../../doctors/C.arabica.blake.system/ToBlakeSystem/C.arabica.geneDescriptions.tab",,header =F)

colnames(geneDescriptions) <- c("geneID","TranscriptID","Description","GO","EC")

#PCDand21

 PCDand21.descriptions <- geneDescriptions[geneDescriptions$geneID %in% PCDand21$PCD,]

grep("resistance",PCDand21.descriptions$Description )

length(grep("resistance",PCDand21.descriptions$Description ))

#127

length(grep("pentatricopeptide",PCDand21.descriptions$Description ))

#15

length(grep("--NA--",PCDand21.descriptions$Description ))


length(grep("Dicer",PCDand21.descriptions$Description ))

#2

#PCDand24

PCDand24.descriptions <- geneDescriptions[geneDescriptions$geneID %in% PCDand24$PCD,]

grep("resistance",PCDand24.descriptions$Description )

length(grep("resistance",PCDand24.descriptions$Description ))

#uncharacterized
length(grep("uncharacterized",PCDand24.descriptions$Description ))

length(grep("pentatricopeptide",PCDand24.descriptions$Description ))

#3

length(grep("--NA--",PCDand24.descriptions$Description ))


#############################################
#Overlap with TE.

allTEdescriptions <- read.delim("intersected.with.repRegion.tab",header=F)

allTEdescriptions<- allTEdescriptions[,c(1,2,3,12,13)]

colnames(allTEdescriptions) <- c("chr","start","stop","description","ID")

##############
###TE21loci###
##############

 TE21loci.descriptions <- sig.num.data$TE[sig.num.data$X.Locus %in% TE21loci]

TE21loci.allOccur.descriptions <- allTEdescriptions$description[allTEdescriptions$ID %in% TE21loci]

#Helitron

length(grep("Helitron", TE21loci.descriptions))

length(grep("Helitron", TE21loci.allOccur.descriptions))

#Gypsy
length(grep("Gypsy", TE21loci.descriptions))

length(grep("Gypsy", TE21loci.allOccur.descriptions))

#Copia
length(grep("Copia", TE21loci.descriptions))

length(grep("Copia", TE21loci.allOccur.descriptions))

#hAT

length(grep("hAT", TE21loci.descriptions))

length(grep("hAT", TE21loci.allOccur.descriptions))


##############
###TE24loci###
##############

TE24loci.descriptions <- sig.num.data$TE[sig.num.data$X.Locus %in% TE24loci]

TE24loci.allOccur.descriptions <- allTEdescriptions$description[allTEdescriptions$ID %in% TE24loci]

#Helitron

length(grep("Helitron", TE24loci.descriptions))

length(grep("Helitron", TE24loci.allOccur.descriptions))

#Gypsy
length(grep("Gypsy", TE24loci.descriptions))

length(grep("Gypsy", TE24loci.allOccur.descriptions))

#Copia
length(grep("Copia", TE24loci.descriptions))

length(grep("Copia", TE24loci.allOccur.descriptions))

#hAT

length(grep("hAT", TE24loci.descriptions))

length(grep("hAT", TE24loci.allOccur.descriptions))


#MuDR
length(grep("MuDR", TE24loci.descriptions))

length(grep("MuDR", TE24loci.allOccur.descriptions))

#HARB
length(grep("HARB", TE24loci.descriptions))

length(grep("HARB", TE24loci.allOccur.descriptions))

#Ogre
length(grep("Ogre", TE24loci.descriptions))

length(grep("Ogre", TE24loci.allOccur.descriptions))

#AalpV

length(grep("AalpV", TE24loci.descriptions))

length(grep("AalpV", TE24loci.allOccur.descriptions))

#21PHAS overlapping Protein Coding Genes (PCD) and TEs (special focus here)

write.table(as.data.frame(PCDand21$X.Locus[PCDand21$X.Locus %in% overlapTE21$X.Locus]),row.names=F,col.names=F,file="PCD.TE.21sigPHAS.tab")

#21PHAS overlapping only PCD
`%nin%` = Negate(`%in%`)

write.table(as.data.frame(PCDand21$X.Locus[PCDand21$X.Locus %nin% overlapTE21$X.Locus]),row.names=F,col.names=F,file="PCD.21sigPHAS.tab")

#24PHAS overlapping PCD and TEs
write.table(as.data.frame(PCDand24$X.Locus[PCDand24$X.Locus %in% overlapTE24$X.Locus]),row.names=F,col.names=F,file="PCD.TE.24sigPHAS.tab")


#24PHAS overlapping only PCD
write.table(as.data.frame(PCDand24$X.Locus[PCDand24$X.Locus %nin% overlapTE24$X.Locus]),row.names=F,col.names=F,file="PCD.24sigPHAS.tab")


#24PHAS overlapping TE and not PCD

write.table(as.data.frame(
TEand24$X.Locus[TEand24$X.Locus %nin% overlapPCD24$X.Locus]),
row.names=F,
col.names=F,
file="TE.24sigPHAS.tab")


24PHAS overlapping neither TE and PCD

write.table(as.data.frame(intersect(notTE24loci,notPCD24loci)),row.names=F,col.names=F,file="not.PCD.TE.24sigPHAS.tab")


#rownames(comparissonsOfInterestExpressedUP)<- c("sn-RNA","sno-RNA","t-RNA","21-PHAS","24-PHAS","L24P-like","UNK","mi-RNA")

pastelColors <- c(
"#C0C8CE",
"#CAB393",
"#BF7F58",
"#5F6D54"
)

#UNKtype <- as.vector(unlist(read.table(pipe("pbpaste"), sep="\t", header=FALSE)))

#UNK <- data[data$X.Locus %in% UNKtype,]

#write.table(as.data.frame(UNK),row.names=F,col.names=F,file="UNKtype.tab")

overallDATA <- read.delim("overallPHASDATA.txt")

uniqueProportion<- overallDATA$UniqueReads/overallDATA$Reads
overallDATA<- cbind(overallDATA,uniqueProportion )

library("vioplot")
pdf("LengthXtype.pdf")
vioplot(Length ~ type,data=overallDATA,col=pastelColors,horizontal = TRUE)
dev.off()

pdf("uniqueProportionXtype.pdf")
vioplot(uniqueProportion ~ type,data=overallDATA,col=pastelColors,horizontal = TRUE)
dev.off()


pdf("UniqueReadsXtype.pdf")
vioplot(UniqueReads ~ type,data=overallDATA,col=pastelColors,horizontal = TRUE)
dev.off()


pdf("FracTopXtype.pdf")
vioplot(FracTop ~ type,data=overallDATA,col=pastelColors,horizontal = TRUE)
dev.off()


pdf("ComplexityXtype.pdf")
vioplot(Complexity ~ type,data=overallDATA,col=pastelColors,horizontal = TRUE)
dev.off()


pdf("PhaseScoreXtype.pdf")
vioplot(PhaseScore ~ type,data=overallDATA,col=pastelColors,horizontal = TRUE)
dev.off()

pdf("RPMXtype.pdf")
vioplot(log(RPM+1,10) ~ type,data=overallDATA,col=pastelColors,horizontal = TRUE)
dev.off()

UniqueProportion21phas <- overallDATA$uniqueProportion[overallDATA$type == "21-PHAS"]

UniqueProportion24phas <- overallDATA$uniqueProportion[overallDATA$type == "24-PHAS"]

UniqueProportion24Plikephas <- overallDATA$uniqueProportion[overallDATA$type == "24P-like"]

UniqueProportionUNK <- overallDATA$uniqueProportion[overallDATA$type == "UNK"]


TypeMeans <- c(mean(UniqueProportion21phas),
mean(UniqueProportion24phas),
mean(UniqueProportion24Plikephas),
mean(UniqueProportionUNK,na.rm=T))

TypeSD <- c(sd(UniqueProportion21phas),
sd(UniqueProportion24phas),
sd(UniqueProportion24Plikephas),
sd(UniqueProportionUNK,na.rm=T)
)
pdf("UniqueProportionBARPLOT.pdf")

barplot <- barplot(TypeMeans, col=pastelColors,ylim=c(0,1),ylab= "unique mapping rate", names= c("21-PHAS","24-PHAS","24P-like","UNK"))

arrows(x0 = barplot, y0 = TypeMeans - TypeSD, x1 = barplot, y1=TypeMeans + TypeSD ,code=3,angle=90,length=0.05,col="black",lwd=4)

dev.off()

save.image("cofeePhas.RData")


#######################################

##############Search for conserved motifs

sed -i"" 's/NC_039898.1/1/g'  ../PHAS.fa

sed -i"" 's/NC_039899.1/2/g'  ../PHAS.fa

sed -i"" 's/NC_039900.1/3/g'  ../PHAS.fa

sed -i"" 's/NC_039901.1/4/g'  ../PHAS.fa

sed -i"" 's/NC_039902.1/5/g'  ../PHAS.fa

sed -i"" 's/NC_039903.1/6/g'  ../PHAS.fa

sed -i"" 's/NC_039904.1/7/g'  ../PHAS.fa

sed -i"" 's/NC_039905.1/8/g'  ../PHAS.fa

sed -i"" 's/NC_039906.1/9/g'  ../PHAS.fa

sed -i"" 's/NC_039907.1/10/g'  ../PHAS.fa

sed -i"" 's/NC_039908.1/11/g'  ../PHAS.fa

sed -i"" 's/NC_039909.1/12/g'  ../PHAS.fa

sed -i"" 's/NC_039910.1/13/g'  ../PHAS.fa

sed -i"" 's/NC_039911.1/14/g'  ../PHAS.fa

sed -i"" 's/NC_039912.1/15/g'  ../PHAS.fa

sed -i"" 's/NC_039913.1/16/g'  ../PHAS.fa

sed -i"" 's/NC_039914.1/17/g'  ../PHAS.fa

sed -i"" 's/NC_039915.1/18/g'  ../PHAS.fa

sed -i"" 's/NC_039916.1/19/g'  ../PHAS.fa

sed -i"" 's/NC_039917.1/20/g'  ../PHAS.fa

sed -i"" 's/NC_039918.1/21/g'  ../PHAS.fa

sed -i"" 's/NC_039919.1/22/g'  ../PHAS.fa

#zoops	Zero or One Occurrence Per Sequence	MEME assumes that each sequence may contain at most one occurrence of each motif. This option is useful when you suspect that some motifs may be missing from some of the primary sequences. In that case, the motifs found will be more accurate than using the first option. This option takes more computer time than the first option (about twice as much) and is slightly less sensitive to weak motifs present in all of the primary sequences.

#PCD.TE.21sigPHAS

cd /Carabica/PHAS/PCD.TE.21sigPHAS
sed -i"" 's/"//g' PCD.TE.21sigPHAS.tab

grep -A 1 -f PCD.TE.21sigPHAS.tab ../PHAS.fa > PCD.TE.21sigPHAS.fasta

sed -i"" '/^-/d' PCD.TE.21sigPHAS.fasta

~/meme/bin/meme PCD.TE.21sigPHAS.fasta  -nmotifs 5 -minw 17 -maxw 24 -oc PCD.TE.21sigPHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 40

~/meme/bin/meme PCD.TE.21sigPHAS.fasta  -nmotifs 1 -minw 17 -maxw 24 -oc best.PCD.TE.21sigPHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 40

~/meme/bin/meme PCD.TE.21sigPHAS.fasta  -nmotifs 5 -w 18 -oc 18.PCD.TE.21sigPHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 40


#PCD.21sigPHAS

cd /Carabica/PHAS/PCD.21sigPHAS
sed -i"" 's/"//g' PCD.21sigPHAS.tab

grep -A 1 -f PCD.21sigPHAS.tab ../PHAS.fa > PCD.21sigPHAS.fasta

sed -i"" '/^-/d' PCD.21sigPHAS.fasta

~/meme/bin/meme PCD.21sigPHAS.fasta  -nmotifs 5 -minw 17 -maxw 24 -oc PCD.21sigPHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 40

~/meme/bin/meme PCD.21sigPHAS.fasta  -nmotifs 1 -minw 17 -maxw 24 -oc best.PCD.21sigPHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 40

~/meme/bin/meme PCD.21sigPHAS.fasta  -nmotifs 5 -w 18 -oc 18.PCD.21sigPHAS.meme.out -p 8 -dna -revcomp -mod zoops -minsites 40

#######################################
###########Run PHASIS#################
####################################

#/newdata/data1/homes/tcherubino/Carabica.smallRNA/PHAS

mkdir phasis
cd phasis

cat ../../smallRNAsLibraries/*txt > all.libs.Coffee.tab

sort all.libs.Coffee.tab > sorted.all.libs.Coffee.tab

rm all.libs.Coffee.tab

######################################################
#R

df <- read.delim("sorted.all.libs.Coffee.tab",header=F)
summedFragmentsAcrossLibs <- sapply(split.default(df$V2, df$V1), sum)

write.table(file = "summed.All.Libs.Coffee.tab",as.matrix(summedFragmentsAcrossLibs,ncol=2),sep="\t")

q()

 rm sorted.all.libs.Coffee.tab

sed -i"" 's/"//g' summed.All.Libs.Coffee.tab

#'21 PHASING
python ~/bin/phasis/phasis-core.v1.2.py

find . -type f -empty -print -delete


#24 PHASING
python ~/bin/phasis/phasis-core.v1.2.py

find . -type f -empty -print -delete

########################

mkdir 24PHAS
mv *sRNA_24.cluster ./24PHAS

mkdir 21PHAS

mv *sRNA_21.cluster ./21PHAS

cd 21PHAS

 python3 ~/bin/phasis/phasExplorer.py -mode basic -reference 21_PHAS_COFFE_SHORTSTAK.tab -clusters *sRNA_21.cluster -cutoffs 10,5  -pval 0.05 -phasing 21 -strandmode both

#Strange... No overlap found.

cd 24PHAS

python3 ~/bin/phasis/phasExplorer.py -mode basic -reference 24_PHAS_COFFE_SHORTSTAK.tab -clusters *sRNA_24.cluster -cutoffs 10,5  -pval 0.05 -phasing 24 -strandmode both
