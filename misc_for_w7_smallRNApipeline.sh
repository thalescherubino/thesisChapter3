#screen -S misc_w7

screen -r misc_w7

samtools faidx allSmall.fasta


#Create a miRNA precursor file to check the 10 to one proportion of mapped fragments - typical characterist com true miRNAs

cat MiradorPrecursorsFixIDs.fa CArabicaPrecursorMir.fasta > t

mv t CArabicaPrecursorMir.fasta

for file in `ls ../../*fasta`;
do

#all miRNA
#shortStack v 3.8.5
~/bin/ShortStack/ShortStack --g CArabicaPrecursorMir.fasta --bowtie_cores 55 --readfile $file --mmap u --nohp --mismatches 0 --outdir `basename -s fasta $file`noMismatch.ShortStack

done


###
#G1

#ONLY miRNA PRECURSOR

mkdir processBam
cd processBam
for folder in `ls -d ../*noMismatch.ShortStack/`
do
    echo $folder
    cd $folder
    pwd

    samtools sort -@ 55 $(basename -s noMismatch.ShortStack $folder)bam > $(basename -s noMismatch.ShortStack $folder)sorted.bam

#keep only mapped

    samtools view -b --threads 55 -F 4 $(basename -s noMismatch.ShortStack $folder)sorted.bam > $(basename -s noMismatch.ShortStack $folder)sorted.Mapped.bam

    rm $(basename -s noMismatch.ShortStack $folder)sorted.bam

    cp $(basename -s noMismatch.ShortStack $folder)sorted.Mapped.bam ../processBam


done

cd ../processBam

samtools merge -@ 55 miRprecursor.merged.bam *sorted.Mapped.bam


samtools index -@ 55 miRprecursor.merged.bam
######################################
#some overal exploratory R analysis with all the data from all the libs not allowing mismatches.

data <- read.delim("shortSatckAllResults.txt")



data$type <- NA

#create the type entries

data$type[grep(paste(c("miR","Novel"),collapse="|"), data$Locus)] <- "miR"

data$type[grep(paste(c("PHAS"),collapse="|"), data$Locus)] <- "PHAS"

data$type[grep(paste(c("SNORNA"),collapse="|"), data$Locus)] <- "SNORNA"

data$type[grep(paste(c("SNRNA"),collapse="|"), data$Locus)] <- "SNRNA"

data$type[grep(paste(c("TRNA"),collapse="|"), data$Locus)] <- "TRNA"


data$type <- as.factor(data$type)

data$Condition <- NA

data$Condition[grep(paste(c("Flower"),collapse="|"), data$lib)] <- "Flower"

data$Condition[grep(paste(c("G1"),collapse="|"), data$lib)] <- "Node"

data$Condition[grep(paste(c("G2"),collapse="|"), data$lib)] <- "G1"

data$Condition[grep(paste(c("G3"),collapse="|"), data$lib)] <- "G2"

data$Condition[grep(paste(c("G4_Cat_bud36e","G4_Sir_bud36e"),collapse="|"), data$lib)] <- "G3"

data$Condition[grep(paste(c("G4_Cat_bud36l","G4_Sir_bud36l"),collapse="|"), data$lib)] <- "G4"


data$Condition[grep(paste(c("G5"),collapse="|"), data$lib)] <- "G5"

data$Condition[grep(paste(c("G6"),collapse="|"), data$lib)] <- "G6"

data$Condition <- as.factor(data$Condition)

#major mapped ~ type
boxplot(Complexity ~ type, data =data)

boxplot(FracTop ~ type, data =data)

boxplot(PhaseScore ~ type, data =data)


 dicerCallData<- data[data$DicerCall != "N"  & !is.na(data$DicerCall),c(1,2,13,23,24)]



#barplot(as.numeric(as.character(dicerCallData$DicerCall)) ~ dicerCallData$type)

#now lets get only the PHAS data
PHASdicerCallData <-dicerCallData[dicerCallData$type == "PHAS",]

PHASdicerCallData <- PHASdicerCallData[!is.na(PHASdicerCallData$Locus),]


phasingMatrix <- matrix(nrow = 0,ncol = 2)
for( phas in unique(PHASdicerCallData$Locus)){

selected <-  PHASdicerCallData[PHASdicerCallData$Locus == paste(phas),]

print(phas)
print(selected)
print(median(as.numeric(as.character(selected$DicerCall))))

if(nrow(phasingMatrix) == 0) {

    phasingMatrix <- matrix(c(paste(phas), median(as.numeric(as.character(selected$DicerCall)))),ncol=2)
}else{

    phasingMatrix <- rbind(phasingMatrix, matrix(c(paste(phas), median(as.numeric(as.character(selected$DicerCall)))),ncol=2))

    }
}

phasingMatrix <- as.data.frame(phasingMatrix)

colnames(phasingMatrix) <- c("Locus","phasing")

write.table(phasingMatrix,sep="\t",file="phasing.tab")

#save.image("shortStackStats.RData")

##############################################
######Create a table only with the 24mers + strand####
################################################

Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

phasingPlusMatrix <- matrix(nrow = 0,ncol = 3)

for( phas in unique(PHASdicerCallData$Locus)){

selected <-  PHASdicerCallData[PHASdicerCallData$Locus == paste(phas),]

print(phas)
#print(selected)
print(median(as.numeric(as.character(selected$DicerCall))))
print(as.character(selected$Strand))
print(Modes(as.character(selected$Strand)))

if(nrow(phasingPlusMatrix) == 0) {

    if(any(as.character(selected$Strand) == "+") & median(as.numeric(as.character(selected$DicerCall))) == 24 ){
        phasingPlusMatrix <- matrix(c(paste(phas), median(as.numeric(as.character(selected$DicerCall))),"+"),ncol=3)

        }else{print("Jaquaratutibaia")}

    }else{

    if(any(as.character(selected$Strand) == "+") & median(as.numeric(as.character(selected$DicerCall))) == 24 ){
        phasingPlusMatrix <- rbind(phasingPlusMatrix, matrix(c(paste(phas), median(as.numeric(as.character(selected$DicerCall))),"+"),ncol=3))

        }else{print("Jaquaratutibaia")}
    }
}

phasingPlusMatrix <- as.data.frame(phasingPlusMatrix)

colnames(phasingPlusMatrix) <- c("Locus","phasing","Strand")

write.table(phasingPlusMatrix,sep="\t",file="phasingPlusStrand.tab")

save.image("shortStackStats.RData")

#########################
mkdir seachTargets.miRNA-PHAS

cd seachTargets.miRNA-PHAS

cp ../DE/smallReference/collapesed.C.a.Mature.Mir.fasta ./

#alging all reads to the collapsed tRNA fasta file
bowtie-build collapesed.C.a.Mature.Mir.fasta collapesed.C.a.Mature.Mir.fasta

cat ../DE/*fasta > all.reads.fasta

bowtie -x collapesed.C.a.Mature.Mir.fasta -f all.reads.fasta -y -k 1 -v 0 -p 60 --no-unal -S > bowtie.result.miRNA.sam


# reads processed: 1233211059
# reads with at least one alignment: 199,343,851 (16.16%)
# reads that failed to align: 1033867208 (83.84%)
#Reported 199,343,851 alignments


/home/tcherubino/bin/FastQC/fastqc bowtie.result.miRNA.sam


samtools sort -@ 60 bowtie.result.miRNA.sam > bowtie.result.miRNA.firstRound.sorted.bam

cat bowtie.result.miRNA.firstRound.sorted.bam | samtools view -h -F 4 | samtools idxstats -  > each.mature.tab

#retrieve the reads mapped only to tRNAs
samtools fasta -@ 60 bowtie.result.miRNA.firstRound.sorted.bam > mappedTomiRPrecursors.fasta

sed '/^>/d' mappedTomiRPrecursors.fasta | sort | uniq -c > sequence_count


awk '{print $2"\t"$1}' sequence_count > mappedTomiRNAs.tag


####################################
R

data <- read.delim("mappedTomiRNAs.tag",header=F)

#Count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor. Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)

scallinFactor <- sum(data$V2)/1000000

data$V3 <- data$V2/scallinFactor

#filter sequences larger than 24nt
SizeTreshhold <- data[nchar(data$V1)<= 24,]


#lets filter reads with CPM less than 1

CPM.1.treshold <- SizeTreshhold[SizeTreshhold$V3>1,]

CPM.1.treshold$V4 <- paste0(nchar(CPM.1.treshold$V1),"_",CPM.1.treshold$V1,"_",CPM.1.treshold$V2)

CPM.1 <- cbind(CPM.1.treshold$V4,CPM.1.treshold$V1)

write.table(CPM.1, file="processedReadsMappedTomiTNAs.tab",row.names=F,col.names=F,sep="\t",quote=F)


awk -F'\t' '{printf ">%s\n%s\n",$1,$2}' processedReadsMappedTomiTNAs.tab > mappedTomiRPrecursors.fasta


####################################


bowtie -x ~/Carabica.smallRNA/processedGenome/ToBlakeSystem/whole.C.a.formated.genome -f ./mappedTomiRPrecursors.fasta -y -a -v 3 -p 60 -S > bowtie.result.miR.pre.sam

# reads processed: 353
# reads with at least one alignment: 353 (100.00%)
# reads that failed to align: 0 (0.00%)
#Reported 97383 alignments


samtools sort -@ 60 bowtie.result.miR.pre.sam > bowtie.result.miR.3miss.sorted.bam


samtools index -@ 60 bowtie.result.miR.3miss.sorted.bam
rm  all.reads.fasta

samtools coverage bowtie.result.miR.3miss.sorted.bam -m -o bowtie.result.miR.3miss.sorted.coverage

#############################################
#allow 1 mismatches
############################################
bowtie -x collapesed.C.a.Mature.Mir.fasta -f all.reads.fasta -y -k 1 -v 1 -p 60 --no-unal -S > bowtie.result.miRNA.miss1.sam



# reads processed: 1233211059
# reads with at least one alignment: 206965575 (16.78%)
# reads that failed to align: 1026245484 (83.22%)
#Reported 206965575 alignments


/home/tcherubino/bin/FastQC/fastqc bowtie.result.miRNA.miss1.sam


samtools sort -@ 60 bowtie.result.miRNA.miss1.sam > bowtie.result.miRNA.miss1.firstRound.sorted.bam

cat bowtie.result.miRNA.miss1.firstRound.sorted.bam | samtools view -h -F 4 | samtools idxstats -  |  cut -f 1,3 > miss.1.each.mature.tab

#allow 2 mismatches

bowtie -x collapesed.C.a.Mature.Mir.fasta -f all.reads.fasta -y -k 1 -v 2 -p 60 --no-unal -S > bowtie.result.miRNA.miss2.sam


# reads processed: 1233211059
# reads with at least one alignment: 208541090 (16.91%)
# reads that failed to align: 1024669969 (83.09%)
#Reported 208,541,090 alignments


/home/tcherubino/bin/FastQC/fastqc bowtie.result.miRNA.miss2.sam


samtools sort -@ 60 bowtie.result.miRNA.miss2.sam > bowtie.result.miRNA.miss2.firstRound.sorted.bam

cat bowtie.result.miRNA.miss1.firstRound.sorted.bam | samtools view -h -F 4 | samtools idxstats -  |  cut -f 1,3 > miss.2.each.mature.tab
