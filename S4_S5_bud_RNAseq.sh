#This code started back in 2017!
#screen -S B4xB5rnaseq

screen -r B4xB5rnaseq

mkdir CarabicaB4xB5

cd CarabicaB4xB5

###############
#deconddense##
#############
gzip -d -c $ALLDATA/Expression/Raw/Coffee/RNAseq/SBS/DBI/2017-02-27/G4_C_bud36l_r1.fastq.gz  > G4_C_bud36l_r1.fastq &

gzip -d -c $ALLDATA/Expression/Raw/Coffee/RNAseq/SBS/DBI/2017-02-27/G4_C_bud36l_r2.fastq.gz  > G4_C_bud36l_r2.fastq &

gzip -d -c $ALLDATA/Expression/Raw/Coffee/RNAseq/SBS/DBI/2017-02-27/G4_S_bud36l_r1.fastq.gz > G4_S_bud36l_r1.fastq &

gzip -d -c $ALLDATA/Expression/Raw/Coffee/RNAseq/SBS/DBI/2017-02-27/G4_S_bud36l_r2.fastq.gz > G4_S_bud36l_r2.fastq &

gzip -d -c $ALLDATA/Expression/Raw/Coffee/RNAseq/SBS/DBI/2017-02-27/G5_C_green_r1.fastq.gz > G5_C_green_r1.fastq &

gzip -d -c $ALLDATA/Expression/Raw/Coffee/RNAseq/SBS/DBI/2017-02-27/G5_C_green_r2.fastq.gz >G5_C_green_r2.fastq &

gzip -d -c $ALLDATA/Expression/Raw/Coffee/RNAseq/SBS/DBI/2017-02-27/G5_S_green_r1.fastq.gz >G5_S_green_r1.fastq &

gzip -d -c $ALLDATA/Expression/Raw/Coffee/RNAseq/SBS/DBI/2017-02-27/G5_S_green_r2.fastq.gz > G5_S_green_r2.fastq &

 ~/bin/preprocess.seq/FastQC/fastqc -t 8 *fastq &

######################
###Quality Control####
######################

nano adapters

for lib in $(ls *fastq)
do

 java -jar  ~/bin/preprocess.seq/Trimmomatic/trimmomatic-0.33.jar SE -threads 20 $lib $(basename -s fastq $lib).trimmed.fq ILLUMINACLIP:./adapters:3:25:6 SLIDINGWINDOW:4:28 MINLEN:30;

done

~/bin/preprocess.seq/FastQC/fastqc -t 8 *fq &


rm *fastq
######################################################
################  STAR ALIGNER     ###################
######################################################

#PREPARE REFERENCE

#2.7.10a

/home/tcherubino/bin/STAR/source/STAR runThreadN 20 --runMode genomeGenerate --genomeDir /home/tcherubino/Carabica.smallRNA/processedGenome/ToBlakeSystem --sjdbGTFfile /home/tcherubino/Carabica.smallRNA/processedGenome/ToBlakeSystem/C.a.ProteinCoding.genes.gff --sjdbGTFtagExonParentTranscript Parent --genomeFastaFiles /home/tcherubino/Carabica.smallRNA/processedGenome/ToBlakeSystem/whole.C.a.formated.genome.fa --sjdbGTFfeatureExon CDS


/home/tcherubino/bin/STAR/source/STAR runThreadN 20 --genomeDir /home/tcherubino/Carabica.smallRNA/processedGenome/ToBlakeSystem --readFilesIn G4_C_bud36l_r1..trimmed.fq --outFileNamePrefix G4_C_bud36l_r1

/home/tcherubino/bin/STAR/source/STAR runThreadN 20 --genomeDir /home/tcherubino/Carabica.smallRNA/processedGenome/ToBlakeSystem --readFilesIn G4_C_bud36l_r2..trimmed.fq --outFileNamePrefix G4_C_bud36l_r2


/home/tcherubino/bin/STAR/source/STAR runThreadN 20 --genomeDir /home/tcherubino/Carabica.smallRNA/processedGenome/ToBlakeSystem --readFilesIn G4_S_bud36l_r1..trimmed.fq --outFileNamePrefix G4_S_bud36l_r1

/home/tcherubino/bin/STAR/source/STAR runThreadN 20 --genomeDir /home/tcherubino/Carabica.smallRNA/processedGenome/ToBlakeSystem --readFilesIn G4_S_bud36l_r2..trimmed.fq --outFileNamePrefix G4_S_bud36l_r2

/home/tcherubino/bin/STAR/source/STAR runThreadN 20 --genomeDir /home/tcherubino/Carabica.smallRNA/processedGenome/ToBlakeSystem --readFilesIn G5_C_green_r1..trimmed.fq --outFileNamePrefix G5_C_green_r1

/home/tcherubino/bin/STAR/source/STAR runThreadN 20 --genomeDir /home/tcherubino/Carabica.smallRNA/processedGenome/ToBlakeSystem --readFilesIn G5_C_green_r2..trimmed.fq --outFileNamePrefix G5_C_green_r2

/home/tcherubino/bin/STAR/source/STAR runThreadN 20 --genomeDir /home/tcherubino/Carabica.smallRNA/processedGenome/ToBlakeSystem --readFilesIn G5_S_green_r1..trimmed.fq --outFileNamePrefix G5_S_green_r1

/home/tcherubino/bin/STAR/source/STAR runThreadN 20 --genomeDir /home/tcherubino/Carabica.smallRNA/processedGenome/ToBlakeSystem --readFilesIn G5_S_green_r2..trimmed.fq --outFileNamePrefix G5_S_green_r2


######################################################
################  SAMTOOLS    ########################
######################################################

#Convert sam to bam
for i in `ls *out.sam`;
do
samtools view -@ 20 -F 4 -b -S $i > `basename -s .sam $i`;

done

rm *fq
rm *sam
#sort bam files

for i in `ls *Aligned.out`;
do
samtools sort -@ 20 $i -o `basename $i`.sorted.bam
done



######################################################
################ edgeR rmdup #########################
######################################################

for i in `ls *.bam`;
do
htseq-count -a 10 -t CDS -i Parent -f bam --stranded=no $i /home/tcherubino/Carabica.smallRNA/processedGenome/ToBlakeSystem/C.a.ProteinCoding.genes.gff   > `basename $i`_counts.txt; &
done

###############################################
library("edgeR")
targets <- readTargets()
names <- targets$description
matrix_input <- readDGE(targets, comment.char = "!")

#remove meta Tags
MetaTags <- grep("^__", rownames(matrix_input))
matrix_input <- matrix_input[-MetaTags, ]


reads_before <- sum(matrix_input$counts)
#remove low expressed genes
rnaseqmatrix <- matrix_input$counts
rnaseqmatrix <- rnaseqmatrix[rowMeans(rnaseqmatrix) >=5,]

sum(rnaseqmatrix)/reads_before

conditions = matrix_input$samples[,2]
analysis_matrix <- DGEList(counts = rnaseqmatrix,group = conditions)
colnames(analysis_matrix$counts) <- names
design <- model.matrix(~0+group, data=analysis_matrix$samples)
colnames(design) <- levels(analysis_matrix$samples$group)


#NORMALIZATIONS
analysis_matrix <- calcNormFactors(analysis_matrix)

#To estimate common dispersion:
analysis_matrix <- estimateGLMCommonDisp(analysis_matrix, design)
#To estimate trended dispersions:
analysis_matrix <- estimateGLMTrendedDisp(analysis_matrix, design)
#To estimate tagwise dispersions:
analysis_matrix <- estimateGLMTagwiseDisp(analysis_matrix, design)

#Fit a negative binomial generalized log-linear model to the read counts for each gene.
fit <- glmFit(analysis_matrix,design)

pdf(file = "edgeR_BCV.pdf",h=8,w=8)
plotBCV(analysis_matrix)
dev.off()

system("open edgeR_BCV.pdf")

#An MDS plots shows distances, in terms of biological coeficient of variation (BCV) - An MDS plot shows the relative similarities of the samples.

cols1 <- colorRampPalette(c("#c24a3a","#522c28","#7bc437"))(48)

pdf(file = "edgeR_MDS.pdf",h=15,w=15)
plotMDS(analysis_matrix,col = cols1,cex=1.5)
dev.off()
system("open edgeR_MDS.pdf")


#samples_trees
cpm.matrix <- cpm(analysis_matrix,normalized.lib.sizes=F)
colnames(cpm.matrix) <- names
t.cpm.matrix <- t(cpm.matrix)
sampleTree <- hclust(dist(t.cpm.matrix), method = "average");

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
#abline(h = 15, col = "red")
dev.off()

system("open sampleClustering.pdf")


#analyse densisties before normalization - input
pdf("Densities_input.pdf")
plotDensities( log(matrix_input$counts), legend = "topright")
dev.off()

#analyse densisties low expression filter
pdf("Densities_low_expression_fiter.pdf")
plotDensities( log(rnaseqmatrix), legend = "topright")
dev.off()

#analyse densisties after normalization
pdf("Densities_normalization.pdf")
plotDensities( log(cpm.matrix), legend = "topright")
dev.off()

pdf("Densities_log_cpm_fitted_norm.pdf")
plotDensities(log(cpm(fit$fitted.values, normalized.lib.sizes=TRUE)), legend = "topright")
dev.off()

#histogram of densities log10
pdf("Log10_histogram_normilized.pdf", h=10,w=10)
hist(log(cpm.matrix+1,10), col=gray.colors(19, start = 0.9, end = 0.3))
dev.off()

#histogram of densities Log2
pdf("Log2_histogram_normilized.pdf", h=10,w=10)
hist(log(cpm.matrix+1,2), col=gray.colors(19, start = 0.9, end = 0.3))
dev.off()

#histogram of densities no log
pdf("histogram_normilized.pdf", h=10,w=10)
hist(cpm.matrix, col=gray.colors(19, start = 0.9, end = 0.3))
dev.off()

# Load library for pheatmap
cpm.matrix.corr <- cor(cpm.matrix, method="spearman",use="pairwise.complete.obs")

library("pheatmap")
library(gplots)
pdf(file = "sampleClusteringHeatmap.pdf", width = 50, height = 50);
pheatmap(cpm.matrix.corr, fontsize=50,cellwidth=50, cellheight=50, treeheight_col= 450, treeheight_row=450,angle_col=90, legend=F, cex.legend=1.2)
#heatmap.2(cpm.matrix.corr*10000,trace = "none",margins = c(5, 11),keysize =1 , key.title="",col ="bluered",density.info="none")
dev.off()

system("open sampleClusteringHeatmap.pdf")

system("open *pdf")


#get gene data
genData <- read.delim("C.arabica.geneDescriptions.tab",sep="\t",header=F,fill=T)

colnames(genData) <- c("geneName","transcriptName" ,"description","GO","EC")

#G4 X G5
#7283

lrt_G4_x_G5<- glmLRT(fit, contrast = c(-1,1))
tTags_lrt_G4_x_G5 <- topTags(lrt_G4_x_G5, n= NULL)
keep_sig <- matrix(tTags_lrt_G4_x_G5$table$FDR <= 0.05 & abs(tTags_lrt_G4_x_G5$table$logFC) >= 1)

sig_kept_G4_x_G5<- tTags_lrt_G4_x_G5$table[keep_sig[,1],]

sig_kept_G4_x_G5 <- sig_kept_G4_x_G5[order(sig_kept_G4_x_G5$logFC),]

sig_kept_G4_x_G5 <-cbind(sig_kept_G4_x_G5,genData[match(rownames(sig_kept_G4_x_G5),genData$geneName),])

write.table(sig_kept_G4_x_G5, file = "DE_G4_x_G5.tab", row.names = T, quote = F, sep = "\t")


 #find all mRNAs DE expressed
dedpuplicated_keeped_names <- rownames(sig_kept_G4_x_G5)

save.image("B4_x_B5.RData")
q()

########################
####All go Single line###
########################

 python3 ~/Dropbox/Bioinformatica/pipelines/triks/GOtermUncolpaser.py GO.tab > ~/Dropbox/Bioinformatica/doctors/C.arabica.blake.system/ToBlakeSystem/uncolapsedGO.terms.tab

##########GO#############

awk '{if ($2 != "") print $0 }' UP.G4> UP.G4.GO

python3 ~/Dropbox/Bioinformatica/pipelines/triks/GOtermUncolpaser.py UP.G4.GO > UP.G4.GO.tab

awk '{if ($2 != "") print $0 }' UP.G5> UP.G5.GO

python3 ~/Dropbox/Bioinformatica/pipelines/triks/GOtermUncolpaser.py UP.G5.GO > UP.G5.GO.tab


############################
############OLD#############
############################
######################################################
######################  GO  ##########################
######################################################

#get the IDs for the DE genes
for file in ls *;
	do for wd in `awk '{print $1}' $file` ;
		do grep $wd ../../../proteins/IDs.tab >> `basename $file`.IDs.txt;
	done;
done

#convert a multi-line fasta to single line

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ../../../proteins/coffea_pep.faa > single_line_coffea_pep.faa

#retrieve proteins sequneces


for file in ls *IDs.txt;
	do for wd in `awk '{print $2}' $file` ;
		do grep -A 1 $wd single_line_coffea_pep.faa >> `basename -s .IDs.txt $file`.fasta;
	done;
done
