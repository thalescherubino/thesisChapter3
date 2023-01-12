#Server Dashi
screen -S C.arabica.miRador


mkdir Carabica.smallRNA

cd Carabica.smallRNA

tar -zxvf  Carabica.genome.tar.gz
 mkdir processedGenome

 cp *fa ./processedGenome

 rm *fa

 cd processedGenome

 #The genomic files will be processed to be renamed
 #according to the chrIDs from file ChromossomeDescriptions.xlsx
 #This was done because the mySQL system from Blake's lab
 #requires number as primary keys

cat *fa > whole.C.a.formated.genome.fa

#compact the separated files
tar -zcvf processed.Carabica.genome.tar.gz canephoraChromossomes.fa eugenioidesChromossomes.fa unplacedContigs.fa
canephoraChromossomes.fa

rm canephoraChromossomes.fa eugenioidesChromossomes.fa unplacedContigs.fa

mkdir smallRNAsLibraries

cd smallRNAsLibraries

cp /alldata/Expression/Processed/Coffee/sRNA/SBS/DBI/2016-03-22/chopped/* ./

cp /alldata/Expression/Processed/Coffee/sRNA/SBS/DBI/2016-03-28/chopped/* ./

cp /alldata/Expression/Processed/Coffee/sRNA/SBS/DBI/2016-04-01/chopped/* ./

cp /alldata/Expression/Processed/Coffee/sRNA/SBS/DBI/2017-02-27/chopped/* ./

cp /alldata/Expression/Processed/Coffee/sRNA/SBS/DBI/2017-02-27/chopped/

cd ../

mkdir reuslt

python3 miRador miRador.ini


################################

cat predicted.C.canephora.gff predicted.C.eugenioides.gff > whole.C.a.gff

#remove coments
sed -i"" '/^#/d' whole.C.a.gff

#translate the chormossome IDs to numbers

sed -i"" 's/NC_039898.1/1/g' whole.C.a.gff

sed -i"" 's/NC_039899.1/2/g' whole.C.a.gff

sed -i"" 's/NC_039900.1/3/g' whole.C.a.gff

sed -i"" 's/NC_039901.1/4/g' whole.C.a.gff

sed -i"" 's/NC_039902.1/5/g' whole.C.a.gff

sed -i"" 's/NC_039903.1/6/g' whole.C.a.gff

sed -i"" 's/NC_039904.1/7/g' whole.C.a.gff

sed -i"" 's/NC_039905.1/8/g' whole.C.a.gff

sed -i"" 's/NC_039906.1/9/g' whole.C.a.gff

sed -i"" 's/NC_039907.1/10/g' whole.C.a.gff

sed -i"" 's/NC_039908.1/11/g' whole.C.a.gff

sed -i"" 's/NC_039909.1/12/g' whole.C.a.gff

sed -i"" 's/NC_039910.1/13/g' whole.C.a.gff

sed -i"" 's/NC_039911.1/14/g' whole.C.a.gff

sed -i"" 's/NC_039912.1/15/g' whole.C.a.gff

sed -i"" 's/NC_039913.1/16/g' whole.C.a.gff

sed -i"" 's/NC_039914.1/17/g' whole.C.a.gff

sed -i"" 's/NC_039915.1/18/g' whole.C.a.gff

sed -i"" 's/NC_039916.1/19/g' whole.C.a.gff

sed -i"" 's/NC_039917.1/20/g' whole.C.a.gff

sed -i"" 's/NC_039918.1/21/g' whole.C.a.gff

sed -i"" 's/NC_039919.1/22/g' whole.C.a.gff

cp /home/tcherubino/Carabica.smallRNA/processedGenome/whole.C.a.gff /home/tcherubino/Carabica.smallRNA/processedGenome/temp.gff


awk '{if( $3  ==  "gene" || $3  ==  "CDS") print $0}' /home/tcherubino/Carabica.smallRNA/processedGenome/temp.gff > teste.gff

awk -F"\t" '{if ($3  ==  "gene")  print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t""ID="$1$9";Name="$1$9;
else
  print $0
}' teste.gff > teste.2.gff

awk -F"\t|;" '{
  if($3  -eq  "CDS")
  print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10
else
  print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9";"$10
}' teste.2.gff > teste.3.gff

sed -i"" 's/"//g' teste.3.gff
sed -i"" 's/gene_id //g' teste.3.gff
sed -i"" 's/;//g' teste.3.gff


awk '{if($3  ==  "CDS")
  print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t""Parent="$1$9
else
  print $0
}' teste.3.gff > teste4.gff

mv teste4.gff final.gff

cp final.gff toMergeWhithAnnot.gff

#sed 's/Name=/;Name=/g' teste2.gff > final.gff

#The version used in SPARTA is siliglly different with "ID=XX;Name=" for $3==gene . In addition the $3==CDS are bugged just showing the incremental gene ID - and not the chromossome data, but that don't seens to cause problens once sPARTA only use the gene entries.


###########################################################
#retrieve the subGenome specific ID to incude in the final gff file, wich will latter help retrivering the annotation of each gene
###########################################################

mkdir process

cd process

 cp ../*Annot.tab  ./

cp ../toMergeWhithAnnot.gff ./

sed -i"" '/CDS/d' toMergeWhithAnnot.gff

while  IFS= read -r line;
do


  chromossome=$(echo $line | awk '{print $1}')
  #echo $chromossome
  geneToo=$(echo $line | awk '{print $9}')
  echo $geneToo >> newGeneIDS.txt
  gene=$(echo $line | awk '{print $9}' | awk -Fg '{print$2}')


  #echo $gene
  if [ $chromossome -eq 1 ] || [ $chromossome -eq 3 ] || [ $chromossome -eq 5 ] || [ $chromossome -eq 7 ] || [ $chromossome -eq 10 ] || [ $chromossome -eq 11 ] || [ $chromossome -eq 13 ] || [ $chromossome -eq 16 ] || [ $chromossome -eq 17 ] || [ $chromossome -eq 20 ] || [ $chromossome -eq 21 ]

  then
    #echo "canephora"
    grep -w SubC.c_g$gene.t1  subCanephoraAnnot.tab | awk -F'\t' '{print $1"\t"$2"\t"$12"\t"$14}' >> geneDescriptions.tab

  elif [ $chromossome -eq 2 ] || [ $chromossome -eq 4 ] || [ $chromossome -eq 6 ] || [ $chromossome -eq 8 ] || [ $chromossome -eq 9 ] || [ $chromossome -eq 12 ] || [ $chromossome -eq 14 ] || [ $chromossome -eq 15 ] || [ $chromossome -eq 18 ] || [ $chromossome -eq 19 ] || [ $chromossome -eq 22 ]
    then
      #echo "eugenioides"
      grep -w SubC.e_g$gene.t1 subEugenoidesAnnot.tab | awk -F'\t' '{print $1"\t"$2"\t"$12"\t"$14}' >> geneDescriptions.tab

  fi

done < toMergeWhithAnnot.gff

sed -i"" 's/Name=//g' newGeneIDS.txt

paste newGeneIDS.txt geneDescriptions.tab > t

mv t C.arabica.geneDescriptions.tab


cd ..

sed -i"" 's/ .*//' whole.C.a.formated.genome.fa


fold whole.C.a.formated.genome.fa > t

mv t whole.C.a.formated.genome.fa

mv final.gff C.a.ProteinCoding.genes.gff

cd ..

tar -zcvf ToBlakeSystem.tar.gz ToBlakeSystem/


###########################################################
#run with only one PARE library

#Coffea Catuai, G4 stage (Buds >3mm and <6mm - late stage) and with Heuristics
cp $ALLDATA/Expression/Processed/Coffee/PARE/SBS/DBI/2017-05-09/chopped/6878_chopped.txt ./

#Coffea Siriema, G4 stage (Buds >3mm and <6mm - late stage)
cp $ALLDATA/Expression/Processed/Coffee/PARE/SBS/DBI/2017-05-09/chopped/6879_chopped.txt ./

#Coffea Catuai, G5 stage (buds from 6>10mm light green color)
cp $ALLDATA/Expression/Processed/Coffee/PARE/SBS/DBI/2017-05-09/chopped/6880_chopped.txt ./

#Coffea Siriema, G5 stage (buds from 6>10mm light green color)
cp $ALLDATA/Expression/Processed/Coffee/PARE/SBS/DBI/2017-05-09/chopped/6881_chopped.txt ./

#Coffea arabica, cultivar Catuai, pre-meiotic anthers
cp $ALLDATA/Expression/Processed/Coffee/PARE/SBS/DBI/2017-05-09/chopped/6882_chopped.txt ./

#Coffea arabica, cultivar Catuai, meiotic anthers
cp $ALLDATA/Expression/Processed/Coffee/PARE/SBS/DBI/2017-05-09/chopped/6883_chopped.txt ./

#Coffea arabica, cultivar Catuai, post-meiotic anthers
cp $ALLDATA/Expression/Processed/Coffee/PARE/SBS/DBI/2017-05-09/chopped/6884_chopped.txt ./

#GENIC
python3 ~/bin/sPARTA/sPARTA.py -genomeFile /home/tcherubino/Carabica.smallRNA/processedGenome/whole.C.a.formated.genome.fa -annoType GFF -annoFile final.gff -genomeFeature 0 -miRNAFile /home/tcherubino/Carabica.smallRNA/sPARTA/collapesed.C.a.Mature.Mir.fasta -libs *chopped.txt -tarPred H -tarScore N --tag2FASTA --map2DD --validate -accel 60 -minTagLen 18

mv output/ genic.output

#INTERGENIC
python3 ~/bin/sPARTA/sPARTA.py -genomeFile /home/tcherubino/Carabica.smallRNA/processedGenome/whole.C.a.formated.genome.fa -annoType GFF -annoFile final.gff -genomeFeature 1 -miRNAFile /home/tcherubino/Carabica.smallRNA/sPARTA/collapesed.C.a.Mature.Mir.fasta -libs *chopped.txt -tarPred H -tarScore N --tag2FASTA --map2DD --validate -accel 60 -minTagLen 18


mv  output/ InterGenic.output/

#Process the data for

#Genic Regions!

####################

#file=G4l_Cat.Sig.Gen.Targets.tab

for file in $(ls *validated*);
do
  echo working on $file

  while  IFS= read -r line;
  do
    #echo $line
    miRNA=$(echo $line | awk -F, '{print $1}')
    #echo $miRNA

    target=$(echo $line | awk -F, '{print $2}')
    #echo $target
    grep -w $target ../../processedGenome/ToBlakeSystem/C.arabica.geneDescriptions.tab >> $file.target.data.tab
  done < $file

done

mkdir process

cp * process

nano header

for file in $(ls *csv.target.data.tab);
do
  cat header $file > t
  mv t $file
done

for file in $(ls *revmapped.csv);
do
  sed 's/,/"\t"/g' $file | sed 's/"//g' > t
  mv t $file
done

for file in $(ls *revmapped.csv);
do
  paste $file.target.data.tab $file > $(basename -s _chopped_validated_revmapped.csv $file).targets.descriptions.tab
done

mkdir queries

cd queries

#"miR060"
query="miR060"

for libraries in $(ls ../*targets.descriptions.tab)
do
  grep -w $query $libraries > t


  cat header t > $query.$(basename $libraries)
done

#miR828a-5p
query="ceu-miR828b-5p"


for libraries in $(ls ../*targets.descriptions.tab)
do
  grep -w $query $libraries > t


  cat header t > $query.$(basename $libraries)
done

#Mirador-Novel-ceu-346958-3p

query="Mirador-Novel-ceu-346958-3p"


for libraries in $(ls ../*targets.descriptions.tab)
do
  grep -w $query $libraries > t


  cat header t > $query.$(basename $libraries)
done

#####################################

grep ">" matheusPred.C.arabica.mature.fasta | sed "s/>//g" > Matheus.IDs.txt

#For veen and other statistics let's remove miRador know miRNAs that are present in the matheus predictions - keeping only the novel miRNAs from MiRador

sed '/Mirador-ccp/d' All.libs.canephora.miRNA.W.mirador.novel.validated.txt	| sed '/Mirador-ceu/d' > temp

mv temp All.libs.canephora.miRNA.W.mirador.novel.validated.txt

#################

sed '/Mirador-ccp/d' All.libs.eugenioides.miRNA.W.mirador.novel.validated.txt | sed '/Mirador-ceu/d' > temp

mv temp All.libs.eugenioides.miRNA.W.mirador.novel.validated.txt

################

sed '/Mirador-ccp/d' All.libs.intergenic.miRNA.W.mirador.novel.validated.txt |
sed '/Mirador-ceu/d' > temp

mv temp All.libs.intergenic.miRNA.W.mirador.novel.validated.txt

#Now lets do the inverse, remove Matheus-predictions and Mirador kown miRNAs - keeping only the novel miRNAs from MiRador

sed '/ceu-miR/d' All.libs.canephora.miRNA.W.mirador.novel.validated.txt | sed '/ccp-miR/d' | sed '/unp-miR/d' > temp

mv temp All.libs.canephora.only.mirador.novel.validated.txt

##########################
sed '/ceu-miR/d' All.libs.eugenioides.miRNA.W.mirador.novel.validated.txt | sed '/ccp-miR/d' | sed '/unp-miR/d' > temp

mv temp All.libs.eugenioides.only.mirador.novel.validated.txt

##########################
sed '/ceu-miR/d' All.libs.intergenic.miRNA.W.mirador.novel.validated.txt | sed '/ccp-miR/d' | sed '/unp-miR/d' > temp

mv temp All.libs.intergenic.only.mirador.novel.validated.txt


########################################
##########Send genome to sequoia########
########################################

scp /home/tcherubino/Carabica.smallRNA/processedGenome/ToBlakeSystem.tar.gz  tcherubino@sequoia.datasci.danforthcenter.org:/home/tcherubino/C.arabica

############################
########on sequoia##########
############################
#screen -S RepeatMask.coffee
screen -r RepeatMask.coffee

 cd C.arabica/

 tar -zxvf ToBlakeSystem.tar.gz


RepeatMasker -species eudicots -nolow -no_is -parallel 40 -gff -dir ./ ToBlakeSystem/whole.C.a.formated.genome.fa

mv whole.C.a.formated.genome.fa.out whole.C.a.formated.genome.fa.RepeatMasker.out

mv whole.C.a.formated.genome.fa.out.gff whole.C.a.formated.genome.fa.RepeatMasker.out.gff

mv whole.C.a.formated.genome.fa.tbl whole.C.a.formated.genome.fa.RepeatMasker.tbl

tar -zcvf  whole.C.a.formated.genome.fa.RepeatMasker.tar.gz whole.C.a.formated.genome.fa.RepeatMasker.out whole.C.a.formated.genome.fa.RepeatMasker.out.gff whole.C.a.formated.genome.fa.RepeatMasker.tbl whole.C.a.formated.genome.fa.masked
###################################################
#Transforme the sPARTA output to a gff to be visualised with IGV

#Genic
ls

cd /Users/johnjoyce/Dropbox/Bioinformatica/Carabica.miRNAs/PAREsPARTA/genic.output

awk -F "," '{gsub("w","+",$17); print $0 }' All.libs.validated.uniq.csv > t

 awk -F "," '{gsub("c","-",$17); print $0}' t > All.libs.validated.uniq.tab


cat All.libs.validated.uniq.tab | tr ',' '\t' > t

cat t | tr ' ' '\t' > All.libs.validated.uniq.tab

awk '{print $16"\t""SPARTA""\t""Clevage_Site""\t"$19"\t"$20"\t"$6"\t"$17"\t"".""\t""interferingRNA="$1";""TARGET="$2";""gene="$1";""score="$6";""Clevage_Site="$18 }' All.libs.validated.uniq.tab  >  genic.targets.gff

######################################################

#Revese map all the miRNA precursors back to the genome to find their coordinates
mkdir revMAPmiRprecursors

cd revMAPmiRprecursors

blat ../processedGenome/ToBlakeSystem/whole.C.a.formated.genome.fa CArabicaPrecursorMir.fasta -t=dna -q=dna precursors.rev.map.tab -minIdentity=100 -minScore=66

#did some manual curation...

#now lets make the precursor fasta file match with the the curation

cd /Users/johnjoyce/Dropbox/Bioinformatica/Carabica.miRNAs/bamFiles/process

nano manualCuration.tab


R


ccpFamilies <- data$FamiIy[data$subGenome == "ccp"]

ceuFamilies <- data$FamiIy[data$subGenome == "ceu"]
