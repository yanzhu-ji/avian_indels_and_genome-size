#!/bin/bash
### Description: this script is used to analyze orthologous CR1 among 48 avian species
### Yanzhu Ji yanzhuji20@gmail.com 
### last modified: 09/2016

### 0.1 download genomes from NCBI: 14 downloaded manually by, e.g.,
#     wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000699145.1_ASM69914v1/GCA_000699145.1_ASM69914v1_genomic.fna.gz
#   Others are downloaded automatically.  
#     For scripts, see ../blast-db/wget.sh and ../blast-db/automated_wget.sh in 
#   0.2 Then, make blastdb for each genome assembly.  
#     For scripts, see ../blast-db/mkblastdb.sh 


### 1.1 blast CR1-Y2_Aves.fasta against each species (or selected ones)

species48=(acaChl anaPla apaVit aptFor balReg
               bucRhi calAnn capCar carCri catAur
               chaPel chaVoc chlMac colLiv colStr
               corBra cucCan egrGar eurHel falPer
               fulGla galGal gavSte geoFor halAlb
               halLeu lepDis manVit melGal melUnd
               merNub mesUni nesNot nipNip ophHoa
               pelCri phaCar phaLep phoRub picPub
               podCri pteGut pygAde strCam taeGut
               tauEry tinGut tytAlb)

cd bed-and-fmt6

for s in "${species48[@]}"; do
  echo '#!/bin/bash
  #PBS -l walltime=04:00:00
  #PBS -l nodes=1:ppn=48:b

  module load blast

  cd $PBS_O_WORKDIR' > CR1-Y2-Aves_${s}_blast.sh

  echo "blastn -task blastn \
   -db $RCAC_SCRATCH/CR1-blast/blast-db/${s}_genomic.fna \
   -query $RCAC_SCRATCH/CR1-blast/blast_CR1-Y2_Aves_48-species/CR1-Y2_Aves.fasta \
   -out CR1-Y2-Aves_${s}.fmt6 \
   -outfmt 6 \
   -num_threads 48 \
   -evalue 0.00001" >> CR1-Y2-Aves_${s}_blast.sh
done

for s in "${species48[@]}"; do
  qsub CR1-Y2-Aves_${s}_blast.sh
done


### 1.2 select hits that matched to the 3' end
for s in "${species48[@]}"; do
  awk '$8 > 3300' CR1-Y2-Aves_${s}.fmt6 > CR1-Y2-Aves_${s}_3300.fmt6
  perl ~/perl/format-transfer.pl -i fmt6 CR1-Y2-Aves_${s}_3300.fmt6
  perl ~/perl/clip_bed.pl -i CR1-Y2-Aves_${s}_3300.bed \
    -s $RCAC_SCRATCH/CR1-blast/blast-db/${s}_genomic.fna \
    -b 500 -a 500 -e full > CR1-Y2-Aves_${s}_3300_f500.fasta
done

# stats of hits across species
for s in "${species48[@]}"; do
  allHits=`wc -l CR1-Y2-Aves_${s}.fmt6 | cut -f1 -d' '`
  hits3300=`wc -l CR1-Y2-Aves_${s}_3300.fmt6 | cut -f1 -d' '`
  clip3300=`grep -c ">" CR1-Y2-Aves_${s}_3300_f500.fasta`
  echo  -e "$s\t$allHits\t$hits3300\t$clip3300"
done > CR1-Y2_Aves_blast-hits_stats.txt

cd ..

### 1.3 make blast database for CR1 part of each species 
for s in "${species48[@]}"; do
  mkdir $s
  cd $s
  mv ../CR1-Y2-Aves_${s}_3300_f500.fasta .

  makeblastdb -in CR1-Y2-Aves_${s}_3300_f500.fasta -dbtype nucl

  cd ..
done

# manually moved *.fmt6, *.sh (for each species), and *.bed to "./bed-and-fmt6".


### CR1 in each species except for the query species (downy woodpecker, picPub) are blasted against using the query species.
querySpecies=(picPub)

# first batch and second batch to make this faster
species_1=(acaChl anaPla apaVit aptFor balReg  
               bucRhi calAnn capCar carCri catAur
               chaPel chaVoc chlMac colLiv colStr
               corBra cucCan egrGar eurHel falPer
               fulGla galGal gavSte geoFor halAlb)

species_2=(halLeu lepDis manVit melGal melUnd
               merNub mesUni nesNot nipNip 
               pelCri phaCar phaLep phoRub picPub
               podCri pteGut pygAde strCam taeGut
               tauEry tinGut tytAlb)

for q in "${querySpecies[@]}"; do
  for s in "${species_2[@]}"; do  # "${species_1[@]}" for the first time running
    if [ $q != $s ]; then
      cd $q
      # echo "Do something..."
      bash ../CR1-against-all.sh CR1-Y2-Aves_${q}_3300_f500.fasta ../${s}/CR1-Y2-Aves_${s}_3300_f500.fasta
      cd ..
    else
      echo "Query species equals to database species...skipping $s."
    fi
  done
done

cd ${querySpecies} # cd picPub

species47=(acaChl anaPla apaVit aptFor balReg  
               bucRhi calAnn capCar carCri catAur
               chaPel chaVoc chlMac colLiv colStr
               corBra cucCan egrGar eurHel falPer
               fulGla galGal gavSte geoFor halAlb
               halLeu lepDis manVit melGal melUnd
               merNub mesUni nesNot nipNip ophHoa
               pelCri phaCar phaLep phoRub 
               podCri pteGut pygAde strCam taeGut
               tauEry tinGut tytAlb)

# concatenate two fmt6 files from splitted fasta files

for s in "${species47[@]}"; do
  cat CR1-Y2-Aves_${querySpecies}_3300_f500_0_${s}_400.fmt6 CR1-Y2-Aves_${querySpecies}_3300_f500_1_${s}_400.fmt6 \
      > CR1-Y2-Aves_${querySpecies}_3300_f500_${s}_400.fmt6 
done

### for each fmt6 from the same species, select unique hits in the database species

# check if file exists; if so then remove it as later ">>" is used (which is dangerous!)
if [ -f ${querySpecies}.uniq.list ]; then
    rm ${querySpecies}.uniq.list
    echo ${querySpecies}.uniq.list removed
fi

for i in `ls CR1-Y2-Aves_${querySpecies}_3300_f500_??????_400.fmt6`; do 
  base=${i%%\.fmt6}
  cut -f1 $i | uniq -u >> ${querySpecies}.uniq.list
  cut -f1 $i | uniq -u  | grep -F -f - $i > ${base}_uniq.fmt6
  sort ${querySpecies}.uniq.list | uniq -c > ${querySpecies}_uniq-count.txt
done

# from ${querySpecies}_uniq-count.txt, sort count from large to small, then select only loci with
#      at least 17 number of species out of 47 have a blast hit of this locus.
#      then save this list to test_loci.list
for locus in `cat test_loci.list`; do
  for s in "${species47[@]}"; do
    length=`grep $locus *_0_${s}_400_uniq.fmt6 | awk '{print $4}' `
    echo $s $length
  done
done 

# will have to clip sequences and align again, as the lengths of matches are not the same across species...
#      steps: fmt6 -> bed (w/ strand info)
#             for each locus, get the fasta seqs from all possible species using clip_bed_strand;
#             align each cluster w/ CR1-Y2_Aves.fasta (Muscle)
#             see previous files for more info.


for s in "${species47[@]}"; do
  if grep --quiet -F -f test_loci.list CR1-Y2-Aves_${querySpecies}_3300_f500_${s}_400_uniq.fmt6; then
    grep -F -f test_loci.list CR1-Y2-Aves_${querySpecies}_3300_f500_${s}_400_uniq.fmt6 > all-loci_${s}.fmt6
    perl ~/perl/format-transfer.pl -inputformat fmt6 -outputformat bed all-loci_${s}.fmt6
    perl 
  fi
done

mkdir fasta-3 
# ./fasta-3 is the final directory.  ./fasta and ./fasta-2 are probably experimental directories.

for s in "${species47[@]}"; do
  loci=`cut -f4 all-loci_${s}.bed`
  for locus in "${loci[@]}"; do
    perl ~/perl/clip_bed_strand.pl -i all-loci_${s}.bed -s ../${s}/CR1-Y2-Aves_${s}_3300_f500.fasta \
         -b 100 -a 100 -d >> ${locus}-${s}.fasta
  done
done

for s in "${species47[@]}"; do
  for locus in `cut -f4 all-loci_${s}.bed`; do
      grep $locus all-loci_${s}.bed > ${locus}_${s}_temp.bed
      perl ~/perl/clip_bed_strand.pl -i ${locus}_${s}_temp.bed -s ../${s}/CR1-Y2-Aves_${s}_3300_f500.fasta \
           -b 100 -a 100 -d >> ./fasta/${locus}_${querySpecies}.fasta
      rm ${locus}_${s}_temp.bed
  done
  echo $s done!
done

# cat in sequence from query species (picPub) and muscle align them.
# note: this step is strictly dependent on the format of fasta being catted from, i.e., sequence part makes up 1 line but not multiple lines.

for i in *_${querySpecies}.fasta; do 
  locus=${i%_*}
  grep -A 1 $locus ../CR1-Y2-Aves_${querySpecies}_3300_f500.fasta | cat - ${i} > ${locus}_${querySpecies}_w-${querySpecies}.fasta
  muscle -in ${locus}_${querySpecies}_w-CR1.fasta -out ${locus}_${querySpecies}_w-${querySpecies}_muscle.fasta
done
 
## manually compile a list of consensus sequences resulting from MUSCLE alignments in Jalview...
#    substituted "+" to "N" (plus-2-N).  Takes a while (1.5 hours)
#    headerline of consensus sequences must be formatted as ${locus}_Jalview

# blast-2-seq
blastn -task blastn \
  -query jalview_cons_plus-2-N.fasta \
  -subject ../../CR1-Y2_Aves.fasta \
  -evalue 0.00001 \
  -outfmt 6 \
  -out consensus_CR1-Y2-Aves_E5.fmt6 \

# for loci that have two hits to CR1, temporarily move them to another file and saved the uniq ones to consensus_CR1-Y2-Aves_E5_uniq.fmt6

## clip MSAs and label them w/ species abbreviation
cat consensus_CR1-Y2-Aves_E5_uniq.fmt6 | while read line; do 
  locus=`echo $line | cut -f1 -d" " -`
  locus=${locus%_Jalview}
  cutStart=`echo $line | cut -f7 -d" " -`
  cutEnd=`echo $line | cut -f8 -d" " -`

  ## trim MSAs according to the blast2seq information (start and end in the query/consensus sequences)
  for file in ${locus}_${querySpecies}_w-${querySpecies}_muscle.fasta; do
    grep ">" $file | cut -f1  |sed 's/>//' | sed "s/$/\t$cutStart\t$cutEnd/" > ${locus}_clip_temp.bed
    perl ~/perl/clip_bed.pl -i ${locus}_clip_temp.bed -s $file > ${locus}_trimmed.fasta
    rm ${locus}_clip_temp.bed

    ## for each trimmed fasta file, add species label at the end of sequence ID to facilitate subsequent parsing and stats
    grep ">" ${locus}_trimmed.fasta | while read dbLoci; do
      dbLocus=`echo $dbLoci | cut -f1,2,3 -d"_" - | sed 's/^>//' - `

      if ls ${dbLocus}_trimmed.fasta; then
#t       echo db equals to query
         sed -i "s/${dbLoci}/${dbLoci}_x_y_${querySpecies}/" ${locus}_trimmed.fasta # so that the format is consistent w/ other species and easier to be parsed next step.
      else 
        for s in "${species47[@]}"; do
          if grep --quiet $dbLocus $RCAC_SCRATCH/CR1-blast/blast_CR1-Y2_Aves_48-species/${querySpecies}/all-loci_${s}.bed; then
#t          echo find $s
#t          echo $locus
            sed -i "s/${dbLoci}/${dbLoci}_${s}/" ${locus}_trimmed.fasta 
          fi
        done
      fi
    done
  done
done

for f in *_trimmed.fasta; do
  base=${f%%\.fasta}
  perl ~/perl/getLength.pl $f | cut -f8 -d"_" > ${base}_length.txt
done

rm *.index

# get lengths of consensus sequences...10.21.2015
for i in *_trimmed.fasta; do
  cons -sequence $i -outseq ${i}.cons -name $i
  perl ~/perl/getLength.pl ${i}.cons >> consensus-length.txt
done
# manually corrected "trimmed.fasta" to "trimmed_length.txt" via vi, so that the names match those in dataframe (in R).

### 3.17.2016 started: find gene in and flanking all loci in pairs of species

##  0.  preparation
##  0.1 xenoRefGene tables are downloaded and prepared in bed format (sorted, with #chr); see other directory for commands
##  0.2 files with chromosome sizes are created in other directories.

##  1. examine 8 of 50 loci that are shared among all 5 species with UCSC xenoRefGene (or xeno-p-RefGene)
## first, copy and paste all loci shared among the five species into "all-five-shared_picPub.txt"

sed -i 's/_length\.txt/\.fasta/' all-five-shared_picPub.txt 

sps=(galGal taeGut melGal geoFor melUnd)

for sp in "${sps[@]}"; do
  cat all-five-shared_picPub.txt| while read line; do
    grep $sp ../fasta-3/$line >> ${sp}_7-loci.txt
  done 
done


## 2. for each list of loci in each species, figure out the exact location (without flanking regions)

## 2.0 prepare fasta file with all trimmed loci
cd ../fasta-3

cat KL*_trimmed.fasta > cat_trimmed.fasta  

cd -

## 2.1 species 1
module load blast

for sp in "${sps[@]}"; do
  perl ~/perl/fetchFasta.pl  ../fasta-3/cat_trimmed.fasta ${sp}_7-loci.txt > ${sp}_7-loci-trimmed.fasta
  sed -i 's/-//g' ${sp}_7-loci-trimmed.fasta
done

for sp in "${sps[@]}"; do
  blastn -task blastn \
   -db $RCAC_SCRATCH/CR1-blast/blast-db/${sp}_genomic.fna \
   -query ${sp}_7-loci-trimmed.fasta \
   -out ${sp}_7-loci-trimmed_${sp}.fmt6 \
   -outfmt 6 \
   -evalue 0.00001
done

# check if numbers of loci are match
  
for sp in "${sps[@]}"; do
  awk '$3 == 100 { print }' ${sp}_7-loci-trimmed_${sp}.fmt6 > ${sp}_7-loci-trimmed_${sp}_id100.fmt6
  perl ~/perl/format-transfer.pl -i fmt6 -o bed ${sp}_7-loci-trimmed_${sp}_id100.fmt6
done

# in order that these accession IDs match those in RefGene (or xenoRefGene) downloaded from UCSC, manually changed:
# (1) Genbank ID to chrX in melGal
# (2) removed ".1" of ******.1 in melUnd, geoFor

perl ~/perl/column-transfer.pl -i taeGut_7-loci-trimmed_taeGut_id100.bed,1,1 \
  -r $RCAC_SCRATCH/taeGut/chr2acc_chrAdded.txt,2,1 \
  > taeGut_7-loci-trimmed_taeGut_id100_chr.bed

perl ~/perl/column-transfer.pl -i galGal_7-loci-trimmed_galGal_id100.bed,1,1 \
  -r $RCAC_SCRATCH/galGal/chr2acc_chrAdded.txt,2,1 \
  > galGal_7-loci-trimmed_galGal_id100_chr.bed

# finally, sort bed files

sps=(melGal melUnd geoFor)
for sp in "${sps[@]}"; do
  bedtools sort -i ${sp}_7-loci-trimmed_${sp}_id100.bed > ${sp}_7-loci-trimmed_${sp}_id100_sorted.bed
done

sps=(galGal taeGut)
for sp in "${sps[@]}"; do
  bedtools sort -i ${sp}_7-loci-trimmed_${sp}_id100_chr.bed > ${sp}_7-loci-trimmed_${sp}_id100_sorted.bed
done

## 3. run bedtools or bash script to get overlapping genes (if any) between (1) exact loci and exons (2) extended loci and genes

## 3.1 exact loci with exons----should not have any overlapping (file size = 0)

cd overlap_exon

## first link files for galGal and taeGut
ln -s /scratch/hansen/j/ji20/galGal4/galGal_UCSC_custom_table_xeno-p-RefGene_sorted.bed \
      /scratch/hansen/j/ji20/CR1-blast/blast_CR1-Y2_Aves_48-species/galGal/galGal_UCSC_custom_table_xenoRefGene_sorted.bed

ln -s /scratch/hansen/j/ji20/galGal4/galGal_UCSC_custom_table_xeno-p-RefGene_exon_sorted.bed \
      /scratch/hansen/j/ji20/CR1-blast/blast_CR1-Y2_Aves_48-species/galGal/galGal_UCSC_custom_table_xenoRefGene_exon_sorted.bed

ln -s /scratch/hansen/j/ji20/taeGut/taeGut_UCSC_custom_table_xeno-p-RefGene_sorted.bed \
      /scratch/hansen/j/ji20/CR1-blast/blast_CR1-Y2_Aves_48-species/taeGut/taeGut_UCSC_custom_table_xenoRefGene_sorted.bed

ln -s /scratch/hansen/j/ji20/taeGut/taeGut_UCSC_custom_table_xeno-p-RefGene_exon_sorted.bed \
      /scratch/hansen/j/ji20/CR1-blast/blast_CR1-Y2_Aves_48-species/taeGut/taeGut_UCSC_custom_table_xenoRefGene_exon_sorted.bed

sps=(melGal geoFor melUnd galGal taeGut)

for sp in "${sps[@]}"; do
  bedtools intersect -wa -wb -a ${sp}_7-loci-trimmed_${sp}_id100_sorted.bed  \
   -b $RCAC_SCRATCH/CR1-blast/blast_CR1-Y2_Aves_48-species/${sp}/${sp}_UCSC_custom_table_xenoRefGene_exon_sorted.bed \
   > ${sp}_exon_full.bed
  echo "done $sp"
done

## 3.2 now check if the loci, left-x-loci and right-x-loci have any overlap with genes

# first link chrSizes files for galGal and taeGut

# to be complemented

mkdir overlap-l-r
cd overlap-l-r

for sp in "${sps[@]}"; do
  bash ../gene-loci_intersect.sh $sp 500000 ../${sp}_7-loci-trimmed_${sp}_id100_sorted.bed
done

## 3.3   it is time to summarize the results

# first prepare the bed file with information of machable names of each loci 
## 3.3.0 exons
cd ../overlap_exon

# if there is few genes, this step can be done manually

sed 's/_plus_strand//' ${sp1}_exon_full.bed | sed 's/_minus_strand//' > ${sp1}_exon_full_original-loci.bed
sed 's/_plus_strand//' ${sp2}_exon_full.bed | sed 's/_minus_strand//' > ${sp2}_exon_full_original-loci.bed

perl ~/perl/column-transfer.pl -i ../galGal-taeGut_loci.txt,1,4 \
  -r ${sp}_exon_full_original-loci.bed,4,5 -p \
  > ${sp1}-${sp2}_exon_full_${sp1}-gene.txt

perl ~/perl/column-transfer.pl -i ${sp1}-${sp2}_loci_${sp1}-gene.txt,2,5 \
  -r ${sp2}_exon_full_original-loci.bed,4,5 -p \
  > ${sp1}-${sp2}_exon_full_both-gene.txt

# no overlap with exons in these 7 loci

## 3.3.1 loci itself
cd ..
paste galGal_7-loci.txt \
  melGal_7-loci.txt \
  taeGut_7-loci.txt \
  geoFor_7-loci.txt \
  melUnd_7-loci.txt > 7-loci.txt

sed -i 's/>//g' 7-loci.txt
cd -

perl ~/perl/column-transfer.pl -i ../7-loci.txt,1,6 \
  -r galGal_distinct.bed,4,5 -p \
  > 7-loci_gene1.txt

perl ~/perl/column-transfer.pl -i 7-loci_gene1.txt,2,7 \
  -r melGal_distinct.bed,4,5 -p \
  > 7-loci_gene2.txt

perl ~/perl/column-transfer.pl -i 7-loci_gene2.txt,3,8 \
  -r taeGut_distinct.bed,4,5 -p \
  > 7-loci_gene3.txt

perl ~/perl/column-transfer.pl -i 7-loci_gene3.txt,4,9 \
  -r geoFor_distinct.bed,4,5 -p \
  > 7-loci_gene4.txt

perl ~/perl/column-transfer.pl -i 7-loci_gene4.txt,5,10 \
  -r melUnd_distinct.bed,4,5 -p \
  > 7-loci_gene5.txt

rm 7-loci_gene1.txt
rm 7-loci_gene2.txt
rm 7-loci_gene3.txt
rm 7-loci_gene4.txt

## 3.3.2 x bp right to loci
distance=500000

perl ~/perl/column-transfer.pl -i ../7-loci.txt,1,6 \
  -r galGal_r-${distance}_distinct.bed,4,5 -p \
  > 7-loci_r-${distance}_1.txt

perl ~/perl/column-transfer.pl -i 7-loci_r-${distance}_1.txt,2,7 \
  -r melGal_r-${distance}_distinct.bed,4,5 -p \
  > 7-loci_r-${distance}_2.txt

perl ~/perl/column-transfer.pl -i 7-loci_r-${distance}_2.txt,3,8 \
  -r taeGut_r-${distance}_distinct.bed,4,5 -p \
  > 7-loci_r-${distance}_3.txt

perl ~/perl/column-transfer.pl -i 7-loci_r-${distance}_3.txt,4,9 \
  -r geoFor_r-${distance}_distinct.bed,4,5 -p \
  > 7-loci_r-${distance}_4.txt

perl ~/perl/column-transfer.pl -i 7-loci_r-${distance}_4.txt,5,10 \
  -r melUnd_r-${distance}_distinct.bed,4,5 -p \
  > 7-loci_r-${distance}_5.txt

rm 7-loci_r-${distance}_1.txt
rm 7-loci_r-${distance}_2.txt
rm 7-loci_r-${distance}_3.txt
rm 7-loci_r-${distance}_4.txt



## 3.3.3 x bp left to loci

perl ~/perl/column-transfer.pl -i ../7-loci.txt,1,6 \
  -r galGal_l-${distance}_distinct.bed,4,5 -p \
  > 7-loci_l-${distance}_1.txt

perl ~/perl/column-transfer.pl -i 7-loci_l-${distance}_1.txt,2,7 \
  -r melGal_l-${distance}_distinct.bed,4,5 -p \
  > 7-loci_l-${distance}_2.txt

perl ~/perl/column-transfer.pl -i 7-loci_l-${distance}_2.txt,3,8 \
  -r taeGut_l-${distance}_distinct.bed,4,5 -p \
  > 7-loci_l-${distance}_3.txt

perl ~/perl/column-transfer.pl -i 7-loci_l-${distance}_3.txt,4,9 \
  -r geoFor_l-${distance}_distinct.bed,4,5 -p \
  > 7-loci_l-${distance}_4.txt

perl ~/perl/column-transfer.pl -i 7-loci_l-${distance}_4.txt,5,10 \
  -r melUnd_l-${distance}_distinct.bed,4,5 -p \
  > 7-loci_l-${distance}_5.txt

rm 7-loci_l-${distance}_1.txt
rm 7-loci_l-${distance}_2.txt
rm 7-loci_l-${distance}_3.txt
rm 7-loci_l-${distance}_4.txt

### 4. for each locus that is covered by the five species, see whether they overlap with exons/introns/non-genic regions.
mkdir /scratch/hansen/j/ji20/CR1-blast/blast_CR1-Y2_Aves_48-species/picPub/all-five-loci
cd /scratch/hansen/j/ji20/CR1-blast/blast_CR1-Y2_Aves_48-species/picPub/all-five-loci

sps=(galGal melGal taeGut geoFor melUnd)

for sp in "${sps[@]}"; do
  grep $sp ../fasta-3/KL*.fasta | cut -f2 -d '>' >> ${sp}_all-loci.txt
done

for sp in "${sps[@]}"; do
  perl ~/perl/fetchFasta.pl  ../fasta-3/cat_trimmed.fasta ${sp}_all-loci.txt > ${sp}_all-loci-trimmed.fasta
  sed -i 's/-//g' ${sp}_all-loci-trimmed.fasta
done

for sp in "${sps[@]}"; do
  blastn -task blastn \
   -db $RCAC_SCRATCH/CR1-blast/blast-db/${sp}_genomic.fna \
   -query ${sp}_all-loci-trimmed.fasta \
   -out ${sp}_all-loci-trimmed_${sp}.fmt6 \
   -outfmt 6 \
   -evalue 0.00001
done

# check if numbers of loci are match
  
for sp in "${sps[@]}"; do # && $11 < 1E-50 
  awk '$3 == 100 && $7 != 1 { print }' ${sp}_all-loci-trimmed_${sp}.fmt6 > ${sp}_all-loci-trimmed_${sp}_id100.fmt6
  perl ~/perl/format-transfer.pl -i fmt6 -o bed ${sp}_all-loci-trimmed_${sp}_id100.fmt6
done

# in order that these accession IDs match those in RefGene (or xenoRefGene) downloaded from UCSC, manually changed:
# (1) Genbank ID to chrX in melGal
# (2) removed ".1" of ******.1 in melUnd, geoFor: %s/\.1\t/\t

perl ~/perl/column-transfer.pl -i taeGut_all-loci-trimmed_taeGut_id100.bed,1,1 \
  -r $RCAC_SCRATCH/taeGut/chr2acc_chrAdded.txt,2,1 \
  > taeGut_all-loci-trimmed_taeGut_id100_chr.bed

perl ~/perl/column-transfer.pl -i galGal_all-loci-trimmed_galGal_id100.bed,1,1 \
  -r $RCAC_SCRATCH/galGal/chr2acc_chrAdded.txt,2,1 \
  > galGal_all-loci-trimmed_galGal_id100_chr.bed

# finally, sort bed files

sps=(melGal melUnd geoFor)
for sp in "${sps[@]}"; do
  bedtools sort -i ${sp}_all-loci-trimmed_${sp}_id100.bed > ${sp}_all-loci-trimmed_${sp}_id100_sorted.bed
done

sps=(taeGut galGal)
for sp in "${sps[@]}"; do
  bedtools sort -i ${sp}_all-loci-trimmed_${sp}_id100_chr.bed > ${sp}_all-loci-trimmed_${sp}_id100_sorted.bed
done

# intersect exons

sps=(galGal melGal taeGut geoFor melUnd)
for sp in "${sps[@]}"; do
  bedtools intersect -wa -wb -a ${sp}_all-loci-trimmed_${sp}_id100_sorted.bed  \
   -b $RCAC_SCRATCH/CR1-blast/blast_CR1-Y2_Aves_48-species/${sp}/${sp}_UCSC_custom_table_xenoRefGene_exon_sorted.bed \
   > ${sp}_exon_full.bed
  echo "done $sp"
done

# intersect genes
for sp in "${sps[@]}"; do
  bedtools intersect -wa -wb -a ${sp}_all-loci-trimmed_${sp}_id100_sorted.bed  \
   -b $RCAC_SCRATCH/CR1-blast/blast_CR1-Y2_Aves_48-species/${sp}/${sp}_UCSC_custom_table_xenoRefGene_sorted.bed \
   > ${sp}_all-loci_gene_full.bed
  bedtools merge -c 4,10 -o distinct -i ${sp}_all-loci_gene_full.bed \
  > ${sp}_all-loci_gene_distinct.bed

  echo "done $sp"
done

### gene matrix for loci as long as they are present in one of the five species
mkdir /scratch/hansen/j/ji20/CR1-blast/blast_CR1-Y2_Aves_48-species/picPub/present-five-loci

sps=(galGal melGal taeGut geoFor melUnd)

for sp in "${sps[@]}"; do
   grep ${sp} ../fasta-3/KL*_trimmed.fasta | sed 's/\.\.\/fasta-3\///' | sed 's/:>/\t/' > ${sp}_picPub_ref.txt
done


## get a list of picPub loci, then substitute columns 2-6 with loci in species 1-5

ls ../fasta-3/KL*_trimmed.fasta | sed 's/\.\.\/fasta-3\///' > all-loci_picPub.txt

perl ~/perl/column-transfer.pl -i all-loci_picPub.txt,1,2 \
  -r galGal_picPub_ref.txt,1,2 -p \
  > all-loci_sp1.txt

perl ~/perl/column-transfer.pl -i all-loci_sp1.txt,1,3 \
  -r melGal_picPub_ref.txt,1,2 -p \
  > all-loci_sp2.txt

perl ~/perl/column-transfer.pl -i all-loci_sp2.txt,1,4 \
  -r taeGut_picPub_ref.txt,1,2 -p \
  > all-loci_sp3.txt

perl ~/perl/column-transfer.pl -i all-loci_sp3.txt,1,5 \
  -r geoFor_picPub_ref.txt,1,2 -p \
  > all-loci_sp4.txt

perl ~/perl/column-transfer.pl -i all-loci_sp4.txt,1,6 \
  -r melUnd_picPub_ref.txt,1,2 -p \
  > all-loci_sp5.txt

## similarly, substitute gene names in columns 7-11 based on loci names in columns 2-6.

perl ~/perl/column-transfer.pl -i all-loci_sp5.txt,2,7 \
  -r ../all-five-loci/galGal_all-loci_gene_distinct.bed,4,5 -p \
  > all-gene_sp1.txt

perl ~/perl/column-transfer.pl -i all-gene_sp1.txt,3,8 \
  -r ../all-five-loci/melGal_all-loci_gene_distinct.bed,4,5 -p \
  > all-gene_sp2.txt

perl ~/perl/column-transfer.pl -i all-gene_sp2.txt,4,9 \
  -r ../all-five-loci/taeGut_all-loci_gene_distinct.bed,4,5 -p \
  > all-gene_sp3.txt

perl ~/perl/column-transfer.pl -i all-gene_sp3.txt,5,10 \
  -r ../all-five-loci/geoFor_all-loci_gene_distinct.bed,4,5 -p \
  > all-gene_sp4.txt

perl ~/perl/column-transfer.pl -i all-gene_sp4.txt,6,11 \
  -r ../all-five-loci/melUnd_all-loci_gene_distinct.bed,4,5 -p \
  > all-gene_sp5.txt

# remove intermediate files...

### select fasta files for those with annotation to put into IGV viewer
### for galGal and taeGut, genomes are covered by IGV by default.  Need to load melGal, geoFor and melUnd.
sps=(geoFor melUnd)

for sp in "${sps[@]}"; do
  cut -f1 -d"_" ${sp}_all-loci.txt | sort | uniq > ${sp}_all-loci.prefa
  perl ~/perl/fetchFasta.pl ../../../blast-db/${sp}_genomic.fna ${sp}_all-loci.prefa > ${sp}_all-loci_genomic.fasta
  sed -i 's/\.1//' ${sp}_all-loci_genomic.fasta
done

sp=melGal

cut -f1 -d"_" ${sp}_all-loci.txt | sort | uniq > ${sp}_all-loci.prefa
perl ~/perl/fetchFasta.pl ../../../blast-db/${sp}_genomic.fna ${sp}_all-loci.prefa > ${sp}_all-loci_genomic.fasta

sed -i 's/CM000962\.1/chr1/' melGal_all-loci_genomic.fasta
sed -i 's/CM000963\.1/chr2/' melGal_all-loci_genomic.fasta
sed -i 's/CM000964\.1/chr3/' melGal_all-loci_genomic.fasta
sed -i 's/CM000965\.1/chr4/' melGal_all-loci_genomic.fasta
sed -i 's/CM000969\.1/chr8/' melGal_all-loci_genomic.fasta
