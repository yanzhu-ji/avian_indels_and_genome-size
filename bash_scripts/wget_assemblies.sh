#!/bin/bash

## for ostrich

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000698965.1_ASM69896v1/GCA_000698965.1_ASM69896v1_genomic.fna.gz
gunzip GCA_000698965.1_ASM69896v1_genomic.fna.gz
mv GCA_000698965.1_ASM69896v1_genomic.fna strCam_genomic.fna

## emperor penguin
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000699145.1_ASM69914v1/GCA_000699145.1_ASM69914v1_genomic.fna.gz

#### 6.15.15 downloading genome assemblies with genome size in genomesize.com
## Anna's hummingbird: Calypte anna
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000699085.1_ASM69908v1/GCA_000699085.1_ASM69908v1_genomic.fna.gz
mv *.gz calAnn_genomic.fna.gz
gunzip calAnn_genomic.fna.gz


## Grey-crowned crane: Balearica regulorum
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000709895.1_ASM70989v1/GCA_000709895.1_ASM70989v1_genomic.fna.gz
mv *.gz balReg_genomic.fna.gz
gunzip balReg_genomic.fna.gz

## barn owl: Tyto alba
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000687205.1_ASM68720v1/GCA_000687205.1_ASM68720v1_genomic.fna.gz
mv *.gz tytAlb_genomics.fna.gz
gunzip tytAlb_genomics.fna.gz

## bald eagle: Haliaeetus leucocephalus
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000737465.1_Haliaeetus_leucocephalus-4.0/GCA_000737465.1_Haliaeetus_leucocephalus-4.0_genomic.fna.gz
mv *.gz halLeu_genomic.fna.gz
gunzip halLeu_genomic.fna.gz

## budgerigar: Melopsittacus undulatus
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000238935.1_Melopsittacus_undulatus_6.3/GCA_000238935.1_Melopsittacus_undulatus_6.3_genomic.fna.gz
mv *.gz melUnd_genomic.fna.gz
gunzip melUnd_genomic.fna.gz

## American crow: Corvus brachyrhynchos
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000691975.1_ASM69197v1/GCA_000691975.1_ASM69197v1_genomic.fna.gz
mv GCA_000691975.1_ASM69197v1_genomic.fna.gz corBra_genomic.fna.gz

## Peregrine falcon: Falco peregrinus
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000337955.1_F_peregrinus_v1.0/GCA_000337955.1_F_peregrinus_v1.0_genomic.fna.gz
mv *.gz falPer_genomic.fna.gz

## red-legged seriema: Cariama cristata
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000690535.1_ASM69053v1/GCA_000690535.1_ASM69053v1_genomic.fna.gz
mv *.gz carCri_genomic.fna.gz
gunzip carCri_genomic.fna.gz


## pigeon: Columba livia
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000337935.1_Cliv_1.0/GCA_000337935.1_Cliv_1.0_genomic.fna.gz
mv *.gz colLiv_genomic.fna.gz
gunzip colLiv_genomic.fna.gz

## turkey: Meleagris gallopavo (although it has been downloaded, but I am not sure if this file is the same...last time downloaded the one with chromosomes)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000146605.2_Turkey_2.01/GCA_000146605.2_Turkey_2.01_genomic.fna.gz
mv GCA_000146605.2_Turkey_2.01_genomic.fna.gz melGal_genomic.fna.gz
gunzip melGal_genomic.fna.gz

## finally, peking duck: Anas platyrhynchos
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000355885.1_BGI_duck_1.0/GCA_000355885.1_BGI_duck_1.0_genomic.fna.gz
mv *.gz anaPla_genomic.fna.gz
gunzip anaPla_genomic.fna.gz

################################################
# 9/18, adding more species

# emperior penguin
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000699145.1_ASM69914v1/GCA_000699145.1_ASM69914v1_genomic.fna.gz

# is there an automatic way?
#     getting ftp site -> download *_genomic.fna.gz -> unzip and rename
module load edirect

nfile=`ls -lt *_genomic.fna.gz | wc -l`

if [ $nfile -gt 0 ]; then
  echo "fna.gz file already existed (n=$nfile)! Please double check."
else
  for line in `cat species-34_Abb-BioProjectID.csv`; do
    abb=`echo $line | cut -f1 -d',' ` 
    bioproject=`echo $line | cut -f2 -d',' `
    ftp=`esearch -db assembly -query  $bioproject | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath | sed 's/\s/\n/g' | grep GCA`
    echo $ftp

    wget -q ${ftp}/*_genomic.fna.gz
  
    mv *_genomic.fna.gz ${abb}_genomic.fna.gz
    gunzip ${abb}_genomic.fna.gz
  done
fi
