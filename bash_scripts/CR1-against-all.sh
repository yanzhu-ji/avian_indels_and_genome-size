#!/bin/bash


#query_base=CR1-Y2-Aves_galGal4_3300_f500
#db_base=CR1-Y2-Aves_taeGut_3300_f500

query=$1
db=$2

queryFile=$(basename $query)
queryBase=${queryFile%%.fasta}
queryAbb=`echo $queryBase| cut -f2 -d"_"`

dbFile=$(basename $db)
dbBase=${dbFile%%.fasta}
dbAbb=`echo $dbBase| cut -f2 -d"_"`

#t echo $queryFile $queryBase $queryAbb
#t echo $dbFile $dbBase $dbAbb

# trying to run batch jobs by one script (at best)

cd /scratch/hansen/j/ji20/CR1-blast/blast_CR1-Y2_Aves_48-species/${queryAbb}

# if query fasta has been splitted (assuming the process is run without interruption, so that the 
# presence of xx_0.fasta means the presence of xx_n.fasta), then skip this step.
if [ -f ${queryBase}_0.fasta ]; then 
  echo "Will not run perl this time."
else
  perl ~/perl/split_fasta.pl -i ${queryFile} \
    -n 50000 \
    -o ${queryBase}_
fi

#t pwd

# get the list of files and write jobs in batch
for i in `ls ${queryBase}_*.fasta`; do
  i_base=${i%%.*}
#t echo ${i_base}
  echo '#!/bin/bash
  #PBS -l walltime=04:00:00
  #PBS -l nodes=1:ppn=48:b

  module load blast

  cd $PBS_O_WORKDIR' > ${i_base}_${dbAbb}_blast.sh

  echo "blastn -task blastn \
	-db $RCAC_SCRATCH/CR1-blast/blast_CR1-Y2_Aves_48-species/${dbAbb}/${dbFile} \
	-query $RCAC_SCRATCH/CR1-blast/blast_CR1-Y2_Aves_48-species/${queryAbb}/${i} \
	-out ${i_base}_${dbAbb}.fmt6 \
	-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen length' \
	-num_threads 48 \
	-evalue 0.00001" >> ${i_base}_${dbAbb}_blast.sh 

  echo "awk '\$7 < 400 && \$8 > \$13 - 400' ${i_base}_${dbAbb}.fmt6 > ${i_base}_${dbAbb}_400.fmt6 " >> ${i_base}_${dbAbb}_blast.sh 

#  echo "gzip ${i_base}_${dbAbb}.fmt6" >> ${i_base}_${dbAbb}_blast.sh 

  qsub ${i_base}_${dbAbb}_blast.sh
  echo ${i_base}_${dbAbb}_blast.sh submitted. 
  echo

done

