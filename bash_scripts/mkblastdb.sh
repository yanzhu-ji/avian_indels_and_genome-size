#!/bin/bash
#PBS -l walltime=04:00:00

module load blast

cd $PBS_O_WORKDIR

species=(calAnn balReg tytAlb halLeu melUnd corBra falPer carCri colLiv melGal anaPla)
species34=(   acaChl apaVit aptFor bucRhi capCar 
              catAur chaPel chaVoc chlMac colStr 
              cucCan egrGar eurHel fulGla gavSte 
              geoFor halAlb lepDis manVit merNub 
              mesUni nesNot nipNip ophHoa pelCri 
              phaCar phaLep phoRub picPub podCri 
              pteGut pygAde tauEry tinGut)

for sp in "${species34[@]}"; do
	makeblastdb -in $RCAC_SCRATCH/CR1-blast/blast-db/${sp}_genomic.fna -dbtype nucl
	echo " $sp done!"
done
