#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR

python Excision.py --input RIL275_RelocaTEi.CombinedGFF.characterized.gff > log 2> log2

echo "Done"
