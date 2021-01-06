##########################
## execution example
## export mydir=/shared/homes/12705859/hackflex_libs/ecoli/main
## export ref_genome=/shared/homes/12705859/hackflex_libs/ecoli/main/assembly.fasta
## qsub -V 003_library_alfred.sh
##########################

#!/bin/bash
#PBS -l ncpus=10
#PBS -l walltime=10:00:00
#PBS -l mem=10g
#PBS -N 003_library_alfred
#PBS -M daniela.gaio@student.uts.edu.au

source activate py_3.5

cd $mydir


for file in `ls *dedup.bam`
do
filename=$(basename $file)
N="${filename%.*}"
echo $file
alfred qc -r $ref_genome -o qc_$N $file    # replace ref with $ref_genome
done

for file in `ls qc_*`
do
filename=$(basename $file)
N="${filename%.*}"
echo $file
zgrep ^GC $file | cut -f 2- > GC_$N.tsv
zgrep ^RL $file | cut -f 2- > RL_$N.tsv
zgrep ^ME $file | cut -f 2- > ME_$N.tsv
done