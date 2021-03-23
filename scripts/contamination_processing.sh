##########################
## execution example
## export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_dilution/ecoli
## qsub -V contamination_processing.sh

##########################

# Script outline: 

# extracts unmapped reads from bam and creates a contamination profile with kraken2

##########################

#!/bin/bash
#PBS -l ncpus=10
#PBS -l walltime=10:00:00
#PBS -l mem=20g
#PBS -N contamination_processing.sh
#PBS -M daniela.gaio@student.uts.edu.au

source activate py_3.5

cd $mydir


################################################################


# extract unmapped reads
rm unmapped*
for library in `ls *dedup.bam`
do
filename_lib=$(basename $library)
lib="${filename_lib%.*}"
samtools view -f 0x4 $library | awk 'BEGIN{OFS=""}{print ">",$1,"\n",$10,"\n+\n",$11;}' > unmapped_$lib.fastq
done
cd goal

# get contamination report 
rm kraken_report*
for library in `ls unmapped*fastq`
do
filename_lib=$(basename $library)
lib="${filename_lib%.*}"
singularity exec -B /shared/homes/s1/ docker://quay.io/biocontainers/kraken2:2.0.8_beta--pl526h6bb024c_0 kraken2 --db /shared/homes/s1/minikraken2_v1_8GB $library --use-names --report kraken_report_$lib
done


# move output to directory
mv kraken_report_* out/.


################################################################

# ###
# use: 
# 
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/ecoli
# qsub -V contamination_processing.sh
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/paeruginosa
# qsub -V contamination_processing.sh
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/saureus
# qsub -V contamination_processing.sh
# 
# ###
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_ecoli
# qsub -V contamination_processing.sh
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_paeruginosa
# qsub -V contamination_processing.sh
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_saureus
# qsub -V contamination_processing.sh
# 
# ###
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_barcode
# qsub -V contamination_processing.sh
# 
# ###
# 






#for f in `ls kraken_report*.dedup`; do echo $f ; tail -n 12 $f ; done

#for f in `ls kraken_report*.dedup`; do echo $f ; cat $f | grep "S" ; done



