##########################

# Script outline: 

# PART1: de-interleave the resized libs (done with library_processing.sh)
# PART2: run assembly with shovill
# PART3: get assembly stats 

##########################

#!/bin/bash
#PBS -l ncpus=10
#PBS -l walltime=10:00:00
#PBS -l mem=20g
#PBS -N run_shovill
#PBS -M daniela.gaio@student.uts.edu.au


cd $mydir

source activate py_3.5

################################ PART 1 ################################


# deinterleave the libraries: 
for interleaved in `ls red* | grep -v '\.'`
do
filename_lib=$(basename $interleaved)
lib="${filename_lib%.*}"
reformat.sh in=$lib out1=R1_$lib.fq out2=R2_$lib.fq
done 


# create forward_reverse_fq.txt file based on files just created: 
ls R1*.fq > forward_fq
ls R2*.fq > reverse_fq
paste forward_fq reverse_fq > forward_reverse_fq.txt


################################ PART 2 ################################


# assembly by first subsampling to 200 depth
while read  first  second
do
filename_R1=$(basename $first)
R1="${filename_R1%.*}"
singularity exec /shared/homes/12705859/shovill.sif shovill \
      --R1 $first \
      --R2 $second \
      --depth 200 \
      --outdir shovillsub_200_$R1 \
      --force # force to overwrite the output folder 
done < "forward_reverse_fq.txt"


# assembly by first subsampling to 100 depth
while read  first  second
do
filename_R1=$(basename $first)
R1="${filename_R1%.*}"
singularity exec /shared/homes/12705859/shovill.sif shovill \
      --R1 $first \
      --R2 $second \
      --depth 100 \
      --outdir shovillsub_100_$R1 \
      --force # force to overwrite the output folder 
done < "forward_reverse_fq.txt"


# assembly by first subsampling to 50 depth
while read  first  second
do
filename_R1=$(basename $first)
R1="${filename_R1%.*}"
singularity exec /shared/homes/12705859/shovill.sif shovill \
      --R1 $first \
      --R2 $second \
      --depth 50 \
      --outdir shovillsub_50_$R1 \
      --force # force to overwrite the output folder 
done < "forward_reverse_fq.txt"

# assembly by first subsampling to 20 depth
while read  first  second
do
filename_R1=$(basename $first)
R1="${filename_R1%.*}"
singularity exec /shared/homes/12705859/shovill.sif shovill \
      --R1 $first \
      --R2 $second \
      --depth 20 \
      --outdir shovillsub_20_$R1 \
      --force # force to overwrite the output folder 
done < "forward_reverse_fq.txt"

# assembly by first subsampling to 10 depth
while read  first  second
do
filename_R1=$(basename $first)
R1="${filename_R1%.*}"
singularity exec /shared/homes/12705859/shovill.sif shovill \
      --R1 $first \
      --R2 $second \
      --depth 10 \
      --outdir shovillsub_10_$R1 \
      --force # force to overwrite the output folder 
done < "forward_reverse_fq.txt"

# assembly by first subsampling to 5 depth
while read  first  second
do
filename_R1=$(basename $first)
R1="${filename_R1%.*}"
singularity exec /shared/homes/12705859/shovill.sif shovill \
      --R1 $first \
      --R2 $second \
      --depth 5 \
      --outdir shovillsub_5_$R1 \
      --force # force to overwrite the output folder 
done < "forward_reverse_fq.txt"

mv shovill* out/.


################################ PART 3 ################################

# write some stats about the assemblies: 
export these_assemblies=`ls out/shovillsub*/contigs.fa`
statswrapper.sh addname=t in=$these_assemblies &> assembly_stats.txt
mv assembly_stats.txt out/.


########################################################################


# ###
# 
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/ecoli
# qsub -V run_shovill.sh
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/paeruginosa
# qsub -V run_shovill.sh
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/saureus
# qsub -V run_shovill.sh
# 
# ###


########################################################################




