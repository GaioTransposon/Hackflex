##########################
## execution example
## export mydir=/shared/homes/12705859/hackflex_libs/ecoli/main
## qsub -V 001_library_processing.sh
##########################

#!/bin/bash
#PBS -l ncpus=10
#PBS -l walltime=10:00:00
#PBS -l mem=10g
#PBS -N 001_library_processing
#PBS -M daniela.gaio@student.uts.edu.au

source activate py_3.5

cd $mydir

# create forward_reverse.txt file based on files in directory: 
ls *R1.fastq > forward
ls *R2.fastq > reverse
paste forward reverse > forward_reverse.txt

# interleave the pair-end libraries: 
while read  first  second
do
filename_R1=$(basename $first)
R1="${filename_R1%.*}"
reformat.sh in1=$first in2=$second out=interleaved2_$R1.fastq
done < "forward_reverse.txt"


# get number of reads before quality filtering and trimming: 
rm reads_before_clean*
for library in `ls interleaved*.fastq*`
do
echo "$library" >> reads_before_cleaning_col1
cat $library | wc -l >> reads_before_cleaning_col2
done
paste reads_before_cleaning_col* | column -s $'\t' -t > reads_before_cleaning
# divide the line count by 4 --> number of reads 
cat reads_before_cleaning | awk '{print $1,$2/4}' > reads_before_cleaning.tsv


# cleaning step 1: 
for library in `ls interleaved*.fastq`
do
filename_lib=$(basename $library)
lib="${filename_lib%.*}"
bbduk.sh threads=3 int=t in=$library out=trimmed\_$lib ref=/shared/homes/12705859/MG1655/phix.fasta duk=cleaning1\_$lib
done


# cleaning step 2: 
rm cleaning2*
for trimmed_library in `ls trimmed*`
do
filename_lib=$(basename $trimmed_library)
lib="${filename_lib%.*}"
bbduk.sh threads=3 int=t in=$trimmed_library out=trimmed2\_$lib duk=cleaning2\_$lib ref=/shared/homes/12705859/miniconda3/envs/py_3.5/opt/bbmap-38.22-1/resources/adapters.fa maxgc=0.98 trimq=20 qtrim=r entropy=0.5 minavgquality=25
done
cat cleaning2_trim* > cleaning2.txt


# get number of reads after quality filtering and trimming: 
rm reads_after_clean*
for trimmed_library in `ls trimmed2*`
do
echo "$trimmed_library" >> reads_after_cleaning_col1
cat $trimmed_library | wc -l >> reads_after_cleaning_col2
done
paste reads_after_cleaning_col* | column -s $'\t' -t > reads_after_cleaning
# divide the line count by 4 --> number of reads 
cat reads_after_cleaning | awk '{print $1,$2/4}' > reads_after_cleaning.tsv


# open last file created "reads_after_cleaning"; 
# find smallest library (bases=reads*300) $ save the base count
export smallest_lib_bases=`cat reads_after_cleaning.tsv | awk '{print $2*300}' | sort -n | head -n 1`
for trimmed2_library in `ls trimmed2_*`
do
filename_lib=$(basename $trimmed2_library)
lib="${filename_lib%.*}"
reformat.sh int=t -samplebasestarget=$smallest_lib_bases in=$trimmed2_library out=reduced\_$lib
done



