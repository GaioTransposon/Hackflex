##########################
## execution example
## export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_dilution/ecoli
## export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_dilution/ecoli/source_data/assembly.fasta
## qsub -V library_processing.sh

##########################

# Script outline: 

# PART1: cleaning
# PART2: mapping
# PART3: alfred 
# PART4: gather useful output

##########################

#!/bin/bash
#PBS -l ncpus=10
#PBS -l walltime=10:00:00
#PBS -l mem=20g
#PBS -N library_processing
#PBS -M daniela.gaio@student.uts.edu.au

source activate py_3.5

cd $mydir

# start fresh: remove all except the source files
ls | grep -v source_data | xargs rm -r


################################ PART 1 ################################


# create forward_reverse.txt file based on files in directory: 
ls source_data/*R1.fastq > forward
ls source_data/*R2.fastq > reverse
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


################################ PART 2 ################################


# index the reference genome
bwa index $ref_genome


# map interleaved clean library to reference genome
for library in `ls reduced*`
do
filename_lib=$(basename $library)
lib="${filename_lib%.*}"
bwa mem -p $ref_genome $library > $lib.sam # S. aureus to be replaced with $ref_genome
done


# create .bam
for library in `ls *.sam`
do
filename_lib=$(basename $library)
lib="${filename_lib%.*}"
samtools view -Sb $library | samtools sort -o $lib.bam
done


# samtools index and sort:
for library in `ls *.bam`
do
filename_lib=$(basename $library)
lib="${filename_lib%.*}"
samtools index $library
samtools sort -n -o $lib.namesort $library
done


# samtools fixmate:
for library in `ls *.namesort`
do
filename_lib=$(basename $library)
lib="${filename_lib%.*}"
samtools fixmate -m $library $lib.fixmate
done


# samtools re-sort:
for library in `ls *.fixmate`
do
filename_lib=$(basename $library)
lib="${filename_lib%.*}"
samtools sort -o $lib.positionsort $library
done


# samtools remove dups:
rm dups_stats.txt
for library in *.positionsort ;
do
filename_lib=$(basename $library)
lib="${filename_lib%.*}"
echo "##############" >> dups_stats.txt
echo $lib >> dups_stats.txt
samtools markdup -s $library $lib.dedup.bam &>> dups_stats.txt # -r won't remove but reports on dups
done


# run flagstat & save output 
rm flagstat_out*
for file in `ls *.dedup.bam`
do
filename=$(basename $file)
N="${filename%.*}"
samtools flagstat $file > flagstat_out_$N.txt
done

# samtools mpileup & generate tab file for ALFRED:
for library in `ls *.dedup.bam`
do
filename_lib=$(basename $library)
lib="${filename_lib%.*}"
samtools mpileup $library | cut -f 2,4 > $lib.tsv
done


# get number of reads: 
rm final_read_counts*
for library in `ls *dedup.tsv*`
do
echo "$library" >> final_read_counts_col1
cat $library | wc -l >> final_read_counts_col2
done
paste final_read_counts_col* | column -s $'\t' -t > final_read_counts
cat final_read_counts | awk '{print $1,$2,$2/4}' > final_read_counts.tsv


# get fragment sizes: 
for file in `ls *.dedup.bam`
do
filename=$(basename $file)
N="${filename%.*}"
samtools view $file | cut -f 9 > frag_sizes_$N.txt
done



################################ PART 3 ################################

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

################################ PART 4 ################################


# save all the useful output into a new directory
mkdir out
cp flagstat* out/.
cp GC_* out/.
cp RL_* out/.
cp frag* out/.
cp reduced*.dedup.tsv out/.





# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_dilution/ecoli
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_dilution/ecoli/source_data/assembly.fasta
# qsub -V library_processing.sh 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_dilution/paeruginosa
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_dilution/paeruginosa/source_data/all_p_aeruginosa.contigs.fasta
# qsub -V library_processing.sh 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_dilution/saureus
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_dilution/saureus/source_data/Saureus.fa
# qsub -V library_processing.sh 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_polymerase/ecoli
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_polymerase/ecoli/source_data/assembly.fasta
# qsub -V library_processing.sh 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_KAPA/ecoli
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_KAPA/ecoli/source_data/assembly.fasta
# qsub -V library_processing.sh 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_hackflex/ecoli
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_hackflex/ecoli/source_data/assembly.fasta
# qsub -V library_processing.sh 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_hackflex/paeruginosa
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_hackflex/paeruginosa/source_data/all_p_aeruginosa.contigs.fasta
# qsub -V library_processing.sh 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_hackflex/saureus
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_hackflex/saureus/source_data/Saureus.fa
# qsub -V library_processing.sh 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/ecoli
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/ecoli/source_data/assembly.fasta
# qsub -V library_processing.sh 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/paeruginosa
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/paeruginosa/source_data/all_p_aeruginosa.contigs.fasta
# qsub -V library_processing.sh 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/saureus
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/saureus/source_data/Saureus.fa
# qsub -V library_processing.sh 
