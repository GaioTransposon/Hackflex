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
# PART4: barcode extraction
# PART5: gather useful output

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


## get fragment sizes: 
#for interleaved2 in `ls interleaved2*`
#do
#filename=$(basename $interleaved2)
#N="${filename%.*}"
#cat $interleaved2 | awk '{if(NR%4==2) print length($1)}' > frag_size_$N
#done


# reads BEFORE quality filtering and trimming: library, read avg length, #reads, #bases
rm reads_before_cleaning*
for interleaved2 in `ls interleaved2*`
do
echo $interleaved2 >> reads_before_cleaning_col1 # library
awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' $interleaved2 >> reads_before_cleaning_col2 # read avg length 
grep -E '^[ACTGN]+$' $interleaved2 | wc -l >> reads_before_cleaning_col3 # number of reads
grep -E '^[ACTGN]+$' $interleaved2 | perl -pe 's/[[:space:]]//g' | wc -c >> reads_before_cleaning_col4 # number of bases 
done
paste reads_before_cleaning_col1 reads_before_cleaning_col2 reads_before_cleaning_col3 reads_before_cleaning_col4 > reads_before_cleaning.tsv


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


# reads after quality filtering and trimming: library, read avg length, #reads, #bases
rm reads_after_cleaning*
for trimmed2_library in `ls trimmed2_*`
do
echo $trimmed2_library >> reads_after_cleaning_col1 # library
awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' $trimmed2_library >> reads_after_cleaning_col2 # read avg length 
grep -E '^[ACTGN]+$' $trimmed2_library | wc -l >> reads_after_cleaning_col3 # number of reads
grep -E '^[ACTGN]+$' $trimmed2_library | perl -pe 's/[[:space:]]//g' | wc -c >> reads_after_cleaning_col4 # number of bases 
done
paste reads_after_cleaning_col1 reads_after_cleaning_col2 reads_after_cleaning_col3 reads_after_cleaning_col4 > reads_after_cleaning.tsv


# find smallest library and save the base count and resize libs based on that number 
# int=f because if int(erleaved) is set to true, resizing won't happen as desired. 
# This is probably because of the reads were originally not interleaved. It probably uses header info
# that reformat.sh is not capable of reading correctly. It must be set: int=f. 
export smallest_lib_reads=`cat reads_after_cleaning.tsv | awk '{print $3}' | sort -n | head -n 1`
for trimmed2_library in `ls trimmed2_*`
do
filename_lib=$(basename $trimmed2_library)
lib="${filename_lib%.*}"
reformat.sh int=f samplereadstarget=$smallest_lib_reads in=$trimmed2_library out=reduced\_$lib
done


# get number of reads after library resizing: library, read avg length, #reads, #bases
rm reads_after_resizing*
for reduced in `ls reduced*`
do
echo $reduced >> reads_after_resizing_col1 # library
awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' $reduced >> reads_after_resizing_col2 # read avg length 
grep -E '^[ACTGN]+$' $reduced | wc -l >> reads_after_resizing_col3 # number of reads
grep -E '^[ACTGN]+$' $reduced | perl -pe 's/[[:space:]]//g' | wc -c >> reads_after_resizing_col4 # number of bases 
done
paste reads_after_resizing_col1 reads_after_resizing_col2 reads_after_resizing_col3 reads_after_resizing_col4 > reads_after_resizing.tsv

# concat all stats
paste reads_before_cleaning.tsv reads_after_cleaning.tsv reads_after_resizing.tsv > reads_stats.tsv


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

# extract headers from fastq files
for f in source_data/*.fastq
do filename=$(basename $f)
N="${filename%.*}"
cat $f | grep "@M00" > fastq_files_headers_$N.tsv
done

# add filename as column to each one
for f in fastq_files_headers_*
do filename=$(basename $f)
N="${filename%.*}"
awk '{print FILENAME (NF?" ":"") $0}' $f > fastq2_files_headers_with_file_name_$N.tsv
done

# concatenate all
cat fastq2_files_headers_with_file_name* > all_fastq_headers.tsv 


################################ PART 5 ################################


# save all the useful output into a new directory
mkdir out
cp frag_size* out/.
cp flagstat* out/.
cp GC_* out/.
cp RL_* out/.
cp ME_* out/.
cp reduced*.dedup.tsv out/.
cp reads_stats.tsv out/.
cp all_fastq_headers.tsv out/. # barcodes 


########################################################################



###

#export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_ecoli
#export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_ecoli/source_data/assembly.fasta
#qsub -V library_processing.sh 

#export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_paeruginosa
#export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_paeruginosa/source_data/all_p_aeruginosa.contigs.fasta
#qsub -V library_processing.sh 

#export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_saureus
#export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_saureus/source_data/Saureus.fa
#qsub -V library_processing.sh 

###

#export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/ecoli
#export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/ecoli/source_data/assembly.fasta
#qsub -V library_processing.sh 

#export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/paeruginosa
#export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/paeruginosa/source_data/all_p_aeruginosa.contigs.fasta
#qsub -V library_processing.sh

#export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/saureus
#export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/saureus/source_data/Saureus.fa
#qsub -V library_processing.sh

###

#export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_barcode
#export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_barcode/source_data/assembly.fasta
#qsub -V library_processing.sh 

###