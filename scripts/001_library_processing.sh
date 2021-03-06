##########################
## execution example
## export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_dilution/ecoli
## export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_dilution/ecoli/source_data/E_coli_MG1655_UTS_nanopore_hybrid_500k.fa
## qsub -V library_processing.sh

##########################

# Script outline: 

# PART1: cleaning
# PART2: mapping
# PART3: coverage (bedgraphs) and alfred 
# PART4: barcode extraction (only for goal_barcode libs)
# PART5: contamination analysis (kraken)
# PART6: gather useful output

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
NEU=${R1::-3}
reformat.sh in1=$first in2=$second out=interleaved_$NEU.fastq
done < "forward_reverse.txt"


# reads BEFORE quality filtering and trimming: library, read avg length, #reads, #bp
rm reads_before_cleaning*
for interleaved in `ls interleaved_*`
do
echo $interleaved >> reads_before_cleaning_col1 # library
awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' $interleaved >> reads_before_cleaning_col2 # read avg length 
cat $interleaved | grep "1:N\|2:N" | wc -l >> reads_before_cleaning_col3 # number of reads
cat $interleaved | paste - - - - | cut -f2 | wc -c >> reads_before_cleaning_col4 # number of bp
done
paste reads_before_cleaning_col1 reads_before_cleaning_col2 reads_before_cleaning_col3 reads_before_cleaning_col4 > reads_before_cleaning.tsv


# cleaning step 1: 
rm cleaning1*
rm trimmed*
for library in `ls interleaved_*.fastq`
do
filename_lib=$(basename $library)
lib="${filename_lib%.*}"
bbduk.sh threads=3 int=t in=$library out=trimmed\_$lib ref=/shared/homes/12705859/MG1655/phix.fasta k=31 hdist=1 duk=cleaning1\_$lib
done


# cleaning step 2: 
rm cleaning2*
rm trimmed2*
for trimmed_library in `ls trimmed*`
do
filename_lib=$(basename $trimmed_library)
lib="${filename_lib%.*}"
bbduk.sh threads=3 int=t in=$trimmed_library out=trimmed2\_$lib adapters=/shared/homes/12705859/HACKLEX_LIBS/adapters/merged2_adapters.fa ktrim=r k=21 mink=11 hdist=1 tpe tbo maxgc=0.98 qtrim=rl qtrim=20 entropy=0.5 minavgquality=0 &> cleaning2\_$lib
done 


# reads after quality filtering and trimming: library, read avg length, #reads, #bp
rm reads_after_cleaning*
for trimmed2_library in `ls trimmed2_*`
do
echo $trimmed2_library >> reads_after_cleaning_col1 # library
awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' $trimmed2_library >> reads_after_cleaning_col2 # read avg length 
cat $trimmed2_library | grep "1:N\|2:N" | wc -l >> reads_after_cleaning_col3 # number of reads
cat $trimmed2_library | paste - - - - | cut -f2 | wc -c >> reads_after_cleaning_col4 # number of bp
done
paste reads_after_cleaning_col1 reads_after_cleaning_col2 reads_after_cleaning_col3 reads_after_cleaning_col4 > reads_after_cleaning.tsv


# find smallest library: 
# based on read counts: use `export smallest_lib_reads=`cat reads_after_cleaning.tsv | awk '{print $3/2}' | sort -n | head -n 1` 
# based on base pairs: export smallest_lib_bp=`cat reads_after_cleaning.tsv | awk '{print $4}' | sort -n | head -n 1` 
# ps: read count must be divided by 2 because samplereadstarget is interpreted as the number of read pairs to sample when you set interleaved=t
rm reduced*
export smallest_lib_bp=`cat reads_after_cleaning.tsv | awk '{print $4}' | sort -n | head -n 1` 
echo $smallest_lib_reads
for trimmed2_library in `ls trimmed2_*`
do
filename_lib=$(basename $trimmed2_library)
lib="${filename_lib%.*}"
reformat.sh int=t samplebasestarget=$smallest_lib_bp in=$trimmed2_library out=reduced\_$lib
done


# get number of reads after library resizing: library, read avg length, #reads, #bp
rm reads_after_resizing*
for reduced in `ls reduced*`
do
echo $reduced >> reads_after_resizing_col1 # library
awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' $reduced >> reads_after_resizing_col2 # read avg length 
cat $reduced | grep "1:N\|2:N" | wc -l >> reads_after_resizing_col3 # number of reads
cat $reduced | paste - - - - | cut -f2 | wc -c >> reads_after_resizing_col4 # number of bp
done
paste reads_after_resizing_col1 reads_after_resizing_col2 reads_after_resizing_col3 reads_after_resizing_col4 > reads_after_resizing.tsv


# concat all stats
rm reads_stats.tsv
paste reads_before_cleaning.tsv reads_after_cleaning.tsv reads_after_resizing.tsv > reads_stats.tsv


################################ PART 2 ################################


# index the reference genome
bwa index $ref_genome


# map interleaved clean library to reference genome
for library in `ls reduced*`
do
filename_lib=$(basename $library)
lib="${filename_lib%.*}"
bwa mem -p $ref_genome $library > $lib.sam 
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


################################ PART 3 ################################


# insert size based on .bam
rm picardIS*
for mybam in `ls *.dedup.bam`
do filename=$(basename $mybam)
N="${filename%.*}"
java -jar /shared/homes/12705859/picard/build/libs/picard.jar CollectInsertSizeMetrics \
      I=$mybam \
      O=picardIS_$N.txt \
      H=picardIS_$N.pdf \
      M=0.4
done


# GC bias based on .bam
rm picardGC*
for mybam in `ls *.dedup.bam`
do filename=$(basename $mybam)
N="${filename%.*}"
java -jar /shared/homes/12705859/picard/build/libs/picard.jar CollectGcBiasMetrics \
      I=$mybam \
      O=picardGC_$N.txt \
      CHART=picardGC_$N.pdf \
      S=picardGC_summarymetrics$N.txt \
      R=$ref_genome
done


# Coverage based on .bam
rm picardCO*
for mybam in `ls *.dedup.bam`
do filename=$(basename $mybam)
N="${filename%.*}"
java -jar /shared/homes/12705859/picard/build/libs/picard.jar CollectWgsMetrics \
      I=$mybam \
      O=picardCO_$N.txt \
      READ_LENGTH=300 \
      R=$ref_genome
done


# bedgraphs to report coverage across all regions of chromosome(s)
rm *bedgraph
for mybam in *positionsort
do filename=$(basename $mybam)
N="${filename%.*}"
genomeCoverageBed -d -ibam $mybam > $N.bedgraph
done


# alfred to get stats
for file in `ls *dedup.bam`
do
filename=$(basename $file)
N="${filename%.*}"
echo $file
alfred qc -r $ref_genome -o qc_$N.tsv.gz $file    # replace ref with $ref_genome
alfred qc -r $ref_genome -j qc_$N.json.gz $file    # replace ref with $ref_genome
done

# extract info from alfred's output
for file in `ls qc*tsv.gz`
do
filename=$(basename $file)
N="${filename%.*}"
echo $file
zgrep ^BQ $file | cut -f 2- > BQ_$N.tsv # # Mean base quality (BQ).
zgrep ^RL $file | cut -f 2- > RL_$N.tsv # Read length distribution (RL)
zgrep ^CM $file | cut -f 2- > CM_$N.tsv # Chromosome mapping statistics (CM)
zgrep ^MQ $file | cut -f 2- > MQ_$N.tsv # Mapping quality histogram (MQ)
zgrep ^ME $file | cut -f 2- > ME_$N.tsv # Alignment summary metrics (ME)
done

################################ PART 4 ################################


# do some extra cleaning for the barcodes from the 96 libs (otherwise file is too large to move around)

if [ $mydir == "/shared/homes/12705859/HACKLEX_LIBS/goal_barcode_v0" ] || [ $mydir == "/shared/homes/12705859/HACKLEX_LIBS/goal_barcode_v1" ] then


   # extract headers from fastq files
   for f in source_data/original/*.fastq.gz
   do filename=$(basename $f)
   N="${filename%.*}"
   zcat $f | grep "@M00" > fastq_files_headers_$N.tsv
   done
   
   # add filename as column to each one
   for f in fastq_files_headers_*
   do filename=$(basename $f)
   N="${filename%.*}"
   awk '{print FILENAME (NF?" ":"") $0}' $f > fastq2_files_headers_with_file_name_$N.tsv
   done
   
   # concatenate all
   cat fastq2_files_headers_with_file_name* > all_fastq_headers.tsv 

   awk -F " " '{print $1}' all_fastq_headers.tsv | cut -c20- | rev | cut -c5- | rev > all_fastq_headers_col1.tsv
   awk -F " " '{print $3}' all_fastq_headers.tsv > all_fastq_headers_col2.tsv
   paste all_fastq_headers_col1.tsv all_fastq_headers_col2.tsv > all_fastq_headers_clean.tsv
   
   # Undetermined: 
   zcat source_data/Undet/*.fastq.gz | grep "@M00" > undet_fastq_files_headers.tsv
   awk -F " " '{print $2}' undet_fastq_files_headers.tsv > undet_fastq_files_headers_clean.tsv

   
fi




################################ PART 5 ################################

# extract unmapped reads
rm unmapped*
for library in `ls *dedup.bam`
do
filename_lib=$(basename $library)
lib="${filename_lib%.*}"
samtools view -f 0x4 $library | awk 'BEGIN{OFS=""}{print ">",$1,"\n",$10,"\n+\n",$11;}' > unmapped_$lib.fastq
done

# get contamination report 
rm kraken_report*
for library in `ls unmapped*fastq`
do
filename_lib=$(basename $library)
lib="${filename_lib%.*}"
singularity exec -B /shared/homes/s1/ docker://quay.io/biocontainers/kraken2:2.0.8_beta--pl526h6bb024c_0 kraken2 --db /shared/homes/s1/minikraken2_v1_8GB $library --use-names --report kraken_report_$lib
done

################################ PART 6 ################################

# save all the useful output into a new directory
mkdir out
mv dups_stats.txt out/.
mv cleaning2* out/.
mv flagstat* out/.
mv RL_* out/.
mv CM_* out/.
mv MQ_* out/.
mv BQ* out/.
mv ME_* out/.
mv reads_stats.tsv out/.
mv picard* out/.
mv all_fastq_headers_clean.tsv out/. # barcodes
mv undet_fastq_files_headers_clean.tsv out/. # undet reads
mv qc_*.json.gz out/. # alfred web interactive - input
mv *bedgraph out/. # bedgraphs report coverage across all the regions of the chromosome(s)
mv kraken_report_* out/.

########################################################################


# ###
# 
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/ecoli
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/ecoli/source_data/E_coli_MG1655_UTS_nanopore_hybrid_500k.fa
# qsub -V library_processing.sh
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/paeruginosa
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/paeruginosa/source_data/all_p_aeruginosa.contigs.fasta
# qsub -V library_processing.sh
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/saureus
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_size_selection/saureus/source_data/Saureus.fa
# qsub -V library_processing.sh
# 
# ###
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_ecoli
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_ecoli/source_data/E_coli_MG1655_UTS_nanopore_hybrid_500k.fa
# qsub -V library_processing.sh
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_paeruginosa
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_paeruginosa/source_data/all_p_aeruginosa.contigs.fasta
# qsub -V library_processing.sh
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_saureus
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_saureus/source_data/Saureus.fa
# qsub -V library_processing.sh
# 
# ###
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_barcode_v1
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_barcode_v1/source_data/E_coli_MG1655_UTS_nanopore_hybrid_500k.fa
# qsub -V library_processing.sh
# 
# ###
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_barcode_v0
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_barcode_v0/source_data/E_coli_MG1655_UTS_nanopore_hybrid_500k.fa
# qsub -V library_processing.sh
# 
# ###
# 
# export mydir=/shared/homes/12705859/HACKLEX_LIBS/goal_barcode_v1_top20
# export ref_genome=/shared/homes/12705859/HACKLEX_LIBS/goal_barcode_v0/source_data/E_coli_MG1655_UTS_nanopore_hybrid_500k.fa
# qsub -V library_processing.sh
# 
# ###


#ls -laSh | head -41 | awk '{print $NF}' | xargs -i cp {} ../../goal_barcode_v1_top20/source_data
