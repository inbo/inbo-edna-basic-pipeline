#!/bin/bash -l

#SBATCH --clusters=wice
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --out=demultiplex-%j.out
#SBATCH --job-name=demultiplex-pipeline

# This script is meant to be submitted to a HPC system with SLURM job sheduler.
# SBATCH arguments might need to be adapted to fit your system/account/task

#### INPUT:
## Job array:
# (-d|--data) INPUT_DATA: path to data.csv with the variables as column names and corresponding values as rows (see atools)
# The required columns are: "BAR","FASTQF","FASTQR"

## Command-line arguments take precedence over parameters specified with (-d) INPUT_DATA and atools!
# (-b|--barcodes) BAR: path to barcode file TSV
#     FORMAT: The TSV should be 3 tab-separated columns with no header: <sample_id> <FWD-index> <REV-index>
#     IMPORTANT: <sample_id> should be unique!
# (-F|--fwd) FASTQF: path to forward reads FASTQ of multiplexed run
# (-R|--rev) FASTQR: path to reverse reads FASTQ of multiplexed run

#### RUN description:
# - Preparation: Create CSV file with input for demultiplexing script => "barcode_file.txt, forward.fastq, reverse.fastq"
# - Launch job, see example below.
# - This script creates a unique subfolder in the execution dir to seperate the workers' working dirs.
#    Which is important because folders and files with hardcoded names are generated in each worker -> interference when in same wd

#### EXAMPLE job submission:
# sbatch --array=1-15%5 .../01_demultiplex.slurm -d data.csv
#### You can also submit the script with command line arguments:
# sbatch 01_demultiplex.slurm -F R1.fastq -R R2.fastq -b barcodes.txt

#### References:
# - cutadapt: Marcel Martin. 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. https://doi.org/10.14806/ej.17.1.200
# - seqkit: Wei Shen*, Botond Sipos, and Liuyang Zhao. 2024. SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing. iMeta e191. doi:10.1002/imt2.191.

################################
# Commandline Argument Parsing #
################################
set -e # exit on error

module purge

# in case of no arguments, show help message
if [ $# == 0 ]; then
    printf "
    
[ERROR] No input found!

Usage: sbatch 01_demultiplex.slurm -F R1.fastq -R R2.fastq -b barcodes.txt
Where:
    -F : input FASTQ file which contains the demultiplexed FORWARD sequences
    -R : input FASTQ file which contains the demultiplexed REVERSED sequences
	  -b : path to barcode file CSV

Alternatively, use an array job to process multiple files in parallel (recommended)
Usage: sbatch --array=1-15%5 .../01_demultiplex.slurm -d data.csv
Where:
    -d : CSV file with columns corresponding to variable names and rows to values (see atools)

"
    exit 1
fi

while [[ $# -gt 0 ]]; do
  case $1 in
    -d|--data) #recommended
      INPUT_DATA=$2
      module load atools
      alog --state start
      source <(aenv --data $INPUT_DATA) # reads from a CSV if provided
      shift 2;;
    -b|--barcodes)
      BAR=$2
      shift 2;;
    -F|--fwd)
      FASTQF=$2
      shift 2;;
    -R|--rev)
      FASTQR=$2
      shift 2;;
    -*|--*)
      echo unknown argument $1
      exit 1;;
  esac
done

echo "Using $SLURM_NTASKS threads"
echo "Forward reads: $FASTQF"
echo "Reverse reads: $FASTQR"
echo "Barcode file: $BAR"

# prepare node env.
module load Python/3.11.5-GCCcore-13.2.0
source pre-processing-venv

######################
# Set-up environment #
######################

NEW_DIR=$(basename $BAR | cut -d. -f1) # Name of POOL/Barcode file

mkdir $NEW_DIR #create subdir to work in and cd
cd $NEW_DIR

echo "$(date) - [DEMULTIPLEX] $(basename $FASTQF)"
pwd;

## 0. prepare the barcode FASTA files, 1 for FWD and one for the REV barcodes
FWD_BAR=FWD_barcodes.fasta
REV_BAR=REV_barcodes.fasta
awk -F'\t' -v out1="$FWD_BAR" -v out2="$REV_BAR" '{
  print ">" $1 "\n" $2 > out1
  print ">" $1 "\n" $3 > out2
}' ../$BAR

#######################
### DEMULTIPLEX UDI ###
#######################

# Demultiplexing is done in 2 steps because the amplicon can be in "normal" or reverse complement orientation
## 1. First step is for the "normal" R1|FWD>>>>>REV|R2 reads
##   -> R1: FWD>>>>... and R2: VER<<<<...
echo "$(date) - [START] Demultiplexing in Direct orientation ..."

cutadapt -j $SLURM_NTASKS -e 0 --no-indels --pair-adapters \
  -g ^file:$FWD_BAR -G ^file:$REV_BAR \
  -o '{name}_R1.fastq.gz' -p '{name}_R2.fastq.gz' \
  --rename='{id};sample={adapter_name};' \
  ../$FASTQF ../$FASTQR &> "$NEW_DIR"_demultiplexing.log

echo "$(date) - [STOP] Demultiplexing in Direct orientation"

## 2. Second step is for the "reverse complement" R1|VER<<<<<DWF|R2 reads
## by reversing R1 and R2 we can get the reverse complement R2|FWD>>>>>REV|R1
##   -> R2: FWD>>>>... and R1: VER<<<<...
echo "$(date) - [START] Demultiplexing in Reverse orientation ..."

cutadapt -j $SLURM_NTASKS -e 0 --no-indels --pair-adapters \
  -g ^file:$FWD_BAR -G ^file:$REV_BAR \
  -o '{name}_rcR1.fastq.gz' -p '{name}_rcR2.fastq.gz' \
  --rename='{id};sample={adapter_name};' \
  ../$FASTQR ../$FASTQF &> "$NEW_DIR"_demultiplexing_rc.log

echo "$(date) - [STOP] Demultiplexing in Reverse orientation"

# deactive venv
deactivate

################
### CLEAN UP ###
################

## 1. Concatenate outputs from both directions for each sample
# Loop over all RC files and cat them to the "normal" ones
echo "$(date) - [START] Concatenating NORMAL and REVCOMP demultiplexed reads ..."

RC_SAMPLES_PE=( *_rcR[12].fastq.gz )

for RC_SAMPLE in "${RC_SAMPLES_PE[@]}"
do
 NORMAL_SAMPLE=`echo $RC_SAMPLE | sed -e 's/rc//g'`
 cat $RC_SAMPLE >> $NORMAL_SAMPLE # concatenate the revcomp results to the "normal" results, works on gzipped data as well!
 rm $RC_SAMPLE # remove the revcomp
done

echo "$(date) - [STOP] Concatenating NORMAL and REVCOMP demultiplexed reads"

##  2. Generate md5sum of data
echo "$(date) - [START] Calculating md5 hashes ..."
md5sum [!unknown]*.fastq.gz > "$NEW_DIR"_md5.txt
echo "$(date) - [STOP] Calculating md5 hashes done!"

## 3. Create file with sequence count per sample file
echo "$(date) - [START] Generating reads per sample overview ..."

seqkit stats -e -a -T -o stats.tsv -j $SLURM_NTASKS *.fastq.gz

echo "$(date) - [STOP] Generating reads per sample overview"

echo "$(date) - [FINISHED] $NEW_DIR DONE!"

if [ ! -z "$INPUT_DATA" ]; then
    alog --state end --exit $?
fi

