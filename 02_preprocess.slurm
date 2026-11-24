#!/bin/bash -l

#SBATCH --clusters=wice
#SBATCH --partition=batch
#SBATCH --account=lp_inbo
#SBATCH --ntasks=6
#SBATCH --nodes=1
#SBATCH --time=00:20:00
#SBATCH --out=preprocess-%j.out
#SBATCH --job-name=preprocess-pipeline

# This script is meant to be submitted to a HPC system with SLURM job sheduler.
# SBATCH arguments might need to be adapted to fit your system/account/task

# Author:	Annelies Haegeman - Flanders Research Institute for Agriculture, Fisheries and Food, Melle, Belgium
#  - https://gitlab.com/ahaegeman/eDNA_INBO-ILVO
# Author: Nick Dillen - Research Institute for Nature and Forest (INBO), Belgium
#	Extensively modified by: Io Deflem - Research Institute for Nature and Forest (INBO), Belgium
# Maintainer: Nick Dillen

# This script should be run on an HPC system with SLURM job sheduler. 
# It pre-processes metabarcoding data. Input is Illumina paired-end (demultiplexed) data. 
# It uses several third party software tools, which should be installed on your system:
# 1. vsearch (https://github.com/torognes/vsearch) (Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. doi: 10.7717/peerj.2584)
# 2. cutadapt (https://cutadapt.readthedocs.io/en/stable/guide.html) (https://doi.org/10.14806/ej.17.1.200)	
# 3. seqkit (Wei Shen*, Botond Sipos, and Liuyang Zhao. 2024. SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing. iMeta e191. doi:10.1002/imt2.191.)
# 4. (optional) atools (https://github.com/gjbex/atools)

#### INPUT:
## Job array:
# (-d|--data) INPUT_DATA: path to data.csv with the variables as column names and corresponding values as rows (see atools)
# The required columns are: "FASTQF","FASTQR","CONFIG"

## Command-line arguments take precedence over parameters specified with (-d) INPUT_DATA and atools!
# (-F|--fwd) FASTQF: path to forward reads FASTQ of sample
# (-R|--rev) FASTQR: path to reverse reads FASTQ of sample
# (-c|--config) CONFIG: path to config.json specifying the parameterization for the pipeline (JSON)

#### RUN description
# - Preparation: create data.csv file with input for pre-processing script => "FASTQF,FASTQR,CONFIG"
# - Launch job, see example below.
# Procedure:
#    1. Merging of forward and reverse reads using vsearch --fastq_mergepairs (https://github.com/torognes/vsearch)
#    2. Removal of primers using cutadapt (https://cutadapt.readthedocs.io/en/stable/guide.html)
#    3. Quality filtering using `--fastq_filter` from vsearch
#    4. Dereplication using `--fastx_uniques` from vsearch
# - For each step a dedicated folder will be created to store all output.

#### EXAMPLE job submission:
# sbatch --array=1-200%50 02_preprocess.slurm -d input_data.csv
#### You can also submit the script with command line arguments:
# sbatch 02_preprocess.slurm -F R1.fastq -R R2.fastq -c config.json
## Or overwrite/set one of the args by hand:
# sbatch --array=1-200%50 02_preprocess.slurm -d input_data.csv -c custom_config.json


##########################
### 0. Read user-input ###
##########################
set -e # exit on error
module purge

# in case of no arguments, show help message
if [ $# == 0 ]; then
    printf "
    
[ERROR] No input found!

Usage: sbatch 02_preprocess.slurm -F R1.fastq -R R2.fastq -c config.json
Where:
    -F : input FASTQ file which contains the demultiplexed FORWARD sequences
    -R : input FASTQ file which contains the demultiplexed REVERSED sequences
	  -c : config file to use in pre-processing (JSON)

Alternatively, use an array job to process multiple samples in parallel (recommended)
Usage: sbatch --array=1-200%50 02_preprocess.slurm -d input_data.csv
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
    -c|--config)
      CONFIG=$2
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

# READ parameters from CONFIG FILE
#mergepairs
MERGE_MIN_OVERLAP=$(cat $CONFIG | jq -r '.["fastq_mergepairs"]["fastq_minovlen"] // ""')
MERGE_TR_QUAL=$(cat $CONFIG | jq -r '.["fastq_mergepairs"]["fastq_truncqual"] // ""')
MERGE_MIN_LEN=$(cat $CONFIG | jq -r '.["fastq_mergepairs"]["fastq_minmergelen"] // ""')
MERGE_MAX_LEN=$(cat $CONFIG | jq -r '.["fastq_mergepairs"]["fastq_maxmergelen"] // ""')
MERGE_MAX_N=$(cat $CONFIG | jq -r '.["fastq_mergepairs"]["fastq_maxns"] // ""')
#cutadapt
CUTADAPT_REV_PRIMER=$(cat $CONFIG | jq -r '.["cutadapt"]["adapter"] // ""')
CUTADAPT_FWD_PRIMER=$(cat $CONFIG | jq -r '.["cutadapt"]["front"] // ""')
CUTADAPT_DISCARD_UT=$(cat $CONFIG | jq -r '.["cutadapt"]["discard-untrimmed"] // ""')
CUTADAPT_ERROR_RATE=$(cat $CONFIG | jq -r '.["cutadapt"]["error-rate"] // ""')
#fastq_filter
FILTER_EE=$(cat $CONFIG | jq -r '.["fastq_filter"]["fastq_maxee"] // ""')
FILTER_MAX_N=$(cat $CONFIG | jq -r '.["fastq_filter"]["fastq_maxns"] // ""')
FILTER_MIN_LEN=$(cat $CONFIG | jq -r '.["fastq_filter"]["fastq_minlen"] // ""')
FILTER_MAX_LEN=$(cat $CONFIG | jq -r '.["fastq_filter"]["fastq_maxlen"] // ""')

######################
# Set-up environment #
######################

# load modules (KUL HPC specific)
module load Python/3.11.5-GCCcore-13.2.0
module load VSEARCH/2.30.0-GCC-13.2.0

# make output dirs envorinment variables 
export MERGED_DIR=./01_merged
export TRIMMED_DIR=./02_trimmed
export QF_DIR=./03_quality_filtered
export DEREP_DIR=./04_dereplicated

# Create DIRs (if not exist)
mkdir -p $MERGED_DIR
mkdir -p $TRIMMED_DIR
mkdir -p $QF_DIR
mkdir -p $DEREP_DIR

# Define sample name for worker -> basename for all intermediate files
SAMPLE_B=`basename $FASTQF`
export SAMPLE="${SAMPLE_B%_*.fast*}"

echo "Using $SLURM_NTASKS threads"
echo "Forward reads: $FASTQF"
echo "Reverse reads: $FASTQR"
echo "Config file: $CONFIG"
echo "[optional] Input data file: $INPUT_DATA"

echo "--- [START] pre-processing of sample: $SAMPLE - $(date +"%T")"

# $SAMPLE variable from input or export in job-script
printf "\nPROCESSING sample $SAMPLE \n"

#################################
### 1. Merge paired-end reads ###
#################################

printf "$(date +"%T") - [START] merging reads for $SAMPLE....\n"

## Decide pear or vsearch: If both PEAR and vsearch mergepairs params are given -> vsearch gets priority
if [[ -n $(cat $CONFIG | jq -r '.["fastq_mergepairs"] // ""') ]];
  then
  printf "[INFO] using VSEARCH --fastq_mergepairs based on config input (priority)\n"
  
  vsearch --fastq_mergepairs "$FASTQF" --reverse "$FASTQR" --fastqout - --fastq_allowmergestagger \
    ${MERGE_MIN_OVERLAP:+--fastq_minovlen="$MERGE_MIN_OVERLAP"} \
    ${MERGE_TR_QUAL:+--fastq_truncqual="$MERGE_TR_QUAL"} \
    ${MERGE_MIN_LEN:+--fastq_minmergelen="$MERGE_MIN_LEN"} \
    ${MERGE_MAX_LEN:+--fastq_maxmergelen="$MERGE_MAX_LEN"} \
    --threads $SLURM_NTASKS 2> $MERGED_DIR/"$SAMPLE".mergepairs.log | 
    seqkit seq -o $MERGED_DIR/"$SAMPLE".merged.fastq.gz -w 0
else
printf >&2 "[ERROR] No config.json fields found for paired-end read merging (REQUIRED)\n"; exit 1;
fi

printf "$(date +"%T") - [STOP] done merging reads for $SAMPLE.\n"

#######################
### 2. Trim primers ###
#######################

# load python venv for cutadapt
source pre-processing-venv

printf "$(date +"%T") - [START] removing primers for $SAMPLE....\n"
#remove both FWD and REV primers using cutadapt.
shopt -s extglob        # Enable extended globs (bash)

# -g: Specify required linked adapters
# -e: This sets the error tolerance used when searching for adapters.
# --discard-untrimmed: Discard reads in which no adapter match was found.
cutadapt -g "^$CUTADAPT_FWD_PRIMER...$CUTADAPT_REV_PRIMER$" \
  ${CUTADAPT_DISCARD_UT:+--discard-untrimmed} \
	${CUTADAPT_ERROR_RATE:+--error-rate="$CUTADAPT_ERROR_RATE"} \
  --cores $SLURM_NTASKS	-o $TRIMMED_DIR/"$SAMPLE".fullytrimmed.fastq.gz \
  $MERGED_DIR/"$SAMPLE"@(.assembled|.merged).fastq.gz &> $TRIMMED_DIR/"$SAMPLE".fullytrimmed.log

printf "$(date +"%T") - [STOP] done removing primers for $SAMPLE.\n"

########################################
### 3. Quality statistics and filter ###
########################################

printf "$(date +"%T") - [START] calculating quality statistics and filtering $SAMPLE....\n"
# -fastq_stats: generate statistics.
# -fastq_filter: filter reads.

vsearch --fastq_stats $TRIMMED_DIR/"$SAMPLE".fullytrimmed.fastq.gz --log $QF_DIR/"$SAMPLE".fastqstats.log

vsearch --fastq_filter $TRIMMED_DIR/"$SAMPLE".fullytrimmed.fastq.gz --fastaout - \
	${FILTER_EE:+--fastq_maxee="$FILTER_EE"} \
	${FILTER_MAX_N:+--fastq_maxns="$FILTER_MAX_N"} \
	${FILTER_MIN_LEN:+--fastq_minlen="$FILTER_MIN_LEN"} \
	${FILTER_MAX_LEN:+--fastq_maxlen="$FILTER_MAX_LEN"} \
	2> $QF_DIR/"$SAMPLE".fastqfilter.log |
	seqkit seq -o $QF_DIR/"$SAMPLE".qf.fasta.gz -w 0

printf "$(date +"%T") - [STOP] done calculating quality statistics and trimming $SAMPLE.\n"

###################################
### 3. Per Sample Dereplication ###
###################################

printf "$(date +"%T") - [START] dereplicating  $SAMPLE....\n"

# --relabel: to add file/samplename in header
vsearch --fastx_uniques $QF_DIR/"$SAMPLE".qf.fasta.gz --fastaout - \
  --sizeout --relabel "sample=$SAMPLE;uniq_id=" 2> $DEREP_DIR/"$SAMPLE".fastxuniques.log | 
  seqkit seq -o $DEREP_DIR/"$SAMPLE".uniq.fasta.gz -w 0

printf "$(date +"%T") - [STOP] dereplicating $SAMPLE.\n"

###################################
# deactivate venv
deactivate

printf "\nDONE PROCESSING SAMPLE $SAMPLE.\n"

echo "--- [STOP] pre-processing of sample: $SAMPLE - $(date +"%T")"

echo "[DONE] $(date)"

if [ ! -z "$INPUT_DATA" ]; then
    echo $INPUT_DATA
    alog --state end --exit $?
fi

