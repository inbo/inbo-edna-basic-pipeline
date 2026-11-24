#!/bin/bash -l

#SBATCH --clusters=wice
#SBATCH --partition=batch 
#SBATCH --nodes=1
#SBATCH --ntasks=18
#SBATCH --time=01:00:00
#SBATCH --out=vsearch-OTU-%j.out
#SBATCH --job-name=vsearch-OTU-pipeline

# This script is meant to be submitted to a HPC system with SLURM job sheduler.
# SBATCH arguments might need to be adapted to fit your system/account/task

#### INPUT:
# (-i|--input) IN_FASTA: all_samples.fasta, should have ";sample=" column in the headers (see vsearch manual)
# (-d|--id) ID: similarity to cluster on. If empty -> no OTU clustering is performed and you get a zOTU table, default
# (-l|--min-len) MINL: minimum length of the amplicon, filter
# (-L|--max-len) MAXL: maximum length of the amplicon, filter
# (-r|--min-reads) MINR: minimum number of reads to retain UNIQUE sequence, filter

#### EXAMPLE job submission:
# sbatch .../03_clustering_vsearch.slurm --input input.fasta.gz -l 100 -L 200 -r 2
# sbatch .../03_clustering_vsearch.slurm --input input.fasta.gz -l 100 -L 200 -r 2 -d 0.95

### References:
# Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. doi: 10.7717/peerj.2584
# https://github.com/torognes/vsearch#references

################################
# Commandline Argument Parsing #
################################

#defaults
MINL=0
MAXL=99999
MINR=1
IN_FASTA="input.fasta"

while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--input)
      IN_FASTA=$2
      shift 2;;
    -d|--id)
      ID=$2
      shift 2;;
    -l|--min-len)
      MINL=$2
      shift 2;;
    -L|--max-len)
      MAXL=$2
      shift 2;;
    -r|--min-reads)
      MINR=$2
      shift 2;;
    -*|--*)
      echo unknown argument $1
      exit 1;;
  esac
done

######################
# Set-up environment #
######################
set -e

module purge
module load VSEARCH/2.30.0-GCC-13.2.0

printf "[START] $(date)\n"
printf "[INFO] Using input file: $IN_FASTA \n"

printf "[INFO] Using $SLURM_NTASKS cores \n"

INBASENAME=$(basename $IN_FASTA | cut -d. -f1)
BASEN="$INBASENAME"_minL"$MINL"_maxL"$MAXL"_minR"$MINR"

UNIQ_OUT="$BASEN.uniques.fasta"
DENOISED_OUT="$BASEN.denoised.fasta"
ZOTUS_OUT="$BASEN.zotus.fasta"
OTU_OUT="$BASEN.otus_$ID.fasta"

### Per sample dereplication
# check if derep input:
if [[ $(seqkit head $IN_FASTA | seqkit grep -r -C -p "size=") -gt 1 && $(seqkit head $IN_FASTA | seqkit grep -r -C -p "sample=") -gt 1 ]]; then
  printf "\n[INFO] input is already dereplicated, proceed! ";date 
else
  printf "\n[ERROR] Incorrect input file format!
  Input file needs to be a FASTA file, dereplicated per sample, using 'sample=' and 'size=' tags in the FASTA header to specify unique seqs per sample.\n";date
  exit 1
fi

###########################
### Clustering pipeline ###
###########################

# First dereplicate the input sequences into unique seqs
printf "\n[START] GLOBAL dereplicating ... ";date
vsearch --fastx_uniques $IN_FASTA --fastaout $UNIQ_OUT --sizein --sizeout --relabel_sha1 \
  --minuniquesize "$MINR" --minseqlength "$MINL" --maxseqlength "$MAXL" --fasta_width 0
printf "[DONE] GLOBAL dereplicating finished! ";date

# Then apply denoising algorithm (see UNOISE3)
printf "\n[START] denoising ... ";date
vsearch --cluster_unoise $UNIQ_OUT --centroids $DENOISED_OUT --sizein --sizeout \
  --minsize "$MINR" --minseqlength 1 --threads $SLURM_NTASKS --fasta_width 0
printf "[DONE] denoising finished! ";date

# Need to manually launch chimera detection (in usearch this is contained in UNOISE3)
printf "\n[START] detecting chimeras ... ";date
vsearch --uchime3_denovo $DENOISED_OUT --nonchimeras $ZOTUS_OUT --sizein --sizeout --fasta_width 0
printf "[DONE] chimeras removed! ";date

# Check from user input if need to cluster into OTUs or not

if ! [[ -z "${ID}" || "${ID}" == "1" ]]; then
    printf "\n[START] Clustering zOTUs into OTUs with $ID similarity ";date

    # create OTUs with similarity threshold = $ID and write to output table
    vsearch --cluster_size $ZOTUS_OUT --id "$ID" --centroids $OTU_OUT --sizein --sizeout \
      --threads $SLURM_NTASKS --minseqlength 1 --fasta_width 0

    printf "[DONE] Clustering zOTUs into OTUs with $ID similarity ";date
    printf "\n[START] generating OTU table with --usearch_global ";date

    # match OTUs to original sequence data by alignment -> get table with OTU distributions per sample
    TAB_OUT="$BASEN.otutab_"$ID".table"
    vsearch --usearch_global $IN_FASTA --db $OTU_OUT --id "$ID" --notrunclabels --strand plus --otutabout $TAB_OUT --threads $SLURM_NTASKS

    printf "[DONE] OTU table created! ";date

else
    printf "\n[INFO] Skipping OTU clustering as no ID was provided"
    printf "\n[START] generating zOTU table with --search_exact ";date

    # Finally generate a zOTU table by looking up how many times each zOTU is present per sample (search_exact) in the original input fasta
    TAB_OUT="$BASEN.zotutab.table"
    vsearch --search_exact $IN_FASTA --db $ZOTUS_OUT --otutabout $TAB_OUT --notrunclabels --threads $SLURM_NTASKS

    printf "[DONE] zOTU table created! ";date
fi

printf "[INFO] Output file: "$TAB_OUT" \n"

printf "\n[DONE] $(date)"
