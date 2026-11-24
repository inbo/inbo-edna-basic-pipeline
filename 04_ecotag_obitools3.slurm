#!/bin/bash -l

#SBATCH --clusters=wice
#SBATCH --partition=batch 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:15:00
#SBATCH --out=obi-ecotag-%j.out
#SBATCH --job-name=obi-ecotag

# This script is meant to be submitted to a HPC system with SLURM job sheduler.
# SBATCH arguments might need to be adapted to fit your system/account/task

#### INPUT:
# (-i|--input) IN_FASTA: sequences to be taxonomically assigned, zotus/ASV.fasta (FASTA)
# (-R|--refdb) DB: REF_DB (OBITools3 DMS)

#### EXAMPLE job submission:
# sbatch .../04_ecotag_obitools3.slurm -i otus.fasta -R ../MY_REF_DB.dms

#### References:
# - OBITools: Boyer F, Mercier C, Bonin A, Le Bras Y, Taberlet P, Coissac E: OBITools: a Unix-inspired software package for DNA metabarcoding. Mol Ecol Resour, 2016: 176-182

################################
# Commandline Argument Parsing #
################################

while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--input)
      IN_FASTA=$2
      shift 2;;
    -R|--refdb)
      DB=$2
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
module load Python/3.9.5-GCCcore-10.3.0
source obitools3-venv

echo "--- [START] OBITools3 pipeline part 2 (TAXONOMY) - $(date +"%T")"

DMS=OBI3_$SLURM_JOB_ID

REFDB_FILE=$(basename $DB)
TABOUT=$(basename $IN_FASTA .fasta)."${REFDB_FILE%.fa*}".table

######################################
### OBITools3 Taxonomic Annotation ###
######################################

# Import zOTU table in obi DMS
obi import $IN_FASTA $DMS/0_query_data

# clean ref_DB -> if a job fails for wathever reason this gets corrupted so just clean before using
obi clean_dms $DB

# Tax assignment
obi ecotag --taxonomy $DB/taxonomy/taxdump -R $DB/ecopcr_final_0.995 $DMS/0_query_data $DMS/1_ecotag

# Generate some stats in STDOUT
obi stats -c SCIENTIFIC_NAME $DMS/1_ecotag

# Decorate with taxonomy
obi annotate --with-taxon-at-rank phylum --with-taxon-at-rank class --with-taxon-at-rank order --with-taxon-at-rank family \
    --with-taxon-at-rank genus --with-taxon-at-rank species --taxonomy $DB/taxonomy/taxdump \
    $DMS/1_ecotag $DMS/2_ecotag_clean

# Export
obi export --tab-output $DMS/2_ecotag_clean -o $TABOUT

echo "--- [STOP] OBITools3 pipeline part 2 (DONE) - $(date +"%T")"
obi ls $DMS

# clean ref_DB -> if a job fails for wathever reason this gets corrupted so also clean after using
obi clean_dms $DB
obi clean_dms $DMS
