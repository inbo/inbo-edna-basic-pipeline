# inbo-edna-basic-pipeline


This repository contains the basic bioinformatics steps used in most INBO eDNA pipelines to process raw Illumina paired-end sequencing data to Operational Data (OTU table with taxonomy) for use in downstream data analysis.

This pipeline is part of the INBO eDNA metabarcoding workflow(s) described in https://inbo.github.io/edna-workflows/

> [!IMPORTANT]
> Currently this repository is targeted to INBO users with access to the Tier-2 HPC of the KUL-VSC.
> In principle the described pipelines can be easily adapted to run on other machines: the methodology is straightforward and the used software is well documented, free and open source.
> See [Dependencies](#dependencies) and [References](#references)


## Workflow description

The basic eDNA workflow consists of:

`01_demultiplex.slurm`

- Demultiplexing

`02_preprocess.slurm`

- Paired-end read merging
- Primer trimming
- Basic read QC filter
- Dereplication (per sample)

`03_cluster_vsearch.slurm`

- Dereplication (global)
- Clustering/denoising
- Chimera removal

`04_ecotag_obitools3.slurm`

- Taxonomic annotation

## Preparation and set-up (INBO users)

- Clone this repo to an appropriate location

```
git clone https://github.com/inbo/inbo-edna-basic-pipeline.git
```

- Add the the INBO shared software folder to your `$PATH` (this also contains the needed venvs)

- Submit the `.slurm` files as a batch job with SLURM, e.g.

``` 
sbatch --account=<credit_acc> <other SLURM args> <path_to>/inbo-edna-basic-pipeline/<myjob>.slurm <CMD line args>
```


> [!NOTE]
> Output files (results and log-files) will be created in the directory you launch the job from!

## Run examples

### 1. Demultiplexing

- Per pooled library, create a "barcode file" (TSV) that contains the barcode indexes to identify each sample. 
The TSV should have 3 tab-separated columns with no header: `<sample_id> <FWD-index> <REV-index>`

```
Sample1    AACAA   GGTCGGTGTA
Sample2    TTGTT   ACCTCCAAGC
Sample3    AATCC   CGTAGGTGTT
[...]
```

> [!NOTE]
> You can use the `report.samplesheet()` function from https://github.com/inbo/gendiver

- Create a `data.csv` with 3 columns -> "BAR","FASTQF","FASTQR". The data in the rows correspond to the pools (FASTQ files) you need to demultiplex.

```
"BAR","FASTQF","FASTQR"
P1_barcodes.tsv,P1_R1.fastq.gz,P1_R2.fastq.gz
P2_barcodes.tsv,P2_R1.fastq.gz,P2_R2.fastq.gz
P3_barcodes.tsv,P3_R1.fastq.gz,P3_R2.fastq.gz
[...]
```

- Submit as a job array (e.g. for 5 pooled libs):

```
sbatch --array=1-5 01_demultiplex.slurm -d data.csv
```

> [!NOTE]
> In the default INBO metabarcoding workflow fragments are **not** cloned directionally. 
> During demultiplexing reads are checked in both directions for the indexes.

### 2. Preprocessing

- Choose and copy a configuration file from the `config/` dir
- Create a `data.csv` with 3 columns -> "FASTQF","FASTQR","CONFIG". The data in the rows correspond to the samples to be processed.

```
"FASTQF","FASTQR","CONFIG"
S1_R1.fastq.gz,S1_R2.fastq.gz,config.json
S2_R1.fastq.gz,S2_R2.fastq.gz,config.json
S3_R1.fastq.gz,S3_R2.fastq.gz,config.json
[...]
```

- Submit as a job array (e.g. for 96 samples):

```
sbatch --array=1-96 02_preprocess.slurm -d data.csv
```

### 3. Clustering/denoising

- Pre-processed files need to be concatenated, e.g. use `cat *.uniq.fasta.gz > all_uniqs.fasta.gz`
- Then denoise into zOTUs (and remove singletons prior)

```
sbatch 03_cluster_vsearch.slurm -i all_uniqs.fasta.gz -r 2
```

### 4. Taxonomic annotation

- Requirement: an OBITools3 DMS with reference sequences

> [!NOTE]
> Methodology to create a custom reference database with OBITools: https://github.com/inbo/prj_edna_refdb

```
sbatch 04_ecotag_obitools3.slurm -i zotus.fasta -R refdb.obidms
```

### Extra: follow along in R

In https://github.com/inbo/gendiver there are some utility functions to interact with the pipeline output.

#### Quality Control

You can generate QC reports for each step of the pipeline using the `gendiver::report.*` functions

#### Operational Data

You can read the OTU table and taxonomy file to a data.frame in R using the `gendiver::read.otu_table*` and `gendiver::read.taxonomy*` functions

### Dependencies

#### Software stack

The basic bioinformatics pipeline is run on a HPC system with pre-installed software modules. 
The main used software modules at VSC KUL Tier-2 are:

- `atools/1.5.1.lua`
- `Python/3.11.5-GCCcore-13.2.0` (cutadapt)
- `VSEARCH/2.30.0-GCC-13.2.0`
- `Python/3.9.5-GCCcore-10.3.0` (obitools3)

#### Pyhton virtual environments

Based on both python versions we created 2 virtual environments. These are shared in the INBO project directory (`stg_00184`). 
You can also find the required software for the virtual envs in the respective `requirements.txt` files.

#### Shared binaries

INBO users can add the path to the shared software installations in the project dir to their `$PATH` variable to have all other software available.

### References

**VSEARCH**

https://github.com/torognes/vsearch

Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. doi: 10.7717/peerj.2584

See also the references in: https://github.com/torognes/vsearch#references

**OBITools3**

https://forge.metabarcoding.org/obitools/obitools3

Boyer F, Mercier C, Bonin A, Le Bras Y, Taberlet P, Coissac E: OBITools: a Unix-inspired software package for DNA metabarcoding. Mol Ecol Resour, 2016: 176-182

**Cutadapt**

https://github.com/marcelm/cutadapt/

Marcel Martin. 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. https://doi.org/10.14806/ej.17.1.200

**SeqKit** 

https://github.com/shenwei356/seqkit

Wei Shen*, Botond Sipos, and Liuyang Zhao. 2024. SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing. iMeta e191. doi:10.1002/imt2.191.
