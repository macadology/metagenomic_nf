# Metagenomic Pipeline
Metagenomic pipeline for taxonomy and functional classification.

## Modes
The pipeline can be run in 2 modes. `--general true`
1. The first (specific) assumes that metagenomic_nf is in the project folder and local.config contains parameters specific to the project. The pipeline uses these parameters as default.
2. The second (general) assumes that metagenomic_nf is a generic pipeline and ignore local.config.

## Examples
```bash
#Jonai examples :
nextflow run main.nf -profile jonai #show the reads that'll be used as inputs
nextflow run main.nf -profile jonai -w [WORKDIR] # specify where to store work directory
nextflow run main.nf -profile jonai --readtype custom --querydir [Directory] --queryglob [glob pattern] #show reads queryglob and directory
nextflow run main.nf -profile jonai --readtype raw --profilers fastp,decont
nextflow run main.nf -profile jonai --readtype fastp --profilers decont
nextflow run main.nf -profile jonai --readtype decont --profilers kraken2,bracken --krakenKeepOutput true #Save all kraken output
nextflow run main.nf -profile jonai -w [WORKDIR] --krakenMMAP --profilers kraken2,bracken #Preload kraken database. See README.md
nextflow run main.nf -profile jonai --profilers fastp,kraken2,bracken,srst2,humann3
nextflow run main.nf -profile jonai --profilers align --bwaIndex [ref.fasta]
nextflow run main.nf -profile jonai --profilers align --bwaIndexDir [IndexDir] --bwaIndex [ref.fasta]

#AWSBatch examples :
nextflow run main.nf -profile batch -plugins nf-amazon --bucket-dir s3://jon-nextflow-work --profilers humann3

#ACRC examples :
nextflow run main.nf -profile acrc --querydir [Directory] --queryglob "*_{1,2}*{fastq,fastq.gz,fq,fq.gz}" --outputdir [Directory] -w [WORKDIR] --profilers humann3
```

## Kraken preload
For a faster kraken run on a local machine, increase size of `/dev/shm` in `/etc/fstab` and enable `--krakenMMAP`. The pipeline will copy the database to `/dev/shm`, run kraken with `--memory-mapping`, and remove the database after. I recommend not running other profilers while doing this since the ram will be full.
```bash
nextflow run main.nf -w [WORKDIR] --krakenMMAP --profilers kraken2,bracken
nextflow run main.nf --krakenMMAP --profilers kraken2,bracken --krakenKeepOutput true #Save all kraken output
```

## Folder structure
The pipeline assumes that the reads are present in separate folders in the `querydir` that are named by the sample prefix.

Example A
```bash
querydir
├── WEM001
│   ├── WEM001_1.fq.gz
│   └── WEM001_2.fq.gz
├── WEM002
│   ├── WEM002_1.fq.gz
│   └── WEM002_2.fq.gz
└── metadata.xlsx
```
Example B
```bash
querydir
├── WEM001
│   ├── fastp
│   │   ├── fastp_WEM001_1.fq.gz
│   │   └── fastp_WEM001_2.fq.gz
├── WEM002
│   ├── fastp
│   │   ├── fastp_WEM002_1.fq.gz
│   │   └── fastp_WEM002_2.fq.gz
└── metadata.xlsx
```
The output should mirror the input directory structure.
```bash
outputdir
├── WEM001
│   ├── kraken
│       ├── minikraken_8GB_20200312
│       │   ├── WEM001.bracken.F
│       │   ├── WEM001.bracken.G
│       │   ├── WEM001.bracken.P
│       │   ├── WEM001.bracken.S
│       │   ├── WEM001.kraken2.report
│       │   └── WEM001.kraken2.tax
├── WEM002
│   ├── kraken
│       ├── minikraken_8GB_20200312
│       │   ├── WEM002.bracken.F
│       │   ├── WEM002.bracken.G
│       │   ├── WEM002.bracken.P
│       │   ├── WEM002.bracken.S
│       │   ├── WEM002.kraken2.report
│       │   └── WEM002.kraken2.tax
└── metadata.xlsx
```
