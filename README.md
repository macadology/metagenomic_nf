# Metagenomic Pipeline
Metagenomic pipeline for taxonomy and functional classification.

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
