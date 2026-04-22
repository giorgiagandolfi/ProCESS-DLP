# Tutorial ProCESS-DLP 

This is an easy to use tutorial to reproduce single-cell DNA sequencing data. 

## 0. Set-up
In order to perform a DLP-like ProCESS simulation, you need to install a specific ProCESS version from a specific branch.

```r
devtools::install_github("caravagnalab/ProCESS", ref="on-the-fly_mutations")
```

Or either using the docker image with singularity:

```bash
singularity pull docker://giorgiagandolfi97/process_on_the_fly:v2
```

## 1. Simulate tissue and sampling

Run the Rscript in `scripts/0_simulation/ProCESS_dlp_sim7.R`. This scripts will generate a three clone CRC tumour. Each clone is
characterized by a specific growht rate and a fixed SNV/INDEL rate. CN rate is changing across clones. A total of 900 cells are sampled.
You can run the script by submitting the `scripts/0_simulation/run_process_sim.sh` to the SLURM jobs scheduler.

## 2. Simulate sequencing

Tumour cells are sampled as a DLP+ like experiment, so each cell is sequenced separately. On the other hand normal sample is sequenced as bulk. The code is optimized differently:

- tumour: single cell sequencing is a single job (job array lenght = n° of total cells);
- normal: the total bulk coverage is split into _N_ chunks of subcoverage, of 5X. Each job is a bulk normal sequencing of the selected subcoverage (job array lenght = n° of chunks). Expected total coverage is 100x, total of 20 chunks of 5x coverage.

### 2.1. Single cell tumour sequencing

Run the `sh` script that you find in `scripts/1_sequencing/run_sequencing.sh`.This script will perform sequencing of 
sampled cells in parallel as a job array. 

### 2.2. Bulk normal sample sequencing

Run the `sh` script that you find in `scripts/1_sequencing/run_sequencing_normal.sh`.This script will perform sequencing of 
subsampled chunks in parallel as a job array. 

## 3. Preprocess SAM files and generate FASTQs

Similary as the sequencing procedure, also for the preprocessing you have to run two different steps for normal bulk and tumour single cells. 
Run the `sh` script that you find in `scripts/2_preprocessing/run_preprocessing.sh` for tumour sample and the `scripts/2_preprocessing/run_preprocessing_normal.sh`. This script will:

1. merge all SAM files for each chromosome
2. create FASTQs from SAM files
3. save an `rds` of the ground truth sequenced only mutations, with `DP`, `NV`, `VAF` and classification.

The output folder will contain paired end FASTQs files for each sampled cell and for each normal subchunk.

### Output folder

An example of the output for a simulation comprising:

- a bulk 2x normal sample (splitted into 2 chunks of 1x each)
- two single cells from the simulation

```bash
sequencing_dlp7/
├── normal_sample
│   ├── sequencing_chunk01
│   │   ├── chunk01.R1.fastq.gz
│   │   ├── chunk01.R2.fastq.gz
│   │   ├── chunk01.singleton.fastq.gz
│   │   ├── chunk01.unpaired.fastq.gz
│   │   └── seq_res.rds
│   └── sequencing_chunk02
│       ├── chunk02.R1.fastq.gz
│       ├── chunk02.R2.fastq.gz
│       ├── chunk02.singleton.fastq.gz
│       ├── chunk02.unpaired.fastq.gz
│       └── seq_res.rds
├── sequencing_DPL7_A_439-450
│   ├── DPL7_A_439-450.R1.fastq.gz
│   ├── DPL7_A_439-450.R2.fastq.gz
│   ├── DPL7_A_439-450.singleton.fastq.gz
│   ├── DPL7_A_439-450.unpaired.fastq.gz
│   └── seq_res.rds
└── sequencing_DPL7_A_440-450
    ├── DPL7_A_440-450.R1.fastq.gz
    ├── DPL7_A_440-450.R2.fastq.gz
    ├── DPL7_A_440-450.singleton.fastq.gz
    ├── DPL7_A_440-450.unpaired.fastq.gz
    └── seq_res.rds
```
