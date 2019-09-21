[![Build Status](https://travis-ci.com/maxibor/dada2-nf.svg?token=pwT9AgYi4qJY4LTp9WUy&branch=master)](https://travis-ci.com/maxibor/dada2-nf)

# dada2-nf

Simple [Dada2](https://benjjneb.github.io/dada2/index.html) [Nextflow](https://www.nextflow.io) pipeline

## Dependancies

- [conda](https://conda.io/en/latest/) 
- [Nextflow](https://www.nextflow.io/) : `conda install -c bioconda nextflow`

## Usage

```
nextflow run maxibor/dada2-nf --reads "/path/to/paired_end_reads_*.{1,2}.fastq.gz"
```

### Input

#### --reads

Use this to specify the location of your input FastQ files. For example:

- `--reads  'path/to/data/sample_*_{1,2}.fastq'`
### --silva_db   

Silva database for dada2. 

Example:

`--silva_db ./silva_nr_v132_train_set.fa`

Default: Automatically downloads the database [here](https://doi.org/10.5281/zenodo.1172783)


### --silva_species_db  

Silva species db for dada2.

Example:

`--silva_species_db silva_species_assignment_v132.fa`

 Default: Automatically downloads the database [here](https://doi.org/10.5281/zenodo.1172783)

### --rank

Taxonomic rank to keep. Either "Genus" or "Species"

Example:

`--rank Species`

Default: "Genus"

**Please note the following requirements:**

- The path must be enclosed in quotes
- The path must have at least one * wildcard character
- When using the pipeline with paired end data, the path must use {1,2} notation to specify read pairs.


### Output


## Help

```
$ nextflow run maxibor/dada2-nf --help
N E X T F L O W  ~  version 19.07.0
Launching `maxibor/dada2-nf` [cheesy_wright] - revision: 8a4d08f663 [master]
dada2: simple dada2 16s classifier pipeline
 Homepage: https://github.com/maxibor/dada2-nf
 Author: Maxime Borry <borry@shh.mpg.de>
=========================================
Usage:
The typical command for running the pipeline is as follows:
nextflow run maxibor/dada2-nf --reads '/path/to/paired_end_reads_*.{1,2}.fastq.gz'
Mandatory arguments:
  --reads                       Path to input data (must be surrounded with quotes)

Settings:
  --phred                       Specifies the fastq quality encoding (33 | 64). Defaults to 33
  --pairedEnd                   Specifies if reads are paired-end (true | false). Default = true
  --silva_db                    Silva database for dada2. Default =
  --silva_species_db            Silva species db for dada2. Default =
  --rank                        Taxonomic rank to retain (Genus | Species). Default = Genus


Options:
  --results                     The output directory where the results will be saved. Defaults to results
  --help  --h                   Shows this help page
```
