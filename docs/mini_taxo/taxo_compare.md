---
layout: default
title: "Comparing between Kraken2 and metaphlan3"
author: "Ngoc Giang Le"
version: 0.1
#bibliography:
nav_order: 1
description: "Taxonomic tutorial"
permalink: /
output:
    html_document:
      toc: true
      toc_float: true
      toc_depth: 3
      fig_caption: true
      highlight: tango
---


```{=html}
<style>
body {
text-align: justify}
</style>
```

In this tutorial, we will use Bash and Snakemake to compare two different taxonomic typing programs: Kraken2/Bracken vs Metaphlan3.
These programs use ...
look at the relative abundance of raw metagenomic.

For this exercise, we also look at the effect of threads on the analysis time.

## Setup new project

Create new directory for this project.

```
mkdir taxo_compare
```

Move the downloaded SRA called ERR4097276 from [previous tutorial]() to the new directory.

```
mv ERR4097276.sra* taxo_compare
```

## Create running file

We edit the `samples.tsv`, so it contains the path for the forward (r1) and the reverse (r2) files.

|ID       |r1                       |r2                       |
|---------|-------------------------|-------------------------|
|MetaHIT5 |ERR4097276.sra_1.fastq.gz|ERR4097276.sra_2.fastq.gz|


## Create Snakefile

Open your favorite editor to create a file called `Snakefile`.

We will use pandas to process the input `samples.tsv`.


```
import pandas as pd

def parse_samples(samples_tsv):
    # Load in tab separate file
    # Remove samples with incomplete fields
    # Set ID as the index
    return pd.read_csv(samples_tsv, sep ='\t').dropna().set_index("ID", drop=False)

def get_files(sample_df, wildcards, col):
    # Return forward and reverse reads based on the sample's name
    return sample_df.loc[wildcards.sample, [col]]

_samples = parse_samples("samples.tsv")
```


### Kraken2/Bracken database

These rules will download and setup the Kraken2/Bracken database from the indexeszone (citation).


```
rule all:
    input:
        "databases/hash.k2d",
        "databases/opts.k2d",
        "databases/taxo.k2d",
        "databases/seqid2taxid.map",

rule kraken_db:
    output:
        temp("databases/k2_standard_20220607.tar.gz")
    params:
        directory("databases")
    shell:
        """
        wget -P {params} https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20220607.tar.gz
        """
        
rule kraken_unzip:
    input:
        "databases/k2_standard_20220607.tar.gz"
    output:
        "databases/hash.k2d",
        "databases/opts.k2d",
        "databases/taxo.k2d",
        "databases/seqid2taxid.map",
    params:
        directory("databases")
    shell:
        """
        tar xf {input} -C {params}
        """
```

We need to write rule `all` to tell Snakemake that we would like to have the database ready.

Note: This data set is 60 GB in size.


Snakemake automatically creates the `databases/` directory.


### Running Kraken2/Bracken

Let's write rule to run the program.



Note: Kraken loads in databases so for smaller input files it takes longer.
Only when large raw file then kraken is faster.


Snakemake script

```
rule all:
    input:
        expand("01.taxonomy/raw/{sample}_@.mp3.profile", sample = _samples.index),
        expand("01.taxonomy/raw/{sample}_@.bracken", sample = _samples.index)
Running script
```

The `samples.tsv` contains the following information.


Analysis script

Conclusion:



(Rarefraction tutorial ???)



Advance testing of snakemake

different scripts


## Result


## The accuracy between different databases





https://www.biostars.org/p/10756/


https://software.cqls.oregonstate.edu/updates/metaphlan-3.0.14/index.md




Conda fix:


Note:

This was my fix to install metaphlan3 version 3.10

`vim ~/.condarc`

```
channel_priority: disabled
channels:
  - conda-forge
  - defaults

```

### Split the fastq pairs

The raw is around 4.5 Gb per pair, which we will split into 10 sets.


```
# Create conda environment
conda create --name seqkit -c bioconda seqkit

# Activate the environment
conda activate seqkit

# Split into 10 sets
# To split based on number of reads replace -p 10 with -s <number_of_reads>
seqkit split2 -1 ERR4097276.sra_1.fastq.gz -2 ERR4097276.sra_2.fastq.gz -p 10 -O ERR4097276 -f -e .gz

conda deactivate

```

[image]


```
[INFO] flag -1/--read1 and -2/--read2 given, ignore: -
[INFO] split seqs from ERR4097276.sra_1.fastq.gz and ERR4097276.sra_2.fastq.gz
[INFO] split into 10 parts
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_2.part_001.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_2.part_002.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_2.part_003.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_2.part_004.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_2.part_005.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_2.part_006.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_2.part_007.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_2.part_008.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_2.part_009.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_2.part_010.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_1.part_001.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_1.part_002.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_1.part_003.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_1.part_004.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_1.part_005.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_1.part_006.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_1.part_007.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_1.part_008.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_1.part_009.fastq.gz
[INFO] write 5718873 sequences to file: ERR4097276/ERR4097276.sra_1.part_010.fastq.gz
```

