---
layout: default
title: Main
author: "Ngoc Giang Le"
version: 0.1
nav_order: 1
description: "Main page of the sequencing analysis tutorial"
permalink: /
output: html
---


# [Sequencing analysis of the ZymoBIOMICS Mock community](https://github.com/GiangLeN/Zymo-Mock-sequencing.git)

by [Ngoc Giang Le](https://github.com/GiangLeN)


This project shows how to build pipeline to analyze the metagenomic of publicly available zymoBIOMICS Mock community (citation). There are 8 bacteria and 2 fungi at specific concentration in this mock.

The aim of this tutorial is to look at different tools and approaches and compared them for benchmark purposes. Codes, and scripts are also provided to tie various tools together.

There are several topics to be covered in this tutorial:

- Download data from NCBI
- Environment setup
- Build pipelines for analysis
- Reporting using R

To get the most out of this, it is assumed that you have basic knowledge about the command lines, python, R and anaconda. However, explanation of the codes are provided and so following along should be simple.

There are many tutorials on how to install WSL and anaconda and so will not be covered here. 

Found mistakes, have a suggestion or questions, [please submit an issue on GitHub](https://github.com/GiangLeN/Zymo-Mock-sequencing/issues).

It is possible to view this documents locally or online at <https://gianglen.github.io/Zymo-Mock-sequencing/>.

> :warning:
> **Large analysis**: Will consume huge amount of disk space and computational power.

## Contents

### [1. Prepping guides](#1-prepping-guides)

- [Data of interest](#sra-databases-of-zymobiomics)
  * [BGISEQ](#bgiseq)
  
- [Conda environments](#conda-environments)

### [2. Workflow for NCBI's data](#2-workflow-for-downloading-data-from-ncbi)

### [3. Pipeline with Snakemake](02pipeline.md)

- [Pipeline with snakemake](snakemake.md)



## 1. Prepping guides

### SRA databases of zymoBIOMICS
The SRA is a database that contains many samples. We are interested only in the zymoBIOMIC sequencing.
There are many ways to search for this such as ...

The list of SRA files can be found in the `sra_files` directory.

#### BGISEQ
Name, looks at different extraction methods for BGI sequencing.

### Conda environments

We can create an environment with specify name and specific tools/programs as follows:

`conda create --name <environment_name> -c <channel> <tool1> <tool2>`

Create an environment to download SRA files from NCBI.

```
# Create an environment called ncbi_dl from 
# channel bioconda with the program sra-tools
conda create --name ncbi_dl -c bioconda sra-tools

# Activate the environment
conda activate ncbi_dl
```

<details>
  <summary>Environment from yaml</summary>
  
  It is also possible to install conda environment from yaml files in the `envs` folder

```
# Create environment from file
conda env create -f envs/sra-tools.yaml

# Activate the environment
conda activate ncbi_dl

# Deactivate the environment
conda deactivate
```
</details>

## 2. Workflow for downloading data from NCBI

The workflow needs to be understood before setting up the pipeline.
Here, we will download the data from NCBI and then split the file into forward and reverse reads. Finally, these fastq pairs are compressed to save disk space.

Let's try to download the first experimental run from NCBI  
`prefetch ERR4097239`

Once completed, we can see a folder called `ERR4097239` with the file `ERR4097239.sra` inside.

```
# Split file into forward and reverse reads
fasterq-dump --split-files ERR4097239/ERR4097239.sra
```
Congratulation, you have successfully downloaded and split the first sample.

## 3. Pipeline with Snakemake

For multiple files processing, it's much better to use Snakemake.

### Install snakemake

Build the Snakemake's environment using mamba.

```
# Install snakemake via mamba
conda create --name snakemake -c conda-forge mamba

# Activate environment
conda activate snakemake

# Use mamba to install snakemake
# This tutorial use snakemake version 7.8.5
mamba install -c conda-forge -c bioconda snakemake=7.8.5

# Check snakemake version
snakemake --version
```

### Use snakemake to download SRA

Open your favorite editor and create a new file called `Snakefile`.

Creating the first rule to download the next SRA file

```
rule download:
    output:
        "ERR4097239/ERR4097239.sra" 
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        prefetch ERR4097239
        """
```
The output tells Snakemake that we are expecting `ERR4097239/ERR4097239.sra` as the outcome.

You also might have noticed that we are no longer in the `ncbi_dl` environment. This means sra-tools is not installed in here. We could re-install it here. However, multiple tools in the same environment might lead to compatibility conflict. To avoid this, we will create a yaml file contains program information which tie in to the rule.

The file `envs/sra-tools.yaml`

```
name: sra-tools
channels:
 - bioconda
dependencies:
 - sra-tools
```
Under `shell:` is where you type in the commands/scripts.

Above the rule `download` we need to have rule `all` to specify the final outcome.

```
rule all:
    input:
        "ERR4097108/ERR4097108.sra"
```
Check below to see what the snakefile looks like.

<details>
  <summary>Snakefile</summary>

```
rule all:
    input:
        "ERR4097108/ERR4097108.sra"

rule download:
    output:
        "ERR4097108/ERR4097108.sra"
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        prefetch ERR4097108
        """
```
</details>


To run snakemake:

`snakemake --cores 1 --use-conda`


![First snakemake](img/snakemake_st.png?raw=true "Title")

Snakemake will build the environment the first time it runs.

### Extract forward and reverse reads

Add the second rule to split the sra file into forward and reverse files. The output from rule `donwload` is used as the input for the rule `split_raw`.

```
rule split_raw:
    input:
        "ERR4097108/ERR4097108.sra"
    output:
        multiext("ERR4097108/ERR4097108", "_1.fastq", "_2.fastq")
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        fasterq-dump --split-files ERR4097108
        """
```

Update rule `all` to the following:

```
rule all:
    input:
        multiext("ERR4097108.sra", "_1.fastq", "_2.fastq")
```
Snakemake looks at the final output and works backward to identify the input required for each rule.  

<details>
  <summary>Snakefile</summary>

```
rule all:
    input:
        multiext("ERR4097108.sra", "_1.fastq", "_2.fastq")

rule download:
    output:
        "ERR4097108/ERR4097108.sra"
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        prefetch ERR4097108
        """

rule split_raw:
    input:
        "ERR4097108/ERR4097108.sra"
    output:
        multiext("ERR4097108.sra", "_1.fastq", "_2.fastq")
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        fasterq-dump --split-files {input}
        """
```
</details>

### Compress files

Create new rule to compress reads

```
rule compress:
    input:
        multiext("ERR4097108.sra", "_1.fastq", "_2.fastq")
    output:
        multiext("ERR4097108.sra", "_1.fastq.gz", "_2.fastq.gz")
    shell:
        """
        gzip {input}
        """
```

Update rule `all` to the following:

```
rule all:
    input:
        multiext("ERR4097108.sra", "_1.fastq.gz", "_2.fastq.gz")
```

<details>
  <summary>Snakefile</summary>


```
rule all:
    input:
        multiext("ERR4097108.sra", "_1.fastq.gz", "_2.fastq.gz")

rule download:
    output:
        "ERR4097108/ERR4097108.sra"
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        prefetch ERR4097108
        """

rule split_raw:
    input:
        "ERR4097108/ERR4097108.sra"
    output:
        multiext("ERR4097108.sra", "_1.fastq", "_2.fastq")
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        fasterq-dump --split-files {input}
        """

rule compress:
    input:
        multiext("ERR4097108.sra", "_1.fastq", "_2.fastq")
    output:
        multiext("ERR4097108.sra", "_1.fastq.gz", "_2.fastq.gz")
    shell:
        """
        gzip {input}
        """
```
</details>


### Update snakemake to process multiple files


```
SRA = ["ERR4097108", "ERR4097109"]

rule all:
    input:
        expand(["{sra}.sra_1.fastq.gz", "{sra}.sra_2.fastq.gz"], sra = SRA)

rule download:
    output:
        temp("{sra}/{sra}.sra")
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        prefetch {wildcards.sra}
        """

rule split_raw:
    input:
        "{sra}/{sra}.sra"
    output:
        multiext("{sra}.sra", "_1.fastq", "_2.fastq")
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        fasterq-dump --split-files {input}
        """

rule compress:
    input:
        multiext("{sra}.sra", "_1.fastq", "_2.fastq")
    output:
        multiext("{sra}.sra", "_1.fastq.gz", "_2.fastq.gz")
    shell:
        """
        gzip {input}
        """
```

## Repository tree structure

```
Zymo-Mock-sequencing/
├── LICENSE
├── README.md
├── docs
│   └── index.md
└── envs
    └── ncbi_dl.yaml # To download from NCBI

```

- [Heading](#heading)
  * [Sub-heading](#sub-heading)
    + [Sub-sub-heading](#sub-sub-heading)
- [Heading](#heading-1)
  * [Sub-heading](#sub-heading-1)
    + [Sub-sub-heading](#sub-sub-heading-1)

# Heading levels

> This is a fixture to test heading levels

<!-- toc -->

