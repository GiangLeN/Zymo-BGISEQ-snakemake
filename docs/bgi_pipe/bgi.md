---
layout: default
title: "Sequencing analysis of the ZymoBIOMICS Mock community"
author: "Ngoc Giang Le"
version: 0.1.1
date:
#bibliography:
nav_order: 1
description: "Snakemake tutorial to analyse metagenomics"
permalink: /
output:
    html_document:
      toc: true
      toc_float: true
      toc_depth: 4
      fig_caption: true
      highlight: tango
---


```{=html}
<style>
body {
text-align: justify}
</style>
```

## Introduction

BGI-seq is a sequencing platform developed by the BGI (Beijing Genomics Institute) for whole genome sequencing and metagenomic analysis. In this tutorial, we will walk through the steps of processing metagenomic data generated from BGI-seq using a combination of software tools and pipelines.

Preprocessing
Quality Control
The first step in processing BGI-seq metagenomic data is to perform quality control on the raw reads. This is typically done using software tools such as FastQC, which generates quality reports and identifies potential issues with the sequencing data.

bash
Copy code
fastqc -o output_directory input_file.fastq.gz
Trimming
After performing quality control, we may want to trim the raw reads to remove low-quality regions and adapter sequences. This can be done using software tools such as Trimmomatic or Cutadapt.

bash
Copy code
trimmomatic PE -phred33 input_file_R1.fastq.gz input_file_R2.fastq.gz output_file_R1.fastq.gz output_file_R2.fastq.gz ILLUMINACLIP:adapter_file.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
Host Filtering
If we are interested in analyzing the microbial content of a sample, we may want to remove any reads that originate from the host organism. This can be done using tools such as Bowtie or BWA to align the reads to a reference genome of the host organism.

bash
Copy code
bowtie2-build host_genome.fa host_genome_index
bowtie2 -x host_genome_index -1 input_file_R1.fastq.gz -2 input_file_R2.fastq.gz -S host_aligned.sam
samtools view -bS host_aligned.sam | samtools sort -o host_aligned.bam
samtools index host_aligned.bam
samtools view -bS -f 12 -F 256 host_aligned.bam > non_host.bam
Assembly
Co-assembly
After preprocessing, we can perform assembly of the metagenomic data using software tools such as MEGAHIT or IDBA-UD. Co-assembly of multiple samples can be performed to improve the quality and completeness of the resulting contigs.

bash
Copy code
megahit -1 input_file_R1.fastq.gz -2 input_file_R2.fastq.gz -o output_directory
Gene Prediction
After assembly, we can predict the genes present in the metagenomic data using tools such as Prodigal or MetaGeneMark. This can provide information on the functional content of the metagenome.

bash
Copy code
prodigal -i contigs.fa -o genes.faa -a proteins.fna -p meta
Taxonomic Classification
After predicting genes, we may want to assign taxonomic classifications to the genes and contigs using software tools such as Kraken or MetaPhlAn. This can provide insight into the taxonomic composition of the metagenomic data.

bash
Copy code
kraken2 --db database --threads 4 --output output_file --report report_file genes.faa
Functional Annotation
After assigning taxonomic classifications, we can annotate the functions of the predicted genes using software tools such as eggNOG-mapper or KEGG. This can provide insight into the functional content of the metagenomic data.

bash
Copy code
emapper.py -i proteins.faa -o output_directory --cpu 4 --data_dir data_directory
Conclusion
In this tutorial, we have walked through the steps of processing BGI-seq metagenomic data, including quality control, trimming, host filtering, assembly, gene prediction, taxonomic classification

There are several topics to be covered in multiple tutorials:

-   Setup environment to run programs
-   Download data from NCBI
-   Metagenomic workflow
-   Build metagenomic pipelines using Snakemake
-   Reporting using R


## 1. Preparing sequencing data

We are interested in sequencing data of zymoBIOMICS from different sequencer and project.







Create a new folder for the project.

```
mkdir bgi_download

cd bgi_download
```

To download SRA from ncbi server, we need `sra-tools` from the *bioconda* channel.

There are multiple ways to approach a problem.
Similarly, there are many ways to download SRA samples

The workflow/approach needs to be understood before setting up the pipeline with Snakemake.  

We can create a new environment as so.

`conda create --name sra-tools -c bioconda sra-tools`

Activate the environment to run `sra-tools`.

`conda activate sra-tools`

Let's download the experimental run ERR4097245 from NCBI.

```
prefetch ERR4097245

```

Once completed, we can see `ERR4097245/ERR4097245.sra` as the downloaded file.
Extract forward and reverse reads from the *sra* file.

```
fasterq-dump --split-files ERR4097245/ERR4097245.sra
```
Two newly generated files ***ERR4097245.sra_1.fastq*** and ***ERR4097245.sra_2.fastq*** correspond to forward and reverse, respectively.
Compress these files with:

```
gzip ERR4097245.sra_1.fastq ERR4097245.sra_2.fastq
```

### 1.4 Download workflow with Snakemake

Now that we know about the workflow, we will use Snakemake to help us to download and process SRA files.
Snakemake is a workflow management system for reproducible and scalable data analyses.

#### 1.4.1 Snakemake environments

*Mamba* is the preferred way to install Snakemake as it's solver is much faster than *conda*.

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

#### 1.4.2 Snakemake rules

Open your favorite editor and create a new file called `Snakefile`.

A simple Snakemake rule is as follow:

```
rule <rule_name>:
    input:
        "input_file"
    output:
        "output_file" 
    shell:
        """
        running commands
        """

```

#### 1.4.3 Download rule

We will create a rule called `download` and use it to get ERR4097245 file.


```
rule download:
    output:
        "ERR4097245/ERR4097245.sra" 
    shell:
        """
        prefetch -f yes ERR4097245
        """
```

Lets break down the rule above.  

-   `input:` No input file and so was not specified.
-   `output:` File ***ERR4097245/ERR4097245.sra*** expected as the outcome (same as workflow).
-   `shell:` Running commands/scripts. The `-f yes` tells *sra-tools* to force download the file.


We also need to include a rule to specify the *final* outcome.

```
rule all:
    input:
        "ERR4097245/ERR4097245.sra"
```

Save and exit out of the editor.


<details>
  <summary>Snakefile</summary>
  
```
rule all:
    input:
        "ERR4097245/ERR4097245.sra"
          
rule download:
    output:
        "ERR4097245/ERR4097245.sra" 
    shell:
        """
        prefetch -f yes ERR4097245
        """
```
            
</details>


Start the pipeline with

```
snakemake --cores 1
```

We are running snakemake locally, and so only use one core.

```
Building DAG of jobs...
Nothing to be done (all requested files are present and up to date).
```

Since we just downloaded ERR4097245 as an example earlier, Snakemake recognizes the existence of the final files and so did not run.
Note: If you want to re download ERR4097245 using Snakemake remove `ERR4097245/ERR4097245.sra`.


Let's change the `Snakefile` to download the next SRA sample ERR4097111.


```
rule all:
    input:
        "ERR4097111/ERR4097111.sra"
            
rule download:
    output:
        "ERR4097111/ERR4097111.sra" 
    shell:
        """
        prefetch -f yes ERR4097111
        """
```

![Missing prefetch](img/snakemake_error1.png?raw=true "First error")

Snakemake gives an error as it does not recognize the command *prefetch*.

You might have noticed that we are no longer in the `sra-tools` environment and *sra-tools* is not installed here.
One option is to install *sra-tool* on to the `snakemake` environment.
However, we want to maintain a clean base environment and avoid program conflict causes by different packages.
To solve this, we will tell Snakemake to build a new conda environment for the rule instead.


The yaml file contains information of the programs.
You can use the pre-made `envs/sra-tools.yaml` or create one from scratch.

```
name: sra-tools
channels:
 - bioconda
dependencies:
 - sra-tools=2.11.0
```

Update rule `download` with the line `conda: "envs/sra-tools.yaml"`.


```
rule download:
    output:
        "ERR4097111/ERR4097111.sra" 
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        prefetch -f yes ERR4097111
        """
```

Save your `Snakefile` and include `--use-conda` when execute.

```
snakemake --cores 1 --use-conda
```

Snakemake will download and build the conda environment the first time it runs.  

![Snakemake conda](img/snakemake_st.png?raw=true "Build conda env for first run")

> Congratulation, you have successfully in using Snakemake to downloaded the first sample.


#### 1.4.4 Extract forward and reverse reads

Same as the second step of the workflow, we want to split the `sra` file into forward and reverse files.
Create a new rule called `split_raw`, which use `fasterq-dump`.


```
rule split_raw:
    input:
        "ERR4097111/ERR4097111.sra"
    output:
        multiext("ERR4097111", "_1.fastq", "_2.fastq")
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        fasterq-dump --split-files ERR4097111
        """
```

-   `input:` Output from rule `download`. 
-   `multiext()` tells the program that we are expecting files with multiple extensions.


Update rule `all` as we want new output as the final result.


<details>
  <summary>Snakefile</summary>
  
```
rule all:
    input:
        "ERR4097111/ERR4097111.sra",
        multiext("ERR4097111", "_1.fastq", "_2.fastq")

rule download:
    output:
        "ERR4097111/ERR4097111.sra"
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        prefetch -f yes ERR4097111
        """

rule split_raw:
    input:
        "ERR4097111/ERR4097111.sra"
    output:
        multiext("ERR4097111", "_1.fastq", "_2.fastq")
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        fasterq-dump --split-files {input}
        """
```
            
</details>


*Note:* It is possible to remove line `"ERR4097111/ERR4097111.sra"` from rule `all` and does not impact the pipeline at all.
Snakemake looks at the final output and works backward to identify the input required for each rule.
As ***ERR4097111/ERR4097111.sra*** is the input of `split_raw`, Snakemake will trigger rule `download` to create this input.

#### 1.4.5 Compress files

Use these points to make a rule called `compress`.

- `input:` Use the output from `split_raw`.
- `output:` file ends with *.gz*
- `shell:` Compress with `gzip <file1> <file2>`
- Update rule `all` with the final outcome.

Check here to see the final `Snakefile` should look like.

<details>
  <summary>Snakefile</summary>
  
```
rule all:
    input:
        multiext("ERR4097111", "_1.fastq.gz", "_2.fastq.gz")

rule download:
    output:
        "ERR4097111/ERR4097111.sra"
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        prefetch -f yes ERR4097111
        """

rule split_raw:
    input:
        "ERR4097111/ERR4097111.sra"
    output:
        multiext("ERR4097111", "_1.fastq", "_2.fastq")
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        fasterq-dump --split-files {input}
        """

rule compress:
    input:
        multiext("ERR4097111", "_1.fastq", "_2.fastq")
    output:
        multiext("ERR4097111", "_1.fastq.gz", "_2.fastq.gz")
    shell:
        """
        gzip {input}
        """
       
```
</details>

Run Snakemake to get the compressed files.


#### 1.4.6 Update Snakemake to process multiple files

Snakemake is very efficient in processing multiple files. 
We need to to rewrite rules using wildcards.
The SRA of interest is provided as a python list.

```
SRA = ["ERR4097111", "ERR4097109"]
```

In all the rules, the SRA's name is replaced with `{sra}` variable.

```
rule download:
    output:
        "{sra}/{sra}.sra"
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        prefetch -f yes {wildcards.sra} 
        """
```

To use the variables in the command, we need to specify that we are running `wildcards.sra`.

For rule `all`, we want the final output to be `{sra}_1.fastq.gz` and `{sra}_2.fastq.gz`.
Here we need to specify that the value of `{sra}` comes from the SRA list from above.

`expand(["{sra}_1.fastq.gz", "{sra}_2.fastq.gz"], sra = SRA)`

This means, we are looking for 4 files as the final outcome:

-   ERR4097111_1.fastq.gz
-   ERR4097111_2.fastq.gz
-   ERR4097109_1.fastq.gz
-   ERR4097109_2.fastq.gz


<details>
  <summary>Snakefile</summary>

```
SRA = ["ERR4097111", "ERR4097109"]

rule all:
    input:
        expand(["{sra}_1.fastq.gz", "{sra}_2.fastq.gz"], sra = SRA)


rule download:
    output:
        temp("{sra}/{sra}.sra")
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        prefetch -f yes {wildcards.sra} 
        """

rule split_raw:
    input:
        "{sra}/{sra}.sra"
    output:
        multiext("{sra}", "_1.fastq", "_2.fastq")
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        fasterq-dump --split-files {input}
        """

rule compress:
    input:
        multiext("{sra}", "_1.fastq", "_2.fastq")
    output:
        multiext("{sra}", "_1.fastq.gz", "_2.fastq.gz")
    shell:
        """
        gzip {input}
        """
```

</details>

Run Snakemake with two cores.
Use the command `-r` to see why the rule was triggered and `-p` to show running commands.

```
snakemake --cores 2 --use-conda -r -p
```

Snakemake downloads and processes two SRA files at the same time.
Number of parallel process depends on the number of cores used.
To avoid stressing out the NCBI server, reduce the number of cores.
Also, multiple cores will consume more computational power and huge amount of hard disk.

> Congratulation, you used Snakemake to download and process multiple files.

Exercise: Try to download the other SRAs for analysis



## 2. Metagenomic workflow




## 3 Taxonomic typing of raw data

We want to know what is present in the raw reads.

```
mkdir bgi_pipeline

cd bgi_pipeline
```

What is the different between kraken2 and metaphlan3?
  [Short tutorial comparing Kraken2 vs metaplan3](mini_taxo/taxo_compare.md).

Let's create a new project to analyze the downloaded and processed SRA files.

```
# Create a new folder
mkdir bgi_pipeline

# Move processed SRA files to the new folder
mv ERR*.fastq.gz SRA_pipeline

# Navigate to the new directory
cd bgi_pipeline

```

Create a tab separated sample file called `sample.tsv`.
This file contains *ID* and the paths for *forward* and *reverse* reads of the samples we want to process.

*How to create this file*

|ID        |r1                         |r2                         |
|----------|---------------------------|---------------------------|
| D6300-27 | ERR4097245.sra_1.fastq.gz | ERR4097245.sra_2.fastq.gz |
| D6300-28 | ERR4097111.sra_1.fastq.gz | ERR4097111.sra_2.fastq.gz |
| D6300-29 | ERR4097243.sra_1.fastq.gz | ERR4097243.sra_2.fastq.gz |
| D6300-30 | ERR4097237.sra_1.fastq.gz | ERR4097237.sra_2.fastq.gz |
| D6300-31 | ERR4097238.sra_1.fastq.gz | ERR4097238.sra_2.fastq.gz |
| D6300-32 | ERR4097276.sra_1.fastq.gz | ERR4097276.sra_2.fastq.gz |

Create a brand new `Snakefile` in the `SRA_pipeline` directory.

Create rule to run bracken on the raw reads.
First we need python packages to handle the `sample.tsv` file.
The following functions parse the file and extract the forward and reverse reads of the sample.


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





Use the following rules to download kraken2 database from the index zone (citation).

> Warning: This database is over 60 GB of space

```
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

The rules above download the kraken/bracken database to the folder called `databases` and extract it. Once done the zipped *tar* file is removed. 




```
rule all:
  input:
  expand("01.taxonomy/raw/{sample}.bracken", sample = _samples.index)

rule bracken:
  input:
  # extract the forward (r1) and reverse (r2) reads
  r1 = lambda wildcards: get_files(_samples, wildcards, 'r1'),
r2 = lambda wildcards: get_files(_samples, wildcards, 'r2'),
kraken_db = "database/hash.k2d"
output:
  kraken = temp("01.taxonomy/raw/{sample}.kraken2"),
bracken = "01.taxonomy/raw/{sample}.bracken"
threads: 80
conda:
  "envs/krabraken.yaml"
params:
  db = "database",
report = "01.taxonomy/raw/{sample}_kraken2.report",
shell:
  """
        kraken2 --use-names --gzip-compressed --db {params.db} --report {params.report} --confidence 0.1 --threads {threads} {input.r1} {input.r2} > {output.kraken}
        bracken -d {params.db} -i {params.report} -l S -o {output.bracken}
```


This group of tutorial aims to show how to use snakemake for reproducible research. 

Accuracy, time and speed


fastq to check for raw reads

Assembly SPADES vs MEGAHIT


Build the pipeline over time



Use the wrapper `temp()` to tell Snakemake to remove the intermediate output.
The file is removed at the end of the run.


