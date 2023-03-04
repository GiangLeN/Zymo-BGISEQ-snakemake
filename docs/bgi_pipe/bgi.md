---
layout: default
title: "Sequencing analysis of the ZymoBIOMICS Mock community"
author: "Ngoc Giang Le"
version: 0.1.1
date:
#bibliography:
nav_order: 1
description: "Snakemake tutorial to analyse BGI-seq metagenomics"
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

## 1. Introduction

BGI-seq is a sequencing platform developed by the BGI (Beijing Genomics Institute) for whole genome sequencing and metagenomic analysis.

In this tutorial, we will walk through the steps of creating a Snakemake pipeline to process metagenomic data generated from BGI-seq.

## 2. Preprocessing

Create a new directory for this analysis.

```
mkdir bgi_pipeline
cd bgi_pipeline
```

We will use the output `sample.tsv` from previous [tutorial](https://gianglen.github.io/Zymo-Mock-sequencing/snakemake_tut/intro_snakemake.html) for this analysis.

The input file contains *id* and the paths of *read1* and *read2*.

|id          |read1                       |read2                       |
|------------|----------------------------|----------------------------|
| ERR4097245 | path/ERR4097245_1.fastq.gz | path/ERR4097245_2.fastq.gz |
| ERR4097111 | path/ERR4097111_1.fastq.gz | path/ERR4097111_2.fastq.gz |
| ERR4097243 | path/ERR4097243_1.fastq.gz | path/ERR4097243_2.fastq.gz |

Move `sample.tsv` to the newly created directory.

Open a text editor to create a new *Snakefile*.
The following functions parse the file and extract the forward and reverse reads of the sample.

```
import pandas as pd

def parse_samples(samples_tsv):
    # Load in tab separate file
    # Remove samples with incomplete fields
    # Set ID as the index
    return pd.read_csv(samples_tsv, sep ='\t').dropna().set_index("id", drop=False)

def get_files(sample_df, wildcards, col):
    # Return forward and reverse reads based on the sample's name
    return sample_df.loc[wildcards.sample, [col]]

_samples = parse_samples("sample.tsv")

```

## 3. Quality Control

The first step in processing BGI-seq metagenomic data is to perform quality control on the raw reads.
This is typically done using FastQC.
However, there are quite a few samples and we want to automate the process.
We will use fastp as it can detect and remove adapter sequences and trim the low-quality reads.
Information regarding the raw and filtered sample are also provided and can be extracted.
Any reads below 60 bp are removed as recommend by BGI.

```
rule trimming:
    input:
        r1 = lambda wildcards: get_files(_samples, wildcards, "read1"),
        r2 = lambda wildcards: get_files(_samples, wildcards, "read2")
    output:
        r1_trim = "01_trimmed/{sample}_trimmed_1.fq.gz",
        r2_trim = "01_trimmed/{sample}_trimmed_2.fq.gz",
        json = "01_trimmed/{sample}_trimmed.json",
        html = "01_trimmed/{sample}_trimmed.html"
    threads: 3
    params:
        min_length = 60
    conda:
        "envs/fastp.yaml"
    shell:
        """
        fastp --detect_adapter_for_pe -w {threads} -i {input.r1} -I {input.r2} -o {output.r1_trim} -O {output.r2_trim} --n_base_limit 0 --cut_front --cut_tail --length_required {params.min_length} -j {output.json} -h {output.html}
        """

rule readsCheck:
    input:
        "01_trimmed/{sample}_trimmed.json",
    output:
        "reads_status/{sample}_01_trimming.txt"
    shell:
        """
        echo sample:{wildcards.sample} > {output}
        grep -w "before_filtering" {input} -A1 | grep "total_reads" | tr -d '\t", ' >> {output}
        grep "filtering_result" {input} -A5 | sed '1d' | tr -d '\t", ' >> {output}
        """
```

## 4. Host Filtering

Any contamination reads from human must be removed.
This can be done using tools such as Bowtie or BWA to align the reads to a reference genome of the host organism.
Mapped read can then be removed.
This will improve on the assembly down stream and reduce the computational power.




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


