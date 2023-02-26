---
layout: default
title: "Intorduction to Snakemake"
author: "Ngoc Giang Le"
version: 0.1.1
date:
#bibliography:
nav_order: 1
description: "Snakemake tutorial to download data"
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

## 1.1 Introduction

Snakemake is a workflow management system used in bioinformatics and data science.
It allows users to write complex workflow in a simple way using syntax which specifies input, output, and commands to be executed.
This makes it easy to write, organize, visualize and reproduce the workflow.
Snakemake supports multiple programming languages.
With its built-in support for parallelization and cluster computing, Snakemake makes it easy to scale up the analysis and/or use in the high-performance computing environments.

## 1.2 Setup Snakemake environments

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

## 1.3 Snakemake rules

Open your favorite editor and create a new file called `Snakefile`.

Structure of a simplified Snakemake rule is as follow:

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

### 1.3.1 Download rule

We will create a rule called `download` and use it to get ERR4097207 files.

```
rule download:
    output:
        "ERR4097207_1.fastq.gz",
        "ERR4097207_2.fastq.gz"
    conda:
        "envs/kingfisher.yaml"
    shell:
        """
        kingfisher get -r ERR4097207 -m ena-ftp
        """
```
[Comparison between Kingfisher and SRA-tools.](https://gianglen.github.io/Zymo-Mock-sequencing/bgi_pipe/ENA_vs_SRA.html)  
[Alternative way to download from ENA.](#alternative-method-to-download-ena)  


We also need to include a rule to specify the *final* outcome.
It is possible to list the two files similar to the input above.
However, we will use the function ***multiext("Common_name", "extension1", "extension2")***.

```
rule all:
    input:
        multiext("ERR4097207", "_1.fastq.gz", "_2.fastq.gz")
```

Save and exit out of the editor.

Lets break down the rule.  

-   `input:` No input file and so was not specified.
-   `output:` File ***ERR4097207_1.fastq.gz*** and ***ERR4097207_2.fastq.gz*** expected as the outcome.
-   `shell:` Running commands/scripts. 
-   `conda:` Conda environment to run the rule

You might have noticed that we in the `snakemake` environment and *kingfisher* is not installed here.
One option is to install *kingfisher* on to the same environment.
However, we want to maintain a clean base environment and avoid program conflict causes by different packages.
To solve this, we will tell Snakemake to build a new conda environment for the rule instead.  
The yaml file contains tools and packages information.  


Create *envs/kingfisher.yaml* from scratch.
You can use other text editor instead of Vim.

```
mkdir -p envs
vim envs/kingfisher.yaml 
```

Paste the text below and save the file.

```
name: kingfisher
channels:
 - bioconda
 - conda-forge
dependencies:
 - kingfisher
```

<details>
  <summary>Snakefile</summary>
  
```
rule all:
    input:
        multiext("ERR4097207", "_1.fastq.gz", "_2.fastq.gz")
          
rule download:
    output:
        "ERR4097207_1.fastq.gz",
        "ERR4097207_2.fastq.gz"
    conda:
        "envs/kingfisher.yaml"
    shell:
        """
        kingfisher get -r ERR4097207 -m ena-ftp
        """
```
            
</details>  


Start the pipeline with

```
snakemake --cores 1 --use-conda
```

We are running Snakemake locally, and so only use one core.

> Congratulation, you have successfully in using Snakemake to downloaded the first sample.


### 1.3.2 Download multiple files

Snakemake is very efficient in processing multiple files. 
We need to to rewrite rules using wildcards.
The SRA of interest is provided as a python list.

```
SRA = ["ERR4097111", "ERR4097109", "ERR4097207"]
```

In all the rules, the SRA's name is replaced with `{sample}` variable.

We will also upgrade our rule for additional settings

```
rule download:
    output:
        multiext("{sample}", "_1.fastq.gz", "_2.fastq.gz")
    conda:
        "envs/kingfisher.yaml"
    threads: 3
    conda:
        "envs/kingfisher.yaml"
    shell:
        """
        kingfisher get -r {wildcards.sample} -f fastq.gz -m ena-ftp --check-md5sums --download-threads {threads}
        """
```

-   `threads:` Number of threads for the rule. Default = 1
-   `shell:` To use the variables in the command, we need to specify that we are running `{wildcards.sample}`.  



For rule `all`, we want the final output to be `{sample}_1.fastq.gz` and `{sample}_2.fastq.gz`.
Here we need to specify that the value of `{sample}` comes from the SRA list from above.

`expand(["{sample}_1.fastq.gz", "{sample}_2.fastq.gz"], sample = SRA)`

This means, we are looking for 6 files as the final outcome:

-   ERR4097111_1.fastq.gz
-   ERR4097111_2.fastq.gz
-   ERR4097109_1.fastq.gz
-   ERR4097109_2.fastq.gz
-   ERR4097207_1.fastq.gz
-   ERR4097207_2.fastq.gz


<details>
  <summary>Snakefile</summary>

```
SRA = ["ERR4097111", "ERR4097109", "ERR4097207"]

rule all:
    input:
        expand(["{sample}_1.fastq.gz", "{sample}_2.fastq.gz"], sample = SRA)

rule download:
    output:
        multiext("{sample}", "_1.fastq.gz", "_2.fastq.gz")
    conda:
        "envs/kingfisher.yaml"
    threads: 3
    conda:
        "envs/kingfisher.yaml"
    shell:
        """
        kingfisher get -r {wildcards.sample} -f fastq.gz -m ena-ftp --check-md5sums --download-threads {threads}
        """
```

</details>  

Run Snakemake with 6 cores.
Use the command `-r` (reason) to see why the rule was triggered and `-p` (print) to show running commands.

```
snakemake --cores 6 --use-conda -r -p
```

Snakemake detects the existent of the ERR4097207 files and only download the missing SRA instead.
The number of parallel process depends on the number of cores provided.
In our case two download rules can be ran at the same time.
To avoid stressing out the NCBI server, reduce the number of input *cores* or increase the number of *threads*.
Also, multiple cores will consume more computational power and huge amount of hard disk.


### 1.3.3 Download from input file

We will update the *Snakemake* file further to download SRA from a file with accessions.
This is the content of the file `accessions.txt`:

```
accession
ERR4097268
ERR4097269
ERR4097270
ERR4097271
ERR4097272
```

First we need python packages to handle the input file.
The following functions parse the file and extract the accessions.

```
import pandas as pd

def parse_samples(accession_list):
    # create pandas Dataframe.
    # Remove empty entries
    # set accession as index and do not remove the original accession column
    return pd.read_csv(accession_list, sep ='\t').dropna().set_index("accession", drop=False)

_samples = parse_samples("accession.txt")

```

<details>
  <summary>Snakefile</summary>

```
import pandas as pd

def parse_samples(accession_list):
    return pd.read_csv(accession_list, sep ='\t').dropna().set_index("accession", drop=False)

_samples = parse_samples("accession.txt")

rule all:
    input:
        expand(["{sample}_1.fastq.gz", "{sample}_2.fastq.gz"], sample = _samples.index)

rule download:
    output:
        multiext("{sample}", "_1.fastq.gz", "_2.fastq.gz")
    conda:
        "envs/kingfisher.yaml"
    threads: 3
    conda:
        "envs/kingfisher.yaml"
    shell:
        """
        kingfisher get -r {wildcards.sample} -f fastq.gz -m ena-ftp --check-md5sums --download-threads {threads}
        """
```

</details>  


### 1.3.4 Create input file for analysis

Our next aim is to create a tab separated file called `sample.tsv` with the following structure for later analysis.
This file contains *id* and the path for *forward* and *reverse* reads of the samples we want to process.

|id          |read1                       |read2                       |
|------------|----------------------------|----------------------------|
| ERR4097245 | path/ERR4097245_1.fastq.gz | path/ERR4097245_2.fastq.gz |
| ERR4097111 | path/ERR4097111_1.fastq.gz | path/ERR4097111_2.fastq.gz |
| ERR4097243 | path/ERR4097243_1.fastq.gz | path/ERR4097243_2.fastq.gz |

Lets create a rule to extract this information.

```
rule pipe_samples:
    input:
        multiext("{sample}", "_1.fastq.gz", "_2.fastq.gz")
    output:
        temp("00_raws/{sample}.txt")
    shell:
        """
        echo {wildcards.sample} $(realpath {input}) | sed 's/ /\t/g' > {output}
        """
```

-   `input`: The output from previous rule is used as the input
-   `output`: Write information to *00_raws/{sample}.txt*. Remove once the file is no longer used.
-   `shell`: Print the SRA name and the paths. Replace space with tab

Update the rule *all* as follow:

```
rule all:
    input:
#        expand(["00_raws/{sample}_1.fq.gz", "00_raws/{sample}_2.fq.gz"], sample = _samples.index),
        expand("00_raws/{sample}.txt", sample = _samples.index)
```

Comment out the line will stop Snakemake from creating that output.
In our case, since we already downloaded the files it would have no impact.

Run Snakemake with

```
snakemake --cores 6 --use-conda -r -p
```

### 1.3.5 Combine results

We need to group all the text files together into one file called `sample.tsv`.

```
rule final_samples:
    input:
        expand("00_raws/{sample}.txt", sample = _samples.index)
    output:
        "sample.tsv"
    shell:
        """
        cat {input} > {output}
        """
```

The `expand` function in the input allows to combine all *00_raws/{sample},txt* together.  

Once again we update the rule *all*.

```
rule all:
    input:
        "sample.tsv"
```

*Note:* It is possible to remove other lines from rule `all` and does not impact the pipeline at all.
This is because, Snakemake looks at the final output and works backward to identify the input required for each rule.
In the *final_samples* rule, the *sample* is already defined and so Snakemake can still trigger the required rules.

Run Snakemake

```
snakemake --cores 6 --use-conda -r -p
```
> Congratulation, you used Snakemake to download and process multiple files.

Try too download other files.
[Here](https://gianglen.github.io/Zymo-Mock-sequencing/accessions/PRJEB38036_mocks.txt) is the full accessions list from the same project. 

## Alternative download method from ENA

Coming soon
