---
layout: default
title: "Sequencing analysis of the ZymoBIOMICS Mock community"
author: "Ngoc Giang Le"
version: 0.1
#bibliography:
nav_order: 1
description: "Metagenomic tutorial"
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

This project shows how to build pipeline to analyze the metagenomic of publicly available zymoBIOMICS Mock community. This sample contains 8 bacteria and 2 fungi at specific concentration.

The aim of this tutorial is to look at different tools and approaches and compared them for benchmark purposes. Codes, and scripts are also provided to tie various tools together.

There are several topics to be covered in multiple tutorials:

-   Download data from NCBI
-   Environment setup to run programs
-   Build pipelines for analysis
-   Reporting using R

To get the most out of this, it is assumed that you have basic knowledge about the command lines, python, R and anaconda. However, explanation of the codes are provided and so following along should be simple.

Note: There are many tutorials on how to install WSL and `conda` and so will not be covered here.

Found mistakes, have a suggestion or questions, [please submit an issue on GitHub](https://github.com/GiangLeN/Zymo-Mock-sequencing/issues).

It is possible to view this documents locally or online at <https://gianglen.github.io/Zymo-Mock-sequencing/>.

> :warning: **Large analysis**: Will consume huge amount of disk space and computational power.


## 1. Prepping guides

### 1.1 zymoBIOMICS data on NCBI

The Sequence Read Archive (SRA) contains many sequencing databases from different projects.
Similar to their aims, we also hope to be able to reproduce the results and discover new finding from data analysis.
Here we are looking at research projects that sequenced the zymoBIOMIC Mock community.

The list of SRA files that are used for the analysis can be found in the `sra_files` directory.

#### 1.1.1 BGISEQ

Name, looks at different extraction methods for BGI sequencing.
Five different extraction methods were tested.


| SampleID   | SRA        | Protocol |
|------------|------------|----------|
| D6300-27   | ERR4097245 | MetaHIT  |
| D6300-28   | ERR4097111 | MetaHIT  |
| D6300-29   | ERR4097243 | MetaHIT  |
| D6300-30   | ERR4097237 | MetaHIT  |
| D6300-31   | ERR4097238 | MetaHIT  |
| D6300-32   | ERR4097276 | MetaHIT  |
| D6300-2-1  | ERR4097261 | MN       |
| D6300-2-2  | ERR4097241 | MN       |
| D6300-2-3  | ERR4097242 | MN       |
| D6300-3-1  | ERR4097244 | MN       |
| D6300-4-1  | ERR4097272 | MN       |
| D6300-4-2  | ERR4097208 | MN       |
| D6300-13   | ERR4097266 | MP       |
| D6300-14   | ERR4097269 | MP       |
| D6300-15   | ERR4097268 | MP       |
| D6300-16   | ERR4097271 | MP       |
| D6300-17   | ERR4097270 | MP       |
| D6300-18   | ERR4097262 | MP       |
| D6300-11   | ERR4097264 | PS       |
| D6300-12   | ERR4097232 | PS       |
| D6300-2-7  | ERR4097239 | PS       |
| D6300-2-8  | ERR4097240 | PS       |
| D6300-7    | ERR4097176 | PS       |
| D6300-9    | ERR4097177 | PS       |
| D6300-37   | ERR4097212 | Q        |
| D6300-38   | ERR4097211 | Q        |
| D6300-39   | ERR4097210 | Q        |
| D6300-40   | ERR4097172 | Q        |
| D6300-41   | ERR4097173 | Q        |
| D6300-42   | ERR4097174 | Q        |
| D6300-4-20 | ERR4097207 | ZYMO     |
| D6300-4-21 | ERR4097206 | ZYMO     |
| D6300-4-22 | ERR4097205 | ZYMO     |
| D6300-4-23 | ERR4097204 | ZYMO     |
| D6300-4-24 | ERR4097171 | ZYMO     |
| D6300-6-19 | ERR4097175 | ZYMO     |


### 1.2 Conda environments

Conda is great for reprduce... as you can install environment with specific tools/programs.
The basic command:

`conda create --name <environment_name> -c <channel> <tool1> <tool2>`

Create an environment called sra-tools from the channel bioconda with the program sra-tools


`conda create --name sra-tools -c bioconda sra-tools`


Activate the environment

`conda activate sra-tools`

To return to your previous environment, deactivate conda

`conda deactivate`

<details>
  <summary>Conda with yaml</summary>

    It is also possible to install conda environment from yaml files in the `envs` folder
    
    `conda env create -f envs/sra-tools.yaml`
    
</details>


## 2. Pipeline with Snakemake

### 2.1 Workflow for downloading data from NCBI

The workflow needs to be understood before setting up the pipeline. Here, we will download the data from NCBI and then split the file into forward and reverse reads.

Let's try to download the first experimental run from NCBI by activating the environment we created earlier.

```
# Acivate the environment we created earlier
conda activate sra-tools

# Download the SRA
prefetch ERR4097239

```

Once completed, we can see a folder called `ERR4097239` with the file `ERR4097239.sra` inside.
Split file into forward and reverse reads.

`fasterq-dump --split-files ERR4097239/ERR4097239.sra`


Congratulation, you have successfully downloaded and split the first sample.

### 2.2 Setup Snakemake

For multiple files processing, it's much better to use Snakemake.
Reproducible.
Everything is treated the same way.

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

### 2.3 Create rule for download SRA

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

You might have noticed that we are no longer in the `sra-tools` environment. This means sra-tools is not installed here. To avoid compatibility conflict caused by multiple tools in the same environment, we will create a new yaml file. This file contains information of programs for the rule.

Create a file called `sra-tools.yaml` and place it in the `envs` directory:

    name: sra-tools
    channels:
     - bioconda
    dependencies:
     - sra-tools

Under `shell:` is where you type in the commands/scripts. In this case its

`prefetch <SRA_name>`

Above the rule `download` we need to have rule `all` to specify the *final* outcome.

```
rule all:
    input:
        "ERR4097108/ERR4097108.sra"
```

Save your `Snakefile`.

<details>
    <summary>Snakefile</summary>
        
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

</details>




We are telling snakemake to use one core and use conda to create environment for the rules. Snakemake will build the environment the first time it runs.

`snakemake --cores 1 --use-conda`

![First Snakemake](img/snakemake_st.png?raw=true "First time running Snakemake")

### Extract forward and reverse reads {#extract-forward-and-reverse-reads}

Add a second rule to split the SRA file into forward and reverse files. The output from rule `donwload` is used as the input for the rule `split_raw`. We also expect files with multiple extensions. This is specified by the `multiext()`.

    rule split_raw:
        input:
            "ERR4097108/ERR4097108.sra"
        output:
            multiext("ERR4097108.sra", "_1.fastq", "_2.fastq")
        conda:
            "envs/sra-tools.yaml"
        shell:
            """
            fasterq-dump --split-files ERR4097108
            """

Update rule `all` to the following as we want the split files as the final result.

    rule all:
        input:
            multiext("ERR4097108.sra", "_1.fastq", "_2.fastq")

Snakemake looks at the final output and works backward to identify the input required for each rule. This is why you do not need to keep `"ERR4097108/ERR4097108.sra"` in the rule `all`.

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

### Compress files {#compress-files}

Create new rule to compress reads to save space. Once again, the output from `split_raw` is used as the input for the `compress` rule. For the output, we want the compressed files.

    rule compress:
        input:
            multiext("ERR4097108.sra", "_1.fastq", "_2.fastq")
        output:
            multiext("ERR4097108.sra", "_1.fastq.gz", "_2.fastq.gz")
        shell:
            """
            gzip {input}
            """

Update rule `all` to the following:

    rule all:
        input:
            multiext("ERR4097108.sra", "_1.fastq.gz", "_2.fastq.gz")

This is how the Snakefile should look like

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

### Update Snakemake to process multiple files {#update-snakemake-to-process-multiple-files}

To run multiple files, we need to modify the Snakefile with wildcards. The input is a python list, where each value will be processed one by one.

Here we will create a list of two SRA samples.

    SRA = ["ERR4097108", "ERR4097109"]

For rule `all`, we want the final output to be `{sra}.sra_1.fastq.gz` and `{sra}.sra_2.fastq.gz`. The SRA name is replaced with `{sra}`, where its value is extracted from the list called SRA from above. This means, we are looking for 4 files as the final outcome:

-   ERR4097108.sra_1.fastq.gz
-   ERR4097108.sra_2.fastq.gz
-   ERR4097109.sra_1.fastq.gz
-   ERR4097109.sra_2.fastq.gz

The new rule `download`:

    rule download:
        output:
            temp("{sra}/{sra}.sra")
        conda:
            "envs/sra-tools.yaml"
        shell:
            """
            prefetch {wildcards.sra} 
            """

For the `output`. We replace the SRA name with `{sra}`. The `temp()` outside indicates temporary file. Snakemake will remove these files once the pipeline is completed.

We need to specify that we are running the `wildcards.sra` for the `prefetch` command.

Replace the SRA name with `{sra}` for the other rules.

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

Run Snakemake with two cores and `-r` for reason and `-p` for print command.

`snakemake --cores 2 --use-conda -r -p`

You can see that Snakemake is downloading and processing two SRA files at the same time. To avoid stressing NCBI server, reduce the number of cores. Also, multiple cores will consume more computational power and large space of hard disk.

Congratulation, you managed to use Snakemake to download and process multiple files. Try to download other SRA samples for analysis

## Repository tree structure

    Zymo-Mock-sequencing/
    ├── LICENSE
    ├── README.md
    ├── Snakefile
    ├── docs
    │   ├── img
    │   │   └── snakemake_st.png
    │   └── index.md
    └── envs
        └── sra-tools.yaml # To download from NCBI

.condarc -default (first)

conda config --set channel_priority flexible

install packages on a new environment instead of the default base environment. This way you can avoid package conflicts and maintain a clean base environment.
