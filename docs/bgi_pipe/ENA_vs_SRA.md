---
layout: default
title: "SRA or ENA for sequencing data?"
author: "Ngoc Giang Le"
version: 0.1.1
date:
#bibliography:
nav_order: 1
description: "Where to download sequencing data?"
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

The NCBI Sequence Read Archive (SRA) and the European Nucleotide Archive (ENA) are two of the largest public repositories of high-throughput sequencing data.
DNA, RNA, and protein sequences from different experiments are stored together with the metadata, which describe the sample information, sequencing methods, and experimental design.
Various tools such as SRA-Tools and Kingfisher can be used to download data from these repositories.
In this tutorial, we will compare these two tools and show you how to use them to download sequencing data from ENA and SRA.


## 1.2 Web browser

To search for a project visit the web browser of SRA or ENA.
For each of the record an unique accession number is given depending on the source database.
SRA sample will start with *SRA* and for European Bioinformatics Institute will start with *ERR*.


For this example we will use the project ***PRJEB38036***

## 1.3 Setup tools

Both Kingfisher and SRA-Tools can be installed using the conda package manager.

### 1.3.1 Conda environments

Through out this tutorial, different programs will be used.
To avoid incompatibility causes by installing of these programs to the main environment, we will use `conda`.
This program can install environment with specific tools/programs, and is great for reproducibility.

To create an environment with specific name use this basic form:

`conda create --name <environment_name> -c <channel> <tool1> <tool2>`

Installing channel for the program can be found on [anaconda web page](https://anaconda.org).

After the environment is installed, activate it using

`conda activate <environment_name>`

To exit and return to your previous environment:

`conda deactivate`


<details>
  <summary>Install conda environment using yaml</summary>

`conda env create -f <yaml_file_path>`
    
***Yaml files for this project are located at the `envs/` directory.***
    
</details>



### 1.3.2 SRA-tools

```
conda install -c bioconda sra-tools
```

To download data using SRA-Tools, you need to use the fasterq-dump command.
Here's an example of how to download multiple sequencing runs:

```
fasterq-dump --split-files <ERR or SRA accessions> --outdir sra_data
```
This command will download the data to a directory called *sra_data*.  

Lets test on one sample:

```
fasterq-dump --split-files ERR4097207 -O sra_data
```

The downloaded file is in fastq format.
 
Check the first 10 lines from ERR4097207_1.fastq. 

```head ERR4097207_1.fastq```


```
@ERR4097207.1 1 length=100
AAGGTCTTGCAGCGATTAAATTGCTGAAAAAAGAAGGCATCGTGACGCTGGGCACCGCCGTCTACAGCGCATCGCAGGGCCTGCTGGCGGCGCTGGCGGG
+ERR4097207.1 1 length=100
A>ECEEFF<DECE;AEEFBDE@85B==EBEEC49AD@DDBE@BE3CDAAE=>EEBEFEDBC5FFBE82EDD?EFD>CD=A?EAEFC=D>=AADC7AAD?B
@ERR4097207.2 2 length=100
ACAACGCGGTTGATTAACACCATACCGGGAATGTTATCCATAAATTGCGCCAGTTCATCGTCACTCAATGCTTTTGAGTGAACAATCAACGCATTACAAC
+ERR4097207.2 2 length=100
FFFFFFFFFFFFFGFFGFFFG>FFFGFFFFFFGFFFFFF@F@GEFGFFFGFFGEFFFFGFFFFFFFEFFGGFFEFFBGFFFGGFGGFFFAFFFGGFFFFG
@ERR4097207.3 3 length=58
CAATACATGTGTCACCTTTGCTACGCCTACGTTGATAACGCAACAATGACGTCCAACA
```
The header of the fastq is replaced with the accession number.


### 1.3.3 Installing Kingfisher

Now we will install and try out kingfisher.

```
conda install -c bioconda kingfisher
```

Create directory to download the same ENA file.

```
mkdir ena_data
cd ena_data
kingfisher get -r ERR4097207 -m ena-ftp
```

Lets check the header of the file.

```
zcat ERR4097207_1.fastq.gz | head
```

```
@ERR4097207.1 CL100062948L2C001R001_20/1
AAGGTCTTGCAGCGATTAAATTGCTGAAAAAAGAAGGCATCGTGACGCTGGGCACCGCCGTCTACAGCGCATCGCAGGGCCTGCTGGCGGCGCTGGCGGG
+
A>ECEEFF<DECE;AEEFBDE@85B==EBEEC49AD@DDBE@BE3CDAAE=>EEBEFEDBC5FFBE82EDD?EFD>CD=A?EAEFC=D>=AADC7AAD?B
@ERR4097207.2 CL100062948L2C001R001_44/1
ACAACGCGGTTGATTAACACCATACCGGGAATGTTATCCATAAATTGCGCCAGTTCATCGTCACTCAATGCTTTTGAGTGAACAATCAACGCATTACAAC
+
FFFFFFFFFFFFFGFFGFFFG>FFFGFFFFFFGFFFFFF@F@GEFGFFFGFFGEFFFFGFFFFFFFEFFGGFFEFFBGFFFGGFGGFFFAFFFGGFFFFG
@ERR4097207.3 CL100062948L2C001R001_46/1
CAATACATGTGTCACCTTTGCTACGCCTACGTTGATAACGCAACAATGACGTCCAACA
```

The original header is retained along with the accession.


## 1.4 Conclusion

Both Kingfisher and SRA-Tools are easy to use and capable of downloading data from ENA and SRA.
However, there are some differences to consider when choosing which tool to use.
For us kingfisher can download files in parallel and so much faster than SRA-tools.
It also retains the original header and can be useful for some downstream analysis.



[Back to main page](https://gianglen.github.io/Zymo-Mock-sequencing/)
