---
layout: default
title: "Sequencing analysis of the ZymoBIOMICS Mock community"
author: "Ngoc Giang Le"
version: 0.1.1
date:
#bibliography:
nav_order: 1
description: "Snakemake tutorial to analyse metagenomics sequencing data"
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


## ZymoBIOMICS Mock Community

The [ZymoBIOMICS Mock Community](https://zymoresearch.eu/collections/zymobiomics-microbial-community-standards) is a synthetic microbial community consists of 8 bacterial and 2 yeast strains that are commonly found in environmental and human microbiomes.
The concentration of these species are known.
ZymoBIOMICS Mock Community is used as a positive control in microbiome studies to assess the accuracy, precision, and sensitivity of microbiome sequencing methods.
It can also be used to compare results across different studies and sequencing platforms.

## Public repositories of nucleotide sequences

The [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) and the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home) are both public repositories of nucleotide sequence data.
Both repositories store a wide range of data such as DNA, RNA, and protein sequences from different organisms and sample types.
Using publicly available data, we will download and build metagenomic pipelines.

## Aims

The tutorials will look at different tools and approaches and compared them for benchmark purposes.
Codes, and scripts are also provided to tie various tools together.
Our aim is to use these data to understand more about the sequencing method and optimize the pipeline.
We hope to reproduce the original results as well as discover new findings.

## Preparation

To get the most out of this, it is assumed that you have basic knowledge about the command lines, Bash, Python, R and anaconda.
However, explanation of the codes are provided and so following along should be simple.

*Note:* There are many tutorials on how to install WSL and `conda` and so will not be covered here.

Found mistakes, have suggestions or questions, [please submit an issue on GitHub](https://github.com/GiangLeN/Zymo-Mock-sequencing/issues).

It is possible to view this documents locally or online at <https://gianglen.github.io/Zymo-Mock-sequencing/>.

> :warning: **Large analysis**: Will consume huge amount of time, disk space and computational power.

## Projects with zymoMock

| Project  | Sequencer | Samples | Files | Ref |
|----------|-----------|---------|-------|-----|
|PRJEB38036|BGI-SEQ    | 25      |       |     |


## Metagenomic pipelines tutorials

- [BGI-SEQ]()
- [Nanopore]()

## Small tutorial/tests

- [ENA vs SRA](bgi_pipe/ENA_vs_SRA.html)
- [Introduction to Snakemake]()