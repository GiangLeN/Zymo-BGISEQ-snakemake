---
layout: default
title: "Sequencing analysis of the ZymoBIOMICS Mock community"
author: "Ngoc Giang Le"
version: 0.1
date:
#bibliography:
nav_order: 1
description: "Snakemake tutorial to analyse metagenomics"
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

This project shows how to build metagenomic pipeline to analyze publicly available data of the zymoBIOMICS Mock community.
The sample contains 8 bacteria and 2 fungi at specific concentration.

The aim of this tutorial is to look at different tools and approaches and compared them for benchmark purposes.
Codes, and scripts are also provided to tie various tools together.

There are several topics to be covered in multiple tutorials:

-   Setup environment to run programs
-   Download data from NCBI
-   Build pipelines for analysis
-   Reporting using R

To get the most out of this, it is assumed that you have basic knowledge about the command lines, Bash, Python, R and anaconda.
However, explanation of the codes are provided and so following along should be simple.

*Note:* There are many tutorials on how to install WSL and `conda` and so will not be covered here.

Found mistakes, have suggestions or questions, [please submit an issue on GitHub](https://github.com/GiangLeN/Zymo-Mock-sequencing/issues).

It is possible to view this documents locally or online at <https://gianglen.github.io/Zymo-Mock-sequencing/>.

> :warning: **Large analysis**: Will consume huge amount of time, disk space and computational power.


## Contents

### [1. Prepping for analysis](01_prepping.html)

### [2. Pipeline with Snakemake](02_workflow.html)




