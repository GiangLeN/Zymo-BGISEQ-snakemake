# Nanopore tutorial: Metagenomic analysis of zymoMock data  

Metagenomic sequencing using nanopore technology has become increasingly popular due to its ability to generate long reads and provide high-resolution sequence data.
In this tutorial, we will use Snakemake to process metagenomic data from nanopore sequencing.

## Snakemake environment

Conda or Mamba can manage your environments and packages by ensuring that all the required dependencies are installed.
This helps to reproduce the workflow on different systems.

We need to create an environment.yaml file that includes *Snakemake* and *pandas*.

```
name: snakemake
channels:
  - bioconda
dependencies:
  - pandas
  - snakemake
```

### Install mamba

Mamba can provide faster performance when working with large and complex environments.

```
conda install -c conda-forge mamba
mamba env create -f environment.yml
conda activate snakemake
```

## Tab-separated values (TSV) input file

TSV files contain rows of data, with each column separated by a tab character. 
The *pandas* library is used to read TSV files as input for the workflow. 
This is the structure of `sample.tsv` file that we will use as input.

| id | fq |
|----|----|
|sample1|/path_to_sample/|
|sample2|/path_to_sample2/|
|sample3|/path_to_sample3/


## Parse TSV as input

```
import pandas

# Read TSV file into a Pandas DataFrame and index it based on a specified column
def parse_samples(samples_tsv):
    return pandas.read_csv(samples_tsv, sep='\t').set_index("id", drop=False)
# Returns the value in the specified column for the row that matches the wildcard variable in the DataFrame. 
def get_sample_id(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample, [col]].dropna()[0]
    
_samples = parse_samples("sample.tsv")

```

## Quality control
The first step in processing metagenomic data from nanopore sequencing is quality control. Quality control is important to ensure that the data generated is of high quality and to filter out any low-quality reads. The quality control process involves removing adapters, filtering out reads with low base-calling quality scores, and trimming reads to remove any low-quality bases.

To perform quality control, there are several tools available such as Porechop, Guppy, and NanoFilt. These tools can remove adapters, filter out low-quality reads, and trim reads to a desired length. It is important to choose the appropriate quality control tool based on the specific needs of your analysis.




## Assembly
Once the quality control step is completed, the next step is to assemble the high-quality reads into contigs. Assembly is the process of merging overlapping reads to create longer sequences. The goal of assembly is to reconstruct the original genome of the organisms present in the sample.

There are several tools available for assembling metagenomic data from nanopore sequencing, such as Canu, Flye, and wtdbg2. These tools can assemble reads into long contiguous sequences.

## Taxonomic classification
The next step in processing metagenomic data from nanopore sequencing is taxonomic classification. Taxonomic classification is the process of identifying the organisms present in the sample.

Several tools are available for taxonomic classification, such as Kraken2, Kaiju, and MetaPhlAn2. These tools use reference databases to identify the organisms present in the sample.

## Functional annotation
The final step in processing metagenomic data from nanopore sequencing is functional annotation. Functional annotation is the process of identifying the functions and pathways encoded in the genomes of the organisms present in the sample.

Several tools are available for functional annotation, such as Prokka, RAST, and MG-RAST. These tools use reference databases to annotate the genes and functions encoded in the genomes.

In conclusion, processing metagenomic data from nanopore sequencing involves several steps, including quality control, assembly, taxonomic classification, and functional annotation. Each of these steps is important to ensure that the data generated is of high quality and can be used to identify the organisms present in the sample and the functions and pathways encoded in their genomes.


## Inspired sites
https://github.com/adlape95/Assembly-methods-nanopore
https://github.com/harvardinformatics/nanopore-assembly/blob/main/Nanopore_quickstart.md
