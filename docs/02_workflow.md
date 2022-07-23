## 2. Pipeline with Snakemake

### 2.1 Workflow for downloading data from NCBI

The workflow needs to be understood before setting up the pipeline.

Here, we will download the data from NCBI and then split the file into forward and reverse reads.

Let's try to download the first experimental run from NCBI by activating the environment we created earlier.

```
# Acivate the environment we created earlier
conda activate sra-tools

# Download the SRA
prefetch ERR4097239

```

Once completed, we can see a folder called `ERR4097239` with the file `ERR4097239.sra` inside.

To extract forward and reverse reads from the *sra* file


```
fasterq-dump --split-files ERR4097239/ERR4097239.sra
```

Two files ***ERR4097239.sra_1.fastq*** and ***ERR4097239.sra_2.fastq*** correspond to forward and reverse, respectively.

Congratulation, you have successfully downloaded and split the first sample.

### 2.2 Setup Snakemake

For multiple files processing, it's much better to use Snakemake.

It is easier to reproduce and everything is treated the same way and so will have less errors.


Use *mamba* to build the Snakemake environment.

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

### 2.3 Create rule to download SRA

Open your favorite editor and create a new file called `Snakefile`.
Creating the first rule to download the next SRA.


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

The output tells Snakemake that we are expecting ***ERR4097239/ERR4097239.sra*** as the outcome.  


You might have noticed that we are no longer in the `sra-tools` environment.
This means sra-tools is not installed here.
We will tell Snakemake to create conda environment for the rule instead of install packages on the default base environment.
This was done to avoid package conflict causes by multiple tools and maintain a clean base environment
The yaml file contains information of the programs.

You can use the pre made `envs/sra-tools.yaml` or create one from scratch:

```
name: sra-tools
channels:
 - bioconda
dependencies:
 - sra-tools
```

Under `shell:` is where you type in the commands/scripts. In this case its for downloading,

` prefetch <SRA_name> `

Above the rule `download` we need to have rule `all` to specify the *final* outcome.

```
rule all:
    input:
        "ERR4097108/ERR4097108.sra"
```

Save your `Snakefile` and it should look like this:

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

To run Snakemake

`snakemake --cores 1 --use-conda`

We are running snakemake locally so we use one core. Snakemake will build the environment the first time it runs.  

![Snakemake screen](img/snakemake_st.png?raw=true "First time running Snakemake")

### 2.4 Extract forward and reverse reads

Add a second rule to split the SRA file into forward and reverse files.


```
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
```

We use the output from the `download` rule as the input for the rule `split_raw`. 

The function `multiext()` tells the program that we are expecting files with multiple extensions.

Update rule `all` as we want the split files as the final result.


```
rule all:
    input:
        "ERR4097108/ERR4097108.sra",
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

It is possible to remove line `"ERR4097108/ERR4097108.sra"` from rule `all` and does not impact the pipeline at all.
Snakemake looks at the final output and works backward to identify the input required for each rule.
***ERR4097108/ERR4097108.sra*** is the input of `split_raw` and so Snakemake will activate rule `download` before running `split_raw`.

### 2.5 Compress files

Create new rule to compress newly split reads to save space.

Follow these points to create `compress` rule.

- Use `split_raw`'s output as `compress`'s input

- Compress file as output (file type .gz)

- The command for compressing files `gzip <file1> <file2>`

- Update rule `all` to the final outcome.


The Snakefile should look like this.

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

### 2.6 Update Snakemake to process multiple files

To run multiple files, we need to modify the Snakefile with wildcards.
We specify the interested SRAs in a python list.

```
SRA = ["ERR4097108", "ERR4097109"]
```

In the rule, the name of SRA is replaced with `{sra}` variable.

```
rule download:
    output:
        temp("{sra}/{sra}.sra") # The `temp()` indicates temporary file.
    conda:
        "envs/sra-tools.yaml"
    shell:
        """
        prefetch {wildcards.sra} 
        """
```

Snakemake will remove temporary files at the end of the rule.


For rule `all`, we want the final output to be `{sra}.sra_1.fastq.gz` and `{sra}.sra_2.fastq.gz`.

Here we need to specify that the value of `{sra}` comes from the `SRA` list from above.

`expand(["{sra}.sra_1.fastq.gz", "{sra}.sra_2.fastq.gz"], sra = SRA)`

This means, we are looking for 4 files as the final outcome:

-   ERR4097108.sra_1.fastq.gz
-   ERR4097108.sra_2.fastq.gz
-   ERR4097109.sra_1.fastq.gz
-   ERR4097109.sra_2.fastq.gz


To use the variables in the command, we need to specify that we are running `wildcards.sra`.


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

Run Snakemake with a two cores.
Use the command `-r` to see why the rule was triggered and `-p` to show running commands.

`snakemake --cores 2 --use-conda -r -p`

Snakemake downloads and processes two SRA files at the same time.
Number of parallel process depends on the number of cores used.
To avoid stressing out the NCBI server, reduce the number of cores.
Also, multiple cores will consume more computational power and huge amount of hard disk.

> Congratulation, you managed to use Snakemake to download and process multiple files.
> Try to download other SRA for analysis

Taxonomic tutorial >>> 

