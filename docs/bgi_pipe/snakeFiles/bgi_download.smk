import pandas as pd

def parse_samples(accession_list):
    return pd.read_csv(accession_list, sep ='\t').dropna().set_index("accession", dro
p=False)

_samples = parse_samples("accession.txt")

rule all:
    input:
#        expand(["00_raws/{sample}_1.fq.gz", "00_raws/{sample}_2.fq.gz"], sample = _samples.index)
#        expand("00_raws/{sample}.txt", sample = _samples.index)
        "sample.tsv"


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

rule pipe_samples:
    input:
        multiext("{sample}", "_1.fastq.gz", "_2.fastq.gz")
    output:
        temp("00_raws/{sample}.txt")
    shell:
        """
        echo {wildcards.sample} $(realpath {input}) | sed 's/ /\t/g' > {output}
        """

rule final_samples:
    input:
        expand("00_raws/{sample}.txt", sample = _samples.index)
    output:
        "sample.tsv"
    shell:
        """
        echo -e "id\tread1\tread2" > {output}
        cat {input} >> {output}
        """
