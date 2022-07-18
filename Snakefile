SRA = ["ERR4097276"]

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
