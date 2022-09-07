import pandas as pd

def parse_samples(samples_tsv):
    return pd.read_csv(samples_tsv, sep ='\t').dropna().set_index("ID", drop=False)

def get_files(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample, [col]]

_samples = parse_samples("samples.tsv")


rule all:
    input:
        "gtdbtk_v2_data.tar.gz",
        expand("05.bins/{sample}/dastool/{sample}_DASTool_summary.tsv", sample = _samples.index)

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

rule trimming:
    input:
        r1 = lambda wildcards: get_files(_samples, wildcards, 'r1'),
        r2 = lambda wildcards: get_files(_samples, wildcards, 'r2'),
    output:
        r1_trim = temp("01.trimmed/{sample}_trimmed_1.fq.gz"),
        r2_trim = temp("01.trimmed/{sample}_trimmed_2.fq.gz"),
        json = "01.trimmed/{sample}_trimmed.json"
    conda:
        "envs/fastp.yaml"
    log:
        "logs/{sample}_trimmed.log"
    benchmark:
        "benchmarks/{sample}_trimmed.time"
    shell:
        """
        fastp --detect_adapter_for_pe  -i {input.r1} -I {input.r2} -o {output.r1_trim} -O {output.r2_trim} -j {output.json} 2> {log}
        """

rule bracken:
    input:
        r1 = "01.trimmed/{sample}_trimmed_1.fq.gz",
        r2 = "01.trimmed/{sample}_trimmed_2.fq.gz",
        kraken_db = "databases/hash.k2d"
    output:
        kraken = temp("02.taxonomy/raw/{sample}.kraken2"),
        bracken = "02.taxonomy/raw/{sample}.bracken"
    threads: 9
    conda:
        "envs/krabraken.yaml"
    params:
        db = "databases",
        report = "02.taxonomy/raw/{sample}_kraken2.report",
    log:
        kraken = "logs/{sample}_kraken.log",
        bracken = "logs/{sample}_bracken.log"
    benchmark:
        "benchmarks/{sample}_bracken.time"
    shell:
        """
        kraken2 --use-names --gzip-compressed --db {params.db} --report {params.report} --confidence 0.1 --threads {threads} {input.r1} {input.r2} > {output.kraken} 2> {log.kraken}
        bracken -d {params.db} -i {params.report} -l S -o {output.bracken} 2> {log.bracken}
        """

rule megahit:
    input:
        r1 = "01.trimmed/{sample}_trimmed_1.fq.gz",
        r2 = "01.trimmed/{sample}_trimmed_2.fq.gz",
    output:
        directory("03.assemblies/{sample}")
    threads: 30
    conda:
        "envs/megahit.yaml"
    log:
        "logs/{sample}_megahitslog"
    benchmark:
        "benchmarks/{sample}_megahit.time"
    shell:
        """
        megahit -1 {input.r1} -2 {input.r2} -o {output} --k-min 25 --k-max 99  --k-step 10 -t {threads} 2> {log}
        """

rule contigs:
    input:
        "03.assemblies/{sample}"
    output:
        "04.contigs/{sample}_contigs.fa"
    shell:
        """
        sed 's/>k.*_/>{wildcards.sample}_/g' {input}/final.contigs.fa > {output}
        """

rule contigIndex:
    input:
        "04.contigs/{sample}_contigs.fa"
    conda:
        "envs/bwa.yaml"
    output:
        temp(multiext("04.contigs/{sample}_contigs.fa",".amb",".ann",".bwt",".pac",".sa"))
    shell:
        """
        bwa index {input} 
        """

rule mapRaws:
    input:
        contigs = "04.contigs/{sample}_contigs.fa",
        read1 = "01.trimmed/{sample}_trimmed_1.fq.gz",
        read2 = "01.trimmed/{sample}_trimmed_2.fq.gz",
        dex = multiext("04.contigs/{sample}_contigs.fa",".amb",".ann",".bwt",".pac",".sa")
    output:
        temp("04.contigs/{sample}.sam")
    conda:
        "envs/bwa.yaml"
    threads: 24
    shell:
        """
        bwa mem -t {threads} {input.contigs} {input.read1} {input.read2} > {output}
        """

rule samBam:
    input:
        "04.contigs/{sample}.sam"
    output:
        bam = "04.contigs/{sample}.bam",
        bai = "04.contigs/{sample}.bam.bai"
    conda:
        "envs/samtools.yaml"
    threads: 24
    shell:
        """
        samtools view -@ {threads} -hb -F4 {input} | samtools sort -o {output.bam}
        samtools index {output.bam}
        """

rule coverage:
    input:
        "04.contigs/{sample}.bam",
    output:
        "04.contigs/{sample}_coverage.txt"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools depth {input} | awk '{{ sum += $3}} END {{print "Average coverage:","\t",sum/NR }}' > {output}
        """ 

rule metabat2:
    input:
        bam = "04.contigs/{sample}.bam",
        bai = "04.contigs/{sample}.bam.bai",
        fa = "04.contigs/{sample}_contigs.fa",
    output:
        list = "05.bins/{sample}/{sample}_metabat_contigs_list",
        depth = "05.bins/{sample}/{sample}_contigs_metabat_depth.txt",
    params:
        bins = "05.bins/{sample}/metabat/{sample}_metabins",
        dir = "05.bins/{sample}/metabat"
    threads: 24
    conda:
        "envs/metabat.yaml"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        metabat2 -i {input.fa} -a {output.depth} -o {params.bins} -t {threads} -v 
        awk '/>/{{sub(">","&"FILENAME"@");sub(/\.fa/,x)}}1'  {params.dir}/*.fa | grep ">" | awk -F"@" '{{print $2,$1}}' | sed 's/ .*\//\t/g' > {output.list}
        """

rule maxbin2:
    input:
        batdone = "05.bins/{sample}/{sample}_metabat_contigs_list",
        fa = "04.contigs/{sample}_contigs.fa",
        depth = "05.bins/{sample}/{sample}_contigs_metabat_depth.txt",
    output:
        list = "05.bins/{sample}/{sample}_maxbin_contigs_list",
        depth = "05.bins/{sample}/{sample}_contigs_maxbin_depth.txt",
    params:
        bins = "05.bins/{sample}/maxbin/{sample}_maxbins",
        dir = "05.bins/{sample}/maxbin"
    threads: 24
    conda:
        "envs/maxbin.yaml"
    shell:
        """
        awk -F"\t" '{{print $1"\t"$4}}' {input.depth} | sed '1d' > {output.depth}
        mkdir -p {params.dir}
        run_MaxBin.pl -thread {threads} -contig {input.fa} -out {params.bins} -abund {output.depth}
        awk '/>/{{sub(">","&"FILENAME"@");sub(/\.fasta/,x)}}1'  {params.dir}/*.fasta | grep ">" | awk -F"@" '{{print $2,$1}}' | sed 's/ .*\//\t/g' > {output.list}
        """

rule das_tools:
    input:
        fa = "04.contigs/{sample}_contigs.fa",
        maxbinBins = "05.bins/{sample}/{sample}_maxbin_contigs_list",
        metabatBins = "05.bins/{sample}/{sample}_metabat_contigs_list",
    output:
        "05.bins/{sample}/dastool/{sample}_DASTool_summary.tsv",
    params:
        "05.bins/{sample}/dastool/{sample}",
    threads: 24
    conda:
        "envs/dastool.yaml"
    shell:
        """
        DAS_Tool -i {input.metabatBins},{input.maxbinBins} -l Metabat,Maxbin -c {input.fa} -o {params} --write_bins -t {threads}
        """
      
rule gtdbtk_db:
    output:
        "gtdbtk_v2_data.tar.gz"
    shell:
        """
        wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_v2_data.tar.gz
        """


