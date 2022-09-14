import pandas as pd
import os

def parse_samples(samples_tsv):
    return pd.read_csv(samples_tsv, sep ='\t').dropna().set_index("ID", drop=False)

def get_files(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample, [col]]

_samples = parse_samples("samples.tsv")



rule all:
    input:
        #expand("05.bins/{sample}/dastool/{sample}_DASTool_summary.tsv", sample = _samples.index)
        expand("05.bins/{sample}/{sample}_bin_taxonomy.txt", sample = _samples.index),
        "zymo_gbk_prepared",
        expand("02.taxonomy/raw/{sample}.bracken", sample = _samples.index),
        expand("{sample}_bin.txt", sample = _samples.index),

      
rule gtdbtk_download:
    output:
        temp("databases/gtdbtk_v2_data.tar.gz")
    shell:
        """
        wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_v2_data.tar.gz
        """

rule gtdbtk_database:
    input:
        "databases/gtdbtk_v2_data.tar.gz"
    output:
        directory("databases/release207_v2")
    shell:
        """
        tar zxvf databases/gtdbtk_v2_data.tar.gz -C databases
        """

rule cgview_dl:
    output:
        directory("cgview_comparison_tool")
    conda:
        "envs/git.yaml"
    shell:
        """
        git clone https://github.com/paulstothard/cgview_comparison_tool.git
        """

rule kraken_db:
    output:
        temp("databases/k2_standard_20220607.tar.gz")
    params:
        directory("databases")
    shell:
        """
        wget -P {params} https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20220607.tar.gz
        """

checkpoint zymoMocks:
    output:
        gnome = temp(directory("databases/BioPool_genomes")),
        zip = temp("ZymoBIOMICS.STD.genomes.ZR160406.zip")
    conda:
        "envs/zip.yaml"
    shell:
        """
        wget https://s3.amazonaws.com/zymo-files/BioPool/ZymoBIOMICS.STD.genomes.ZR160406.zip
        unzip -d databases ZymoBIOMICS.STD.genomes.ZR160406.zip
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
        r1_trim = "01.trimmed/{sample}_trimmed_1.fq.gz",
        r2_trim = "01.trimmed/{sample}_trimmed_2.fq.gz",
        json = "01.trimmed/{sample}_trimmed.json",
        html = "01.trimmed/{sample}_trimmed.html"
    conda:
        "envs/fastp.yaml"
    log:
        "logs/{sample}_trimmed.log"
    benchmark:
        "benchmarks/{sample}_01_trimmed.time"
    shell:
        """
        fastp --detect_adapter_for_pe  -i {input.r1} -I {input.r2} -o {output.r1_trim} -O {output.r2_trim} -j {output.json} -h {output.html} 2> {log}
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
        "benchmarks/{sample}_02_bracken.time"
    shell:
        """
        kraken2 --use-names --gzip-compressed --db {params.db} --report {params.report} --confidence 0.1 --threads {threads} {input.r1} {input.r2} > {output.kraken} 2> {log.kraken}
        bracken -d {params.db} -i {params.report} -l S -o {output.bracken} 2> {log.bracken}
        """

rule megahit:
    input:
        r1 = "01.trimmed/{sample}_trimmed_1.fq.gz",
        r2 = "01.trimmed/{sample}_trimmed_2.fq.gz",
        bracken = "02.taxonomy/raw/{sample}.bracken"
    output:
        directory("03.assemblies/{sample}")
    threads: 30
    conda:
        "envs/megahit.yaml"
    log:
        "logs/{sample}_megahit.log"
    benchmark:
        "benchmarks/{sample}_03_megahit.time"
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

rule mapRaws:
    input:
        contigs = "04.contigs/{sample}_contigs.fa",
        read1 = "01.trimmed/{sample}_trimmed_1.fq.gz",
        read2 = "01.trimmed/{sample}_trimmed_2.fq.gz",
    output:
        index = temp(multiext("04.contigs/{sample}_contigs.fa",".amb",".ann",".bwt",".pac",".sa")),
        sam = temp("04.contigs/{sample}.sam"),
        bam = "04.contigs/{sample}.bam",
        bai = "04.contigs/{sample}.bam.bai"
    conda:
        "envs/mapRaws.yaml"
    benchmark:
        "benchmarks/{sample}_04_rawMap.time"
    shell:
        """
        bwa index {input.contigs} 
        bwa mem -t {threads} {input.contigs} {input.read1} {input.read2} > {output.sam}
        samtools view -@ {threads} -hb -F4 {output.sam} | samtools sort -o {output.bam}
        samtools index {output.bam}
        """

rule coverage:
    input:
        "04.contigs/{sample}.bam",
    output:
        "04.contigs/{sample}_coverage.txt"
    conda:
        "envs/samtools.yaml"
    benchmark:
        "benchmarks/{sample}_05_coverage.time"
    shell:
        """
        samtools depth {input} | awk '{{ sum += $3}} END {{print "Average coverage:","\t",sum/NR }}' > {output}
        """ 

rule metabat2:
    input:
        # keep old result
        bam = ancient("04.contigs/{sample}.bam"),
        bai = ancient("04.contigs/{sample}.bam.bai"),
        fa = ancient("04.contigs/{sample}_contigs.fa"),
    output:
        list = "05.bins/{sample}/{sample}_metabat_contigs_list",
        depth = "05.bins/{sample}/{sample}_contigs_metabat_depth.txt",
    params:
        bins = "05.bins/{sample}/metabat/metabat_bin",
        dir = "05.bins/{sample}/metabat"
    threads: 24
    conda:
        "envs/metabat.yaml"
    benchmark:
        "benchmarks/{sample}_06_metabat.time"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        metabat2 -i {input.fa} -a {output.depth} -o {params.bins} -t {threads} -v 
        grep ">" {params.dir}/*.fa | sed 's/.*\///g;s/\.fa.*>/ /g' | awk '{{print $2"\t"$1}}' > {output.list}
        """

rule maxbin2:
    input:
        fa = "04.contigs/{sample}_contigs.fa",
        depth = "05.bins/{sample}/{sample}_contigs_metabat_depth.txt",
    output:
        list = "05.bins/{sample}/{sample}_maxbin_contigs_list",
        depth = temp("05.bins/{sample}/{sample}_contigs_maxbin_depth.txt"),
    params:
        bins = "05.bins/{sample}/maxbin/maxbin",
        dir = "05.bins/{sample}/maxbin"
    threads: 24
    conda:
        "envs/maxbin.yaml"
    benchmark:
        "benchmarks/{sample}_07_maxbin.time"
    shell:
        """
        awk -F"\t" '{{print $1"\t"$4}}' {input.depth} | sed '1d' > {output.depth}
        mkdir -p {params.dir}
        run_MaxBin.pl -thread {threads} -contig {input.fa} -out {params.bins} -abund {output.depth}
        grep ">" {params.dir}/*.fasta | sed 's/.*\///g;s/\.fa.*>/ /g' | awk '{{print $2"\t"$1}}' > {output.list}
        """

checkpoint das_tools:
    input:
        fa = "04.contigs/{sample}_contigs.fa",
        maxbinBins = "05.bins/{sample}/{sample}_maxbin_contigs_list",
        metabatBins = "05.bins/{sample}/{sample}_metabat_contigs_list",
    output:
        directory("05.bins/{sample}/dastool/{sample}_DASTool_bins")
    params:
        "05.bins/{sample}/dastool/{sample}",
    threads: 24
    benchmark:
        "benchmarks/{sample}_08_dastool.time"
    conda:
        "envs/dastool.yaml"
    shell:
        """
        DAS_Tool -i {input.metabatBins},{input.maxbinBins} -l Metabat,Maxbin -c {input.fa} -o {params} --write_bins -t {threads}
        """

rule checkM:
    input:
        "05.bins/{sample}/dastool/{sample}_DASTool_bins"
    output:
        das = "05.bins/{sample}/{sample}_checkM.txt",
        tmpdas = temp(directory("05.bins/{sample}/tmp_dastool"))
    params:
        checkmFolder = "05.bins/{sample}/checkM",
        checkmMarkers = "05.bins/{sample}/checkM/{sample}_checkM.markers"
    threads: 24
    benchmark:
        "benchmarks/{sample}_09_checkm.time"
    conda:
        "envs/checkm.yaml"
    shell:
        """
        mkdir -p {output.tmpdas}
        checkm tree -t {threads} -x fa {input} {params.checkmFolder}
        checkm lineage_set {params.checkmFolder} {params.checkmMarkers}
        checkm analyze -t {threads} -x fa {params.checkmMarkers} {input} {params.checkmFolder} --tmpdir {output.tmpdas}
        checkm qa -t {threads} {params.checkmMarkers} {params.checkmFolder} -f {output.das}
        """
        
rule gtdbtk:
    input:
        das = "05.bins/{sample}/dastool/{sample}_DASTool_bins",
        db = "databases/release207_v2"
    output:
        "05.bins/{sample}/{sample}_bin_taxonomy.txt"
    conda:
        "envs/gtdbtk.yaml"
    params:
        directory("05.bins/{sample}/taxa_gtdbtk")
    benchmark:
        "benchmarks/{sample}_10_gtdbtk.time"
    threads: 4
    shell:
        """
        cd {input.db}
        export "GTDBTK_DATA_PATH=$(pwd)" 
        cd -
        gtdbtk classify_wf --genome_dir {input.das} --out_dir {params} --extension fa --cpus {threads}
        awk -F"\t" '{{print $1,$2}}' {params}/gtdbtk.bac120.summary.tsv > {output}
        """

rule zymo_gbk:
    input:
        "databases/BioPool_genomes/genomes/{zymo}.fasta"
    output:
        gbk = "databases/zymoMocks/{zymo}/zymo_{zymo}.gbk",
        others = temp(multiext("databases/zymoMocks/{zymo}/zymo_{zymo}", ".err", ".faa", ".ffn", ".fna", ".fsa", ".gff", ".log", ".sqn", ".tbl", ".tsv", ".txt"))
    conda:
        "envs/prokka.yaml"
    params:
        "databases/zymoMocks/{zymo}"
    threads: 2
    shell:
        """
        genus=$(echo {wildcards.zymo} | sed 's/_.*//g')
        species=$(echo {wildcards.zymo} | sed 's/.*_//g')
        prokka {input} --outdir {params} --prefix "zymo_"{wildcards.zymo} --cpus {threads} --force --genus $genus --species $species
        """

def get_zymo(wildcards):
    zymodir = checkpoints.zymoMocks.get(**wildcards).output["gnome"] 
    zymos, = glob_wildcards(os.path.join(zymodir,"genomes/{zymos}.fasta")) 
    return [ "databases/zymoMocks/" + zymo + "/zymo_" + zymo + ".gbk" for zymo in zymos ]

rule zymo_done:
    input:
        get_zymo
    output:
        touch("zymo_gbk_prepared")

rule bin_prokka:
    input:
        "05.bins/{sample}/dastool/{sample}_DASTool_bins/{bin}.fa",
    output:
        gbk = "06.compareBins/{sample}/{bin}/{bin}.gbk",
        dir = directory("06.compareBins/{sample}/{bin}"),
    conda:
        "envs/prokka.yaml"
    benchmark:
        "benchmarks/{sample}_11_prokka_{bin}.time"
    threads: 2
    shell:
        """
        prokka {input} --outdir {output.dir} --prefix {wildcards.bin} --cpus {threads} --force
        """

rule cgView:
     input:
        cgview = "cgview_comparison_tool",
        zymo = "zymo_gbk_prepared",
        bingbk = "06.compareBins/{sample}/{bin}/{bin}.gbk",
        bintaxo = "05.bins/{sample}/{sample}_bin_taxonomy.txt"
     output:
        "06.compareBins/{sample}/{sample}_{bin}_dvd_med.png"
     conda:
        "envs/cgview.yaml"
     params:
        module = "lib/perl_modules",
        taxa = "06.compareBins/{sample}/{bin}.tax",
        cgv = "06.compareBins/{sample}/{bin}_cgv",
        zymo = "databases/zymoMocks/"
     benchmark:
        "benchmarks/{sample}_12_cgview_{bin}.time"
     shell:
        """
        # paths for cgview
        cd {input.cgview}
        export CCT_HOME=$(pwd)  
        cd {params.module}
        export PERL5LIB=$(pwd)  
        cd ../../../
        
        if [ -f "${{CCT_HOME}}"/cog_db/.complete ]; then
            echo "COG BLAST database already created"
        else
            mkdir -p "${{CCT_HOME}}"/cog_db
            cp "${{CCT_HOME}}"/lib/scripts/assign_cogs/db/whog.gz "${{CCT_HOME}}"/cog_db
            cp "${{CCT_HOME}}"/lib/scripts/assign_cogs/db/myva.gz "${{CCT_HOME}}"/cog_db

            gunzip "${{CCT_HOME}}"/cog_db/whog.gz
            gunzip "${{CCT_HOME}}"/cog_db/myva.gz

            echo "Preparing COG BLAST database"
            formatdb -p T -i "${{CCT_HOME}}"/cog_db/myva -o T -l "${{CCT_HOME}}"/cog_db/formatdb.log

            echo "COG BLAST database created"
            touch "${{CCT_HOME}}"/cog_db/.complete
        fi

        # Identify bin taxonomy
        grep $(basename {input.bingbk} ".gbk") {input.bintaxo} | sed 's/.*;s__//g;s/Limosilactobacillus/Lactobacillus/g' > {params.taxa}

        # Locate reference file
        species=$(cat {params.taxa} | sed 's/ /_/g')
        species_found=$(find {params.zymo} -name "zymo_${{species}}.gbk")

        genus=$(echo $species | sed 's/_.*//g')
        genus_found=$(find {params.zymo} -name "zymo_${{genus}}*.gbk")

        cgview_comparison_tool/scripts/build_blast_atlas.sh -i {input.bingbk} -p {params.cgv}

        if [[ -f $species_found ]]; then
            cp $species_found {params.cgv}"/comparison_genomes/"
        elif [[ -f $genus_found ]]; then
            cp $genus_found {params.cgv}"/comparison_genomes/"
        fi

        cgview_comparison_tool/scripts/build_blast_atlas.sh -p {params.cgv} -z medium
        cp {params.cgv}"/maps_for_dna_vs_dna/dna_vs_dna_medium.png"  {output}
        """

def get_bins(wildcards):
    tmpdir = checkpoints.das_tools.get(**wildcards).output[0] 
    return expand("06.compareBins/{sample}/{sample}_{bin}_dvd_med.png",
           sample = wildcards.sample,
           bin = glob_wildcards(os.path.join(tmpdir,"{bin}.fa")).bin) 

rule cgview_done:
     input:
        get_bins
     output:
        # Must specify sample here for snakemake to know the wildcard
        "{sample}_bin.txt"
     shell:
        """
        echo {input} > {output} 
        """

