## index ##

## bwa-mem2 ##

rule bwa_mem2_index:
    input:
        "{results}/{name}/genome/{build}.{release}.genome.fa"
    output:
        "{results}/{name}/index/bwa-mem2/{build}.{release}.0123",
        "{results}/{name}/index/bwa-mem2/{build}.{release}.amb",
        "{results}/{name}/index/bwa-mem2/{build}.{release}.ann",
        "{results}/{name}/index/bwa-mem2/{build}.{release}.bwt.2bit.64",
        "{results}/{name}/index/bwa-mem2/{build}.{release}.pac"
    params:
        out_prefix="{results}/{name}/index/bwa-mem2/{build}.{release}"
    log:
        "{results}/{name}/logs/bwa-mem2.{build}.{release}.log"
    benchmark:
        "{results}/{name}/.benchmark/bwa-mem2.{build}.{release}.tsv"
    conda:
        "../envs/bwa-mem2.yaml"
    resources:
        mem_gb=50
    shell:
        "bwa-mem2 index "
        "-p {params.out_prefix} "
        "{input} "
        ">{log} 2>&1 "


## Bismark ##

rule bismark_input_fastq:
    input:
        "{results}/{name}/genome/{build}.{release}.genome.fa"
    output:
        temp("{results}/{name}/index/bismark/{build}.{release}.genome.fa") # temp
    params:
        in_fasta="../../genome/{build}.{release}.genome.fa",
    shell:
        "[ -f {output} ] && rm {output}; "
        "ln -s {params.in_fasta} {output} "


## Bismark genome preparation
rule bismark_genome_preparation:
    input:
        expand("{{results}}/{{name}}/index/bismark/{build}.{release}.genome.fa",
            build=config["ref"][name]["build"],
            release=config["ref"][name]["gtf"],
        )
    output:
        directory("{results}/{name}/index/bismark/Bisulfite_Genome")
    log:
        "{results}/{name}/logs/bismark_genome_preparation.{name}.log"
    benchmark:
        "{results}/{name}/.benchmark/bismark_genome_preparation.{name}.benchmark.tsv"
    params:
        extra="--verbose"  # optional params string
    conda:
        "../envs/bismark.yaml"
    threads:
        8
    shell:
        "bismark_genome_preparation "
        "--parallel 4 "
        "--bowtie2 "
        "{params.extra} "
        "$(dirname {input}) " # Bismark uses *folder* of fastq as input and writes output to this folder!!!
        ">{log} 2>&1 "

