## GRC ##

## genome ##

rule GRC_get_genome:
    output:
        "{results}/{name}/genome/{build}.{release}.genome.fa"
    params:
        species=lambda wildcards: config["ref"][wildcards.name]["species"],
        datatype="dna",
        build=lambda wildcards: config["ref"][wildcards.name]["build"],
        release=lambda wildcards: config["ref"][wildcards.name]["release"],
    log:
        "{results}/{name}/logs/get_genome.{build}.{release}.log"
    wildcard_constraints:
        build="GRC.+"
    cache:
        True  # save space and time with between workflow caching (see docs)
    resources:
        api_requests=1
    wrapper:
        "v1.4.0/bio/reference/ensembl-sequence"


rule GRC_genome_cdna:
    output:
        "{results}/{name}/genome/{build}.{release}.cdna.fa"
    params:
        species=lambda wildcards: config["ref"][wildcards.name]["species"],
        datatype="cdna",
        build=lambda wildcards: config["ref"][wildcards.name]["build"],
        release=lambda wildcards: config["ref"][wildcards.name]["release"],
    log:
        "{results}/{name}/logs/genome_cdna.{build}.{release}.log"
    cache:
        True  # save space and time with between workflow caching (see docs)
    resources:
        api_requests=1
    wrapper:
        "v1.4.0/bio/reference/ensembl-sequence"


## annotation ##

rule GRC_get_annotation:
    output:
        temp("{results}/{name}/annotation/{build}.{release}.gtf"),
    params:
        species=lambda wildcards: config["ref"][wildcards.name]["species"],
        release=lambda wildcards: config["ref"][wildcards.name]["gtf"],
        build=lambda wildcards: config["ref"][wildcards.name]["build"],
        flavor="",  # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
    log:
        "{results}/{name}/logs/get_annotation.{build}.{release}.log",
    wildcard_constraints:
        build="GRC.+"
    cache: True  # save space and time with between workflow caching (see docs)
    resources:
        api_requests=1
    wrapper:
        "v1.4.0/bio/reference/ensembl-annotation"


## index ##

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

