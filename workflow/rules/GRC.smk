## GRC ##

## genome ##

rule get_genome:
    output:
        "{results}/{name}/genome/{build}.{release}.genome.fa"
    params:
        species=lambda wildcards: config["ref"][wildcards.name]["species"],
        datatype="dna",
        build=lambda wildcards: config["ref"][wildcards.name]["build"],
        release=lambda wildcards: config["ref"][wildcards.name]["release"],
    log:
        "{results}/{name}/logs/get_genome.{build}.{release}.log"
    cache:
        True  # save space and time with between workflow caching (see docs)
    resources:
        api_requests=1
    wrapper:
        "v1.4.0/bio/reference/ensembl-sequence"


rule genome_index:
    input:
        "{results}/{name}/genome/{build}.{release}.genome.fa"
    output:
        "{results}/{name}/genome/{build}.{release}.genome.fa.fai"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input} "


rule genome_gatk_dict:
    input:
        "{results}/{name}/genome/{build}.{release}.genome.fa"
    output:
        "{results}/{name}/genome/{build}.{release}.genome.dict"
    log:
        "{results}/{name}/logs/genome_gatk_dict.{build}.{release}.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk CreateSequenceDictionary "
        "-R {input} "
        ">{log} 2>&1 "


rule genome_cdna:
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

rule get_annotation:
    output:
        temp("{results}/{name}/annotation/{build}.{release}.gtf"),
    params:
        species=lambda wildcards: config["ref"][wildcards.name]["species"],
        release=lambda wildcards: config["ref"][wildcards.name]["gtf"],
        build=lambda wildcards: config["ref"][wildcards.name]["build"],
        flavor="",  # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
    log:
        "{results}/{name}/logs/get_annotation.{build}.{release}.log",
    cache: True  # save space and time with between workflow caching (see docs)
    resources:
        api_requests=1
    wrapper:
        "v1.4.0/bio/reference/ensembl-annotation"


rule bgzip_gtf:
    input:
        "{results}/{name}/annotation/{build}.{release}.gtf"
    output:
        gtf="{results}/{name}/annotation/{build}.{release}.gtf.gz",
        tbi="{results}/{name}/annotation/{build}.{release}.gtf.gz.tbi"
    conda:
        "../envs/tabix.yaml"
    shell:
        "( "
        "grep '^#' {input}; "
        "grep -v '^#' {input} "
        "| sort -k1,1 -k4,4n "
        ") "
        "| bgzip "
        ">{output.gtf} "
        "&& tabix -p gff {output.gtf} "


rule gtf2bed:
    input:
        "{results}/{name}/annotation/{build}.{release}.gtf.gz"
    output:
        "{results}/{name}/annotation/{build}.{release}.bed"
    conda:
        "../envs/ucsc.yaml"
    shell:
        "zcat {input} "
        """| awk '{{if ($3!="gene") print $0}}' """
        "| grep -v '^#' "
        "| gtfToGenePred stdin stdout "
        ##"| gtfToGenePred -genePredExt stdin stdout "
        ##"| gtfToGenePred -geneNameAsName2 -genePredExt stdin stdout "
        ##"""| awk '{{symbol=$12; if (symbol!="none") $1=$1"-"symbol; print $0}}' """
        "| genePredToBed stdin {output} "


rule gtf2genes_bed:
    input:
        "{results}/{name}/annotation/{build}.{release}.gtf.gz"
    output:
        "{results}/{name}/annotation/{build}.{release}.genes.bed"
    shell:
        "zcat {input} "
        """| awk '{{if ($3=="gene") print $0}}' """
        """| awk 'BEGIN{{OFS="\t"}}{{ """
        """match($0,/gene_id "([a-zA-Z0-9_-]+)"/,gene_id); """
        """match($0,/gene_name "([a-zA-Z0-9_-]+)"/,symbol); """
        """name=gene_id[1]; """
        """if (symbol[1]!="") name=name"-"symbol[1]; """
        """print $1, $4-1, $5, name, 0, $7}}' """
        ">{output} "


rule bed2genes_saf:
    input:
        "{results}/{name}/annotation/{build}.{release}.genes.bed"
    output:
        "{results}/{name}/annotation/{build}.{release}.genes.saf"
    shell:
        "( "
        "echo -e 'GeneID\tChr\tStart\tEnd\tStrand'; "
        "cat {input} "
        """| awk 'BEGIN{{OFS="\t"}}{{print $4, $1, $2, $3, $6}}' """
        ") "
        ">{output} "


rule gtf2id_tsv:
    input:
        "{results}/{name}/annotation/{build}.{release}.gtf.gz"
    output:
        "{results}/{name}/annotation/{build}.{release}.ids.tsv"
    conda:
        "../envs/ucsc.yaml"
    shell:
        "zcat {input} "
        """| awk '{{if ($3!="gene") print $0}}' """
        """| awk 'BEGIN{{OFS="\t"}}{{ """
        """match($0,/transcript_id "([a-zA-Z0-9_-]+)"/,transcript_id); """
        """match($0,/gene_id "([a-zA-Z0-9_-]+)"/,gene_id); """
        """match($0,/gene_name "([a-zA-Z0-9_-]+)"/,symbol); """
        """name=gene_id[1]; """
        """print transcript_id[1], gene_id[1], symbol[1]}}' """
        "| sed '/^\s*$/d' " # delete "empty" lines
        "| sort -u "
        ">{output} "


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

