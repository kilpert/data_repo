rule T2T_get_genome:
    output:
        "{results}/{name}/genome/{build}.{release}.genome.fa"
    params:
        url=lambda wildcards: config["ref"][wildcards.name]["genome"]["fasta_url"],
        download="{results}/{name}/genome/download.fa.gz"
    log:
        "{results}/{name}/logs/get_genome.{build}.{release}.log"
    wildcard_constraints:
        build="T2T.+"
    cache:
        True  # save space and time with between workflow caching (see docs)
    resources:
        api_requests=1
    shell:
        "wget "
        ##"-nv " # non-verbose basic output
        "--progress=dot:mega "
        "-O {params.download} "
        "{params.url} "
        ">{log} 2>&1; "
        "zcat {params.download} "
        ">{output}; "
        "rm {params.download} "


rule T2T_get_annotation_gtf:
    output:
        temp("{results}/{name}/annotation/{build}.{release}.gtf")
    params:
        url=lambda wildcards: config["ref"][wildcards.name]["annotation"]["gtf_url"],
        download="{results}/{name}/annotation/annotation.gtf.gz"
    log:
        "{results}/{name}/logs/get_annotation_gtf.{build}.{release}.log"
    wildcard_constraints:
        build="T2T.+"
    conda:
        "../envs/tabix.yaml"
    cache:
        True  # save space and time with between workflow caching (see docs)
    resources:
        api_requests=1
    shell:
        "wget "
        "--progress=dot:mega "
        "-O {params.download} "
        "{params.url} "
        ">{log} 2>&1; "
        "zcat {params.download} "
        ">{output}; "
        "rm {params.download} "


# rule T2T_gtf2genes_bed:
#     input:
#         "{results}/{name}/annotation/{build}.{release}.gtf.gz"
#     output:
#         "{results}/{name}/annotation/{build}.{release}.genes.bed"
#     wildcard_constraints:
#         build="T2T.+"
#     shell:
#         "zcat {input} "
#         """| awk '{{if ($3=="gene") print $0}}' """
#         """| awk 'BEGIN{{OFS="\t"}}{{ """
#         """match($0,/gene_id "([a-zA-Z0-9_-]+)"/,gene_id); """
#         """match($0,/gene_name "([a-zA-Z0-9_-]+)"/,symbol); """
#         """name=gene_id[1]; """
#         """if (symbol[1]!="") name=name"-"symbol[1]; """
#         """print $1, $4-1, $5, name, 0, $7}}' """
#         ">{output} "

