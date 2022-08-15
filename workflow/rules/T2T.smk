rule get_genome_T2T:
    output:
        "{results}/{name}/genome/{build}.{release}.genome.fa"
    params:
        url="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz",
    log:
        "{results}/{name}/logs/get_genome.{build}.{release}.log"
    cache:
        True  # save space and time with between workflow caching (see docs)
    resources:
        api_requests=1
    shell:
        "wget "
        "-O {output} "
        "{params.url} "

