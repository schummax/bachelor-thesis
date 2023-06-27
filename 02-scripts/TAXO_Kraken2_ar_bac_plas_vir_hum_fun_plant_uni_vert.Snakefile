################################################################################
# Taxonomic profiling with Kraken2 against a large database consisting of
# archaeal, bacterial, viral, fungal, unicellular, human, and plant genomes (in
# total: 35,367 species) that Maxime downloaded and built for 35-mer in May
# 2022.
#
# Alex Huebner, 23/06/23
################################################################################

import os

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### SAMPLES ###################################################################
SAMPLES, = glob_wildcards("03-data/processed_data/{sample}_1.fastq.gz")
################################################################################

rule all:
    input:
        "05-results/KRK2_ar_bac_plas_vir_hum_fun_plant_uni_vert.tsv",

#### Prepare auxilliary files for assigning taxonomy from ids ##################

rule download_taxdump:
    output:
        "tmp/taxdump.tar.gz"
    message: "Download NCBI taxonomy information"
    params:
        url = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
    shell:
        """
        wget -O {output} {params.url}
        """

rule extract_ncbi_tax_dumps:
    input:
        "tmp/taxdump.tar.gz"
    output:
        "tmp/taxdump/merged.dmp"
    message: "Extract taxdump files required for taxpasta"
    params:
        dir = "tmp/taxdump"
    shell:
        """
        mkdir -p {params.dir}
        for f in names nodes merged; do
            tar -C {params.dir} -xzf {input} ${{f}}.dmp
        done
        """

################################################################################

#### Screen HMP samples (paired-end data) ######################################

rule kraken2_pairedend:
    output:
        "04-analysis/kraken2/{sample}.txt"
    message: "Classify against ar_bac_plas_vir_hum_fun_plant_uni_vert database using Kraken2: {wildcards.sample}"
    conda: "ENVS_Kraken2_Bracken.yaml"
    resources:
        mem = 800,
        cores = 1p6,
    params:
        fastq = lambda wildcards: " ".join([f"03-data/processed_data/{wildcards.sample}_{i}.fastq.gz" for i in [1, 2]]),
        db = "/home/maxime_borry/SDAG_old/02_db/kraken2/ar_bac_plas_vir_hum_fun_plant_uni_vert",
    threads: 16
    shell:
        """
        kraken2 \
            --db {params.db} \
            --threads {threads} \
            --output - \
            --report {output} \
            --confidence 0.15 \
            --use-names \
            --gzip-compressed \
            {params.fastq}
        """

rule summarise_reports_pe:
    input:
        reads = expand("04-analysis/kraken2/{sample}.txt", sample=SAMPLES),
        taxdump = "tmp/taxdump/merged.dmp"
    output:
        "05-results/KRK2_ar_bac_plas_vir_hum_fun_plant_uni_vert.tsv"
    message: "Summarise the Kraken2 reports of the paired HMP samples"
    conda: "ENVS_taxpasta.yaml"
    resources:
        mem = 32,
        cores = 1
    params:
        reports = " ".join([f"04-analysis/kraken2/{sample}.txt" for sample in SAMPLES]),
        taxdir = "tmp/taxdump"
    shell:
        """
        taxpasta merge -p kraken2 -o {output} \
            --output-format TSV \
            --wide \
            --add-name --add-rank \
            --taxonomy {params.taxdir} \
            {params.reports}
        """

################################################################################
