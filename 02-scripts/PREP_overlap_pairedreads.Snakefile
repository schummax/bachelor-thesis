################################################################################
# Project: Training project
# Part: Preparation of the data
# Step: Infer DNA molecule length distribution using fastp
#
# Dependent on:
#   - PREP_remove_hostDNA.Snakefile
#
# Alex Huebner, 18/04/23
################################################################################

import json
import os

import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")
    
os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### SAMPLES ###################################################################
SAMPLES, = glob_wildcards("03-data/processed_data/{sample}_1.fastq.gz")
################################################################################

rule all:
    input:
        "05-results/QUAL_DNA_molecule_length_dist.tsv"

rule merge_overlap:
    output:
        merged = "04-analysis/fastp_overlap/{sample}.merged.fastq.gz",
        pe1 = "04-analysis/fastp_overlap/{sample}_1.unmerged.fastq.gz",
        pe2 = "04-analysis/fastp_overlap/{sample}_2.unmerged.fastq.gz",
        single = "04-analysis/fastp_overlap/{sample}_0.unmerged.fastq.gz",
        json = "04-analysis/fastp_overlap/{sample}.fastp.json"
    message: "Merge overlapping read pairs using Fastp: {wildcards.sample}"
    conda: "ENVS_fastp.yaml"
    resources:
        mem = 32,
        cores = 4 
    params:
        pe1 = "03-data/processed_data/{sample}_1.fastq.gz",
        pe2 = "03-data/processed_data/{sample}_2.fastq.gz"
    threads: 4
    shell:
        """
        fastp --in1 {params.pe1} --in2 {params.pe2} \
            --merge --merged_out {output.merged} \
            --out1 {output.pe1} --out2 {output.pe2} \
            --unpaired1 {output.single} \
            -A -G -l 30 --overlap_len_require 11 \
            --json {output.json} --html /dev/null
        """
        
rule merge_json:
    input:
        expand("04-analysis/fastp_overlap/{sample}.fastp.json", sample=SAMPLES)
    output:
        "05-results/QUAL_DNA_molecule_length_dist.tsv"
    message: "Extract features from JSON files and combine into table"
    run:
        results = []
        for fn in input:
            with open(fn, "rt") as jsonfile:
                fastp = json.load(jsonfile)
                insertsize = fastp['insert_size']
                insertsize['histogram'].append(insertsize['unknown'])
                results.append(pd.DataFrame.from_dict({'length': list(range(len(insertsize['histogram']))),
                                                       'count': insertsize['histogram']}) \
                            .assign(sample=os.path.basename(fn).replace(".fastp.json", "")))

        pd.concat(results)[['sample', 'length', 'count']] \
            .sort_values(['sample', 'length']) \
            .to_csv(output[0], sep="\t", index=False)