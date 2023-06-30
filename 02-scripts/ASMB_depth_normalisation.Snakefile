################################################################################
# Normalisation of the depth of the k-mers using bbNorm
#
# Due to the sheer size of information, we are not able to efficiently assemble
# the ice core sequencing data. Therefore, I will try to normalise the
# sequencing depth using bbnorm first.
#
# Alex Huebner, 23/06/23
################################################################################

from math import ceil
import os

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### SAMPLES ###################################################################
SAMPLES, = glob_wildcards("03-data/processed_data/{sample}_1.fastq.gz")
################################################################################

#### Auxilliary functions ######################################################

def return_memory_requirement(wildcards):
    return int(open(checkpoints.determine_memory.get(**wildcards).output[0], "rt") \
               .readline().rstrip())

################################################################################

localrules: determine_memory

rule all:
    input:
        expand("tmp/bbnorm/{sample}_{i}.bbnorm.fastq.gz", sample=SAMPLES, i=[1, 2])

#### BBnorm ####################################################################

rule loglog:
    input:
        sample = ["03-data/processed_data/{sample}_1.fastq.gz", "03-data/processed_data/{sample}_2.fastq.gz"]
    output:
        touch("tmp/bbnorm/{sample}.loglog_done")
    message: "Estimate the cardinality of the k-mer space: {wildcards.sample}"
    resources:
        mem = 120,
        mem_gb = 120,
        cores = 1
    threads: 1
    log: "tmp/bbnorm/{sample}.loglog"
    benchmark: "tmp/bbnorm/benchmark.loglog.{sample}.txt"
    wrapper:
        "v2.0.0/bio/bbtools/loglog"

checkpoint determine_memory:
    input:
        "tmp/bbnorm/{sample}.loglog_done"
    output:
        "tmp/bbnorm/{sample}.memory_requirement"
    message: "Calculate the memory requirement for bbnorm: {wildcards.sample}"
    params:
        log = "tmp/bbnorm/{sample}.loglog"
    run:
        with open(params.log, "rt") as infile:
            for line in infile:
                if line.startswith("Cardinality"):
                    cardinality = int(re.search(r'Cardinality: +([0-9]+)',
                                                line.rstrip()).group(1))
                    memory = ceil(12 * cardinality / 1e9)
        with open(output[0], "wt") as outfile:
            outfile.write(f"{memory}\n")

rule bbnorm:
    input:
        fastq = ["03-data/processed_data/{sample}_1.fastq.gz", "03-data/processed_data/{sample}_2.fastq.gz"],
        memory = "tmp/bbnorm/{sample}.memory_requirement"
    output:
        fastq = ["tmp/bbnorm/{sample}_1.bbnorm.fastq.gz", "tmp/bbnorm/{sample}_2.bbnorm.fastq.gz"]
    message: "Apply BBnorm to sequencing data for down-sampling and error-correction: {wildcards.sample}"
    resources:
        mem = lambda wildcards: ceil(return_memory_requirement(wildcards) * 1.25),
        cores = 16,
        mem_gb = lambda wildcards: return_memory_requirement(wildcards),
        tmpdir = "/tmp",
        bbnorm = 1
    params:
        target = 800,
        error_correction = True,
        extra = lambda wildcards: f"k=31 minkmers=15 mindepth=3 prefilter=t"
    threads: 16
    log: "04-analysis/bbnorm/{sample}_bbnorm.log"
    wrapper:
        "file:///mnt/archgen/users/schumacher/bachelorthesis/02-scripts/py_scripts/wrapper.py"   # "file:///home/alexander_huebner/github/snakemake-wrappers/bio/bbtools/bbnorm"

################################################################################
