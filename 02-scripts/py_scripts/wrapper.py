__author__ = "Alexander Hübner"
__copyright__ = "Copyright 2022, Alexander Hübner"
__license__ = "MIT"

import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

java_opts = get_java_opts(snakemake)
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Prepare sequencing data
n = len(snakemake.input.fastq)
assert n == 1 or n == 2, "Either one FastQ file (single-end) or two FastQ files (paired-end) are allowed as input."

if n == 1:  # single-end
    in_fqs = f"in={snakemake.input.fastq[0]}"
    out_fqs = f"out={snakemake.output.fastq[0]}"
else: # paired-end
    in_fqs = f"in={snakemake.input.fastq[0]} in2={snakemake.input.fastq[1]}"
    out_fqs = f"out={snakemake.output.fastq[0]} out2={snakemake.output.fastq[1]}"

# Down-sampling
target = snakemake.params.get("target", 0)
assert int(target), "The parameter target must be an integer value."
if target > 0:
    down_sampling = f"target={target} keepall=f"
else:
    down_sampling = ""

# Error correction
error_correction = snakemake.params.get("error_correction", False)
if error_correction:
    err_corr = "ecc=t"
else:
    err_corr = ""

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        f"(bbnorm.sh {java_opts}"
        " {in_fqs} {out_fqs}"
        " {down_sampling} {err_corr}"
        f" tmpdir={tmpdir} threads={snakemake.threads}"
        " {extra} {log})"
    )

