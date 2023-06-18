################################################################################
# Download BAM files from ENA using aspera
#
# Alex Huebner, 22/02/2023
################################################################################

from pathlib import Path
import os
import sys

import pandas as pd
import enasearch

#### ERRS ######################################################################
ERRS = [line.rstrip()
        for line in open(config['samplelist'], "rt")]
################################################################################

#### Auxilliary functions ######################################################


def determine_url(wildcards, slot, protocol):
    rep = pd.read_csv(checkpoints.download_info.get(**wildcards).output[0],
                         sep="\t", index_col=['run_accession'])
    if rep.at[wildcards.err, f'{slot}_{protocol}'] != "":
        return rep.at[wildcards.err, f'{slot}_{protocol}']
    else:
        print(f'There are no files for the column "{slot}_{protocol}" in the ENA database. '
              'Fallback to fastq_ftp.', file=sys.stderr)
        return rep.at[wildcards.err, 'fastq_ftp']


def determine_md5sum(wildcards, slot):
    rep = pd.read_csv(checkpoints.download_info.get(**wildcards).output[0],
                         sep="\t", index_col=['run_accession'])
    if rep.at[wildcards.err, f'{slot}_md5'] != "":
        return rep.at[wildcards.err, f'{slot}_md5']
    else:
        print(f'There are no entries for the column "{slot}_md5" in the ENA database. '
              'Fallback to fastq_md5.', file=sys.stderr)
        return rep.at[wildcards.err, 'fastq_md5']


################################################################################

rule all:
    input:
        [f"{config['outdir']}/{err}.validated" for err in ERRS]

checkpoint download_info:
    output:
        f"{config['outdir']}/ENA_report.tsv"
    message: "Download the info of the run accessions in ENA"
    run:
        reports = []
        for err in ERRS:
            report = pd.DataFrame([line.rstrip().split("\t")
                                   for line in enasearch
                                   .retrieve_filereport(accession=err,
                                                           result="read_run")
                                   .split("\n")])
            report.columns = report.iloc[0,:].tolist()
            report = report.iloc[1:(report.shape[0] - 1)]
            reports.append(report)
        pd.concat(reports).to_csv(output[0], sep="\t", index=False)

rule download_seqdata:
    input:
        f"{config['outdir']}/ENA_report.tsv"
    output:
        touch("{dir}/{err}.downloaded")
    message: "Download sequencing data file {wildcards.err}"
    resources:
        downloads = 1
    params:
        urls = lambda wildcards: " ".join(determine_url(wildcards, config['slot'], config['protocol']).split(";")),
        aspera_path = config['aspera_path']
    shell:
        """
        cd {wildcards.dir}
        for url in {params.urls}; do
            {params.aspera_path}/bin/ascp -i /usr/local64/opt/aspera/connect/etc/asperaweb_id_dsa.openssh \
                -Tr -Q -l 100m -P33001 -L- \
                era-fasp@${{url}} $(basename ${{url}})
        done
        """

rule compute_md5sum:
    input:
        rep = f"{config['outdir']}/ENA_report.tsv",
        download = "{dir}/{err}.downloaded"
    output:
        touch("{dir}/{err}.md5sum")
    message: "Compute md5sum for sequencing data file of {wildcards.err}"
    params:
        fns = lambda wildcards: " ".join([f"{wildcards.dir}/{os.path.basename(url)}" for url in determine_url(wildcards, config['slot'], config['protocol']).split(";")])
    shell:
        """
        for fn in {params.fns}; do
          md5sum ${{fn}} > ${{fn}}.md5
        done
        """

rule validate_md5sum:
    input:
        rep = f"{config['outdir']}/ENA_report.tsv",
        md5 = "{dir}/{err}.md5sum"
    output:
        "{dir}/{err}.validated"
    message: "Validate sequencing data file {wildcards.err}"
    run:
        urls = determine_url(wildcards, config['slot'], config['protocol']).split(";")
        md5sum = [open(f"{wildcards.dir}/{os.path.basename(url)}.md5", "rt").readline().rstrip().split()[0]
                  for url in urls]
        ena_md5sum = determine_md5sum(wildcards, config['slot']).split(";")
        if ena_md5sum == md5sum:
            Path(output[0]).touch()
        else:
            print(f"Downloaded file is different from file on webserver.",
                  file=sys.stderr)
            sys.exit(1)
