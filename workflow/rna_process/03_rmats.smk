#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import shutil
import os
from collections import defaultdict
from itertools import combinations

## Read in samplesheet
samples = pd.read_csv(config["samplesheet"], sep='\t')
samples = samples.astype(str)

for c in set(config["mergeBy"] + config["groupBy"]):
    samples[c] = samples[c].fillna("").astype(str).str.replace("\r", "", regex=False).str.strip()

bam_dir = str(config["bamDir"]).replace("\r", "").strip()
bam_suffix = str(config["bamSuffix"]).replace("\r", "").strip()

samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)
samples['group'] = samples[config['groupBy']].agg('_'.join, axis=1)
samples["bam"] = samples["mn"].apply(lambda s: os.path.join(bam_dir, f"{s}{bam_suffix}"))

condition_bams = defaultdict(list)
for _, r in samples.iterrows():
    condition_bams[r["group"]].append(r["bam"])

for g in condition_bams:
    condition_bams[g] = sorted(condition_bams[g])

GROUPS = sorted(condition_bams.keys())

# making pairwise
COMPARISONS = [tuple(x) for x in config["comparison_group"]]
#COMPARISONS = []
#if len(GROUPS) >= 2:
#    COMPARISONS = list(combinations(GROUPS, 2))

rule all:
    input:
        expand("rna_output/rmats/lists/{group}.txt", group=GROUPS),
        expand("rna_output/rmats/{g1}_vs_{g2}/summary.txt",zip,g1=[x[0] for x in COMPARISONS],g2=[x[1] for x in COMPARISONS])

rule sample_list:
    output:
        txt = "rna_output/rmats/lists/{group}.txt"
    params:
        bams = lambda wildcards: condition_bams[wildcards.group]
    run:
        os.makedirs("rna_output/rmats/lists", exist_ok=True)
        with open(output.txt, 'w') as f:
            f.write(','.join(params.bams))

rule rmats:
    input:
        b1 = "rna_output/rmats/lists/{g1}.txt",
        b2 = "rna_output/rmats/lists/{g2}.txt"
    output:
        outdir = "rna_output/rmats/{g1}_vs_{g2}/summary.txt"
    params:
        rmatsVersion = config['rmats_turboVer'],
        gtf = config['gtf'],
        rlength = config['readlength'],
        strand = config['t'],
        outdir = "rna_output/rmats/{g1}_vs_{g2}",
        tmpdir = "rna_output/rmats/{g1}_vs_{g2}/tmp",
        novelss = "--novelSS" if config.get('novelSS', 'FALSE').upper() == 'TRUE' else ""
    threads: 16
    log:
        "rna_output/logs/rmats_{g1}_vs_{g2}.log"
    shell:
        """
        module load rmats-turbo/{params.rmatsVersion}

        mkdir -p {params.outdir} {params.tmpdir}

        python /nas/longleaf/rhel9/apps/rmats-turbo/4.3.0/miniconda3-2/bin/rmats.py \
            --b1 {input.b1} \
            --b2 {input.b2} \
            --gtf {params.gtf} \
            --readLength {params.rlength} \
            --variable-read-length \
            -t {params.strand} \
            --nthread {threads} \
            --od {params.outdir} \
            --tmp {params.tmpdir} \
            --libType fr-unstranded \
            --task both \
            {params.novelss} \
            > {log} 2>&1
        """
