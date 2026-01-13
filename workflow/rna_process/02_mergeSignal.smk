#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import shutil
import os
from collections import defaultdict

## Read in samplesheet
samples = pd.read_csv(config["samplesheet"], sep='\t')
# samples = pd.read_csv("samplesheet.txt",sep='\t')
samples = samples.astype(str)

for c in set(config["mergeBy"] + config["groupBy"]):
    samples[c] = samples[c].fillna("").astype(str).str.replace("\r", "", regex=False).str.strip()

bam_dir = str(config["bamDir"]).replace("\r", "").strip()
bam_suffix = str(config["bamSuffix"]).replace("\r", "").strip()

samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)
#mergeBy = ["Proj","Donor","Condition","Treatment"]
#samples['mn'] = samples[mergeBy].agg('_'.join, axis=1)

samples['group'] = samples[config['groupBy']].agg('_'.join, axis=1)
#groupBy = ["Condition","Treatment"]
#samples['group'] = samples[groupBy].agg('_'.join, axis=1)

samples["bam"] = samples["mn"].apply(lambda s: os.path.join(bam_dir, f"{s}{bam_suffix}"))

condition_bams = defaultdict(list)
for _, r in samples.iterrows():
    condition_bams[r["group"]].append(r["bam"])

for g in condition_bams:
    condition_bams[g] = sorted(condition_bams[g])

GROUPS = sorted(condition_bams.keys())

rule all:
    input:
        expand("rna_output/align/merged/{group}_merged.bam.bai", group=GROUPS),
        expand("rna_output/align/merged/{group}_merged.bam", group=GROUPS),
        expand("rna_output/align/merged/{group}_merged_flagstat.txt", group=GROUPS),
        expand("rna_output/signals/merged/{group}.bw", group=GROUPS),
        expand("rna_output/signals/merged_norm/{group}.bw", group=GROUPS),
        expand("rna_output/signals/merged_strand/{group}.fwd.bw", group=GROUPS) if config.get("make_strand_bw", False) else [],
        expand("rna_output/signals/merged_strand/{group}.rev.bw", group=GROUPS) if config.get("make_strand_bw", False) else []

rule mergeAlign:
    input:
        bams=lambda wildcards: condition_bams[wildcards.group]
    output:
        merged_bam="rna_output/align/merged/{group}_merged.bam",
        merged_bai="rna_output/align/merged/{group}_merged.bam.bai",
        stats="rna_output/align/merged/{group}_merged_flagstat.txt"
    params:
        samtools_version=config["samtoolsVers"]
    threads: 8
    log:
        err="rna_output/logs/mergeAlign_{group}.err",
        out="rna_output/logs/mergeAlign_{group}.out"
    shell:
        """
        module load samtools/{params.samtools_version}
        mkdir -p rna_output/align/merged

        samtools merge -@ {threads} -f {output.merged_bam} {input.bams} 1>> {log.out} 2>> {log.err}
        samtools index -@ {threads} {output.merged_bam} 1>> {log.out} 2>> {log.err}
        samtools flagstat {output.merged_bam} > {output.stats} 2>> {log.err}
        """

rule mergeSignal:
    input:
        bam = rules.mergeAlign.output.merged_bam,
        bai = rules.mergeAlign.output.merged_bai
    output:
        bw = "rna_output/signals/merged/{group}.bw"
    threads: 4
    params:
        deeptools_ver = config["deeptoolsVers"],
        binSize = config['bin_size']
    log:
        err = "rna_output/logs/mergeSignal_{group}.err"
    benchmark:
        "rna_output/benchmarks/mergeSignal_{group}.tsv"
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p rna_output/signals/merged
        
        bamCoverage \
          --bam {input.bam} \
          --binSize {params.binSize} \
          --samFlagExclude 256 \
          --numberOfProcessors {threads} \
          -o {output.bw} > {log.err} 2>&1
        """

rule mergeSignal_norm:
    input:
        bam = rules.mergeAlign.output.merged_bam,
        bai = rules.mergeAlign.output.merged_bai
    output:
        bw = "rna_output/signals/merged_norm/{group}.bw"
    threads: 4
    params:
        deeptools_ver = config["deeptoolsVers"],
        binSize = config['bin_size'],
        NormOption=config['normalize_option']
    log:
        err = "rna_output/logs/mergeSignal_norm_{group}.err"
    benchmark:
        "rna_output/benchmarks/mergeSignal_norm_{group}.tsv"
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p rna_output/signals/merged_norm

        bamCoverage --normalizeUsing {params.NormOption} \
          --bam {input.bam} \
          --binSize {params.binSize} \
          --samFlagExclude 256 \
          --numberOfProcessors {threads} \
          -o {output.bw} > {log.err} 2>&1
        """

rule mergeForwardSignal:
    input:
        bam = rules.mergeAlign.output.merged_bam,
        bai = rules.mergeAlign.output.merged_bai
    output:
        bw = "rna_output/signals/merged_strand/{group}.fwd.bw"
    threads: 4
    params:
        deeptools_ver = config["deeptoolsVers"],
        binSize = config['bin_size'],
        NormOption=config['normalize_option']
    log:
        err = "rna_output/logs/mergeSignal_fwd_{group}.err"
    benchmark:
        "rna_output/benchmarks/mergeSignal_fwd_{group}.tsv"
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p rna_output/signals/merged_strand

        bamCoverage --normalizeUsing {params.NormOption} \
          --bam {input.bam} \
          --filterRNAstrand forward \
          --binSize {params.binSize} \
          --samFlagExclude 256 \
          --numberOfProcessors {threads} \
          -o {output.bw} > {log.err} 2>&1
        """

rule mergeReverseSignal:
    input:
        bam = rules.mergeAlign.output.merged_bam,
        bai = rules.mergeAlign.output.merged_bai
    output:
        bw = "rna_output/signals/merged_strand/{group}.rev.bw"
    threads: 4
    params:
        deeptools_ver = config["deeptoolsVers"],
        binSize = config['bin_size'],
        NormOption=config['normalize_option']
    log:
        err = "rna_output/logs/mergeSignal_rev_{group}.err"
    benchmark:
        "rna_output/benchmarks/mergeSignal_rev_{group}.tsv"
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p rna_output/signals/merged_strand

        bamCoverage --normalizeUsing {params.NormOption} \
          --bam {input.bam} \
          --filterRNAstrand reverse \
          --binSize {params.binSize} \
          --samFlagExclude 256 \
          --numberOfProcessors {threads} \
          -o {output.bw} > {log.err} 2>&1
        """


