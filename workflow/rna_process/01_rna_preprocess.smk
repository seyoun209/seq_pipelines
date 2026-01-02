#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import shutil
import os
from utils.namer import namer

## Read in samplesheet
samples = pd.read_csv(config["samplesheet"], sep='\t')
samples = samples.astype(str)

## Concatenate the sequencing directory to Read1 and Read2 for full paths
samples['Read1'] = samples[['Sequencing_Directory', 'Read1']].apply(lambda row: os.path.join(*row), axis=1)
samples['Read2'] = samples[['Sequencing_Directory', 'Read2']].apply(lambda row: os.path.join(*row), axis=1)

samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)

read1 = samples.groupby('mn')['Read1'].apply(list).to_dict()
read2 = samples.groupby('mn')['Read2'].apply(list).to_dict()
runName = namer(samples, config['mergeBy'])

## Success message
onsuccess:
    print("RNA_preprocess completed successfully! Wahoo!")

##### Define rules #####
rule all:
    input:
    # Raw FastQC outputs
        expand("rna_output/QC/raw/{sampleName}_{read}_fastqc.{ext}",sampleName=read1.keys(), read=['R1', 'R2'], ext=['zip', 'html']),
    # Trimming outputs
        expand('rna_output/trim/{sampleName}_R{read}_val_{read}.fq.gz',sampleName=read1.keys(), read=['1', '2']),
    # Trimmed FastQC outputs
        expand("rna_output/QC/trimmed/{sampleName}_R{read}_val_{read}_fastqc.{ext}",sampleName=read1.keys(), read=['1', '2'], ext=['zip', 'html']),
    # MultiQC reports
        "rna_output/QC/raw/multiqc_report.html",
        "rna_output/QC/trimmed/multiqc_report.html",
    # Alignment outputs (genome and transcriptome)
        expand("rna_output/align/{sampleName}_Aligned.out.bam",sampleName=read1.keys()),
        expand("rna_output/align/{sampleName}_Aligned.toTranscriptome.out.bam",sampleName=read1.keys()),
    # Sorted/indexed BAMs
        expand("rna_output/align/{sampleName}_sorted.bam",sampleName=read1.keys()),
        expand("rna_output/align/{sampleName}_sorted.bai",sampleName=read1.keys()),
    # Flagstat
        expand("rna_output/align/{sampleName}_sorted.flagstat",sampleName=read1.keys()),
    # Salmon quantification
        expand('rna_output/quant/{sampleName}/quant.sf',sampleName=read1.keys()),
    # Qualimap reports
        expand("rna_output/QC/qualimap/{sampleName}_bamqc/qualimapReport.html",sampleName=read1.keys()),
        expand("rna_output/QC/qualimap/{sampleName}_rnaseq/qualimapReport.html",sampleName=read1.keys())
#    # VerifyBamID2 outputs
#    expand("rna_output/verifyBamID/{sampleName}.selfSM", sampleName=read1.keys()),
#    # Signal tracks
#    expand("rna_output/signals/no_norm_indv/{sampleName}.bw",sampleName=read1.keys()),
#   expand("rna_output/signals/norm_indv/{sampleName}.bw",sampleName=read1.keys())


rule catReads:
    output:
        R1 = temp('rna_output/fastq/{sampleName}_R1.fastq.gz'),
        R2 = temp('rna_output/fastq/{sampleName}_R2.fastq.gz')
    params:
        r1 = lambda wildcards: read1[wildcards.sampleName],
        r2 = lambda wildcards: read2[wildcards.sampleName]
    log:
        err = 'rna_output/logs/catReads_{sampleName}.err',
        out = 'rna_output/logs/catReads_{sampleName}.out'
    shell:
        """
        mkdir -p rna_output/fastq
        cat {params.r1} > {output.R1} 2> {log.err}
        cat {params.r2} > {output.R2} 2>> {log.err}
        """

rule fastqc_raw:
    input:
        R1 = rules.catReads.output.R1,
        R2 = rules.catReads.output.R2
    output:
        zip1 = "rna_output/QC/raw/{sampleName}_R1_fastqc.zip",
        zip2 = "rna_output/QC/raw/{sampleName}_R2_fastqc.zip",
        html1 = "rna_output/QC/raw/{sampleName}_R1_fastqc.html",
        html2 = "rna_output/QC/raw/{sampleName}_R2_fastqc.html"
    log:
        err = 'rna_output/logs/fastqc_raw_{sampleName}.err',
        out = 'rna_output/logs/fastqc_raw_{sampleName}.out'
    params:
        dir = "rna_output/QC/raw",
        version = config['fastqcVers']
    benchmark:
        'rna_output/benchmarks/fastqc_raw_{sampleName}.tsv'
    shell:
        """
        module load fastqc/{params.version}
        mkdir -p {params.dir}
        fastqc -o {params.dir} {input.R1} {input.R2} 1> {log.out} 2> {log.err}
        """
rule trim:
    input:
        R1 = rules.catReads.output.R1,
        R2 = rules.catReads.output.R2
    output:
        trim1 = 'rna_output/trim/{sampleName}_R1_val_1.fq.gz',
        trim2 = 'rna_output/trim/{sampleName}_R2_val_2.fq.gz',
        report1 = 'rna_output/trim/{sampleName}_R1.fastq.gz_trimming_report.txt',
        report2 = 'rna_output/trim/{sampleName}_R2.fastq.gz_trimming_report.txt'
    threads: 4
    params:
        trim_galore_ver = config['trim_galore'],
        python_ver = config['pythonVers'],
        pigz_ver = config['pgizVers']
    log:
        err = 'rna_output/logs/trim_{sampleName}.err',
        out = 'rna_output/logs/trim_{sampleName}.out'
    benchmark:
        'rna_output/benchmarks/trim_{sampleName}.tsv'
    shell:
        """
        module load trim_galore/{params.trim_galore_ver}
        module load python/{params.python_ver}
        module load pigz/{params.pigz_ver}
        mkdir -p rna_output/trim
        trim_galore -o rna_output/trim --cores {threads} --paired {input.R1} {input.R2} 2> {log.err}
        """

rule fastqc_trimmed:
    input:
        trim1 = rules.trim.output.trim1,
        trim2 = rules.trim.output.trim2
    output:
        zip1 = "rna_output/QC/trimmed/{sampleName}_R1_val_1_fastqc.zip",
        zip2 = "rna_output/QC/trimmed/{sampleName}_R2_val_2_fastqc.zip",
        html1 = "rna_output/QC/trimmed/{sampleName}_R1_val_1_fastqc.html",
        html2 = "rna_output/QC/trimmed/{sampleName}_R2_val_2_fastqc.html"
    log:
        err = 'rna_output/logs/fastqc_trimmed_{sampleName}.err',
        out = 'rna_output/logs/fastqc_trimmed_{sampleName}.out'
    params:
        dir = "rna_output/QC/trimmed",
        version = config['fastqcVers']
    benchmark:
        'rna_output/benchmarks/fastqc_trimmed_{sampleName}.tsv'
    shell:
        """
        module load fastqc/{params.version}
        mkdir -p {params.dir}
        fastqc -o {params.dir} {input.trim1} {input.trim2} 1> {log.out} 2> {log.err}
        """

rule multiqc_raw:
    input:
        expand("rna_output/QC/raw/{sampleName}_{read}_fastqc.zip",
               sampleName=read1.keys(), read=['R1', 'R2'])
    output:
        "rna_output/QC/raw/multiqc_report.html"
    params:
        indir = "rna_output/QC/raw",
        outdir = "rna_output/QC/raw",
        multiqc_ver = config['multiqcVers']
    log:
        err = 'rna_output/logs/multiqc_raw.err',
        out = 'rna_output/logs/multiqc_raw.out'
    shell:
        """
        module load multiqc/{params.multiqc_ver}
        multiqc {params.indir} -o {params.outdir} 1> {log.out} 2> {log.err}
        """

rule multiqc_trimmed:
    input:
        expand("rna_output/QC/trimmed/{sampleName}_R{read}_val_{read}_fastqc.zip",
               sampleName=read1.keys(), read=['1', '2'])
    output:
        "rna_output/QC/trimmed/multiqc_report.html"
    params:
        indir = "rna_output/QC/trimmed",
        outdir = "rna_output/QC/trimmed",
        multiqc_ver = config['multiqcVers']
    log:
        err = 'rna_output/logs/multiqc_trimmed.err',
        out = 'rna_output/logs/multiqc_trimmed.out'
    shell:
        """
        module load multiqc/{params.multiqc_ver}
        multiqc {params.indir} -o {params.outdir} 1> {log.out} 2> {log.err}
        """
rule align:
    input:
        R1 = rules.trim.output.trim1,
        R2 = rules.trim.output.trim2
    output:
        genome_bam = "rna_output/align/{sampleName}_Aligned.out.bam",
        transcriptome_bam = "rna_output/align/{sampleName}_Aligned.toTranscriptome.out.bam",
        log_out = "rna_output/align/{sampleName}_Log.out",
        progress_log = "rna_output/align/{sampleName}_Log.progress.out",
        final_log = "rna_output/align/{sampleName}_Log.final.out",
        sj_tab = "rna_output/align/{sampleName}_SJ.out.tab"
    threads: 8
    log:
        err = 'rna_output/logs/align_{sampleName}.err',
        out = 'rna_output/logs/align_{sampleName}.out'
    params:
        index = config['star'],
        gtf = config['gtf'],
        starVer = config['starVers'],
        dir_align = "rna_output/align/{sampleName}_"
    benchmark:
        'rna_output/benchmarks/align_{sampleName}.tsv'
    shell:
        """
        ml star/{params.starVer}
        mkdir -p rna_output/align
        
        STAR --genomeDir {params.index} \
                --runThreadN {threads} \
                --sjdbGTFfile {params.gtf} \
                --outFileNamePrefix {params.dir_align} \
                --readFilesCommand zcat \
                --outSAMtype BAM Unsorted \
                --readFilesIn {input.R1} {input.R2} \
                --outFilterType BySJout \
                --alignSJoverhangMin 8 \
                --outFilterMultimapNmax 20 \
                --alignSJDBoverhangMin 1 \
                --outFilterMismatchNmax 999 \
                --outFilterMismatchNoverReadLmax 0.04 \
                --alignIntronMin 20 \
                --alignIntronMax 1000000 \
                --alignMatesGapMax 1000000 \
                --quantMode TranscriptomeSAM \
                --outSAMattributes NH HI AS NM MD 1> {log.out} 2> {log.err}
        """
rule sort_index:
    input:
        rules.align.output.genome_bam
    output:
        bam = 'rna_output/align/{sampleName}_sorted.bam',
        bai = 'rna_output/align/{sampleName}_sorted.bai'
    threads: 8
    params:
        samtoolsVersion = config['samtoolsVers']
    log:
        out = 'rna_output/logs/sort_index_{sampleName}.out',
        err = 'rna_output/logs/sort_index_{sampleName}.err'
    benchmark:
        'rna_output/benchmarks/sort_index_{sampleName}.tsv'
    shell:
        """
        module load samtools/{params.samtoolsVersion}
        samtools sort -@ {threads} {input} -o {output.bam} 2> {log.err}
        samtools index -@ {threads} {output.bam} 2>> {log.err}
        """

rule flagstat:
    input:
        rules.sort_index.output.bam
    output:
        "rna_output/align/{sampleName}_sorted.flagstat"
    params:
        samtoolsVersion = config['samtoolsVers']
    log:
        err = 'rna_output/logs/flagstat_{sampleName}.err'
    shell:
        """
        module load samtools/{params.samtoolsVersion}
        samtools flagstat {input} > {output} 2> {log.err}
        """

rule qualimap_bamqc:
    input:
        bam = rules.sort_index.output.bam,
        bai = rules.sort_index.output.bai
    output:
        "rna_output/QC/qualimap/{sampleName}_bamqc/qualimapReport.html"
    params:
        gtf = config['gtf'],
        outdir = "rna_output/QC/qualimap/{sampleName}_bamqc",
        qualimap_ver = config['qualimapVers'],
        java_mem = config.get('java_mem', '16G')
    threads: 4
    log:
        err = 'rna_output/logs/qualimap_bamqc_{sampleName}.err',
        out = 'rna_output/logs/qualimap_bamqc_{sampleName}.out'
    benchmark:
        'rna_output/benchmarks/qualimap_bamqc_{sampleName}.tsv'
    shell:
        """
        module load qualimap/{params.qualimap_ver}
        mkdir -p {params.outdir}
        qualimap bamqc -bam {input.bam} -gff {params.gtf} \
            -outdir {params.outdir} --java-mem-size={params.java_mem} \
            1> {log.out} 2> {log.err}
        """

rule qualimap_rnaseq:
    input:
        bam = rules.sort_index.output.bam,
        bai = rules.sort_index.output.bai
    output:
        "rna_output/QC/qualimap/{sampleName}_rnaseq/qualimapReport.html"
    params:
        gtf = config['gtf'],
        outdir = "rna_output/QC/qualimap/{sampleName}_rnaseq",
        qualimap_ver = config['qualimapVers'],
        java_mem = config.get('java_mem', '16G')
    threads: 4
    log:
        err = 'rna_output/logs/qualimap_rnaseq_{sampleName}.err',
        out = 'rna_output/logs/qualimap_rnaseq_{sampleName}.out'
    benchmark:
        'rna_output/benchmarks/qualimap_rnaseq_{sampleName}.tsv'
    shell:
        """
        module load qualimap/{params.qualimap_ver}
        mkdir -p {params.outdir}
        qualimap rnaseq -bam {input.bam} -gtf {params.gtf} \
            -outdir {params.outdir} --java-mem-size={params.java_mem} \
            1> {log.out} 2> {log.err}
        """

rule quant:
    input:
        transcriptome_bam = rules.align.output.transcriptome_bam
    output:
        quant = 'rna_output/quant/{sampleName}/quant.sf',
        quant_genes = 'rna_output/quant/{sampleName}/quant.genes.sf'
    params:
        salmonVer = config['salmonVers'],
        index = config['salmon'],
        outdir = 'rna_output/quant/{sampleName}',
        lib_type = config.get('salmon_lib_type', 'A')
    log:
        out = 'rna_output/logs/salmon_{sampleName}.out',
        err = 'rna_output/logs/salmon_{sampleName}.err'
    threads: 8
    benchmark:
        'rna_output/benchmarks/salmon_{sampleName}.tsv'
    shell:
        """
        ml salmon/{params.salmonVer}
        mkdir -p rna_output/quant
        salmon quant -t {params.index}/gentrome.fa -l {params.lib_type} \
            -a {input.transcriptome_bam} -o {params.outdir} \
            -p {threads} --seqBias --gcBias 1> {log.out} 2> {log.err}
        """
rule add_read_groups:
    input:
        bam = rules.sort_index.output.bam,
        bai = rules.sort_index.output.bai
    output:
        bam = temp("rna_output/align/{sampleName}_RG.bam"),
        bai = temp("rna_output/align/{sampleName}_RG.bam.bai")
    params:
        picard_version = config["picardVers"],
        samtools_version = config["samtoolsVers"]
    threads: 4
    log:
        err = "rna_output/logs/add_read_groups_{sampleName}.err",
        out = "rna_output/logs/add_read_groups_{sampleName}.out"
    benchmark:
        'rna_output/benchmarks/add_read_groups_{sampleName}.tsv'
    shell:
        """
        module load picard/{params.picard_version}
        module load samtools/{params.samtools_version}

        picard AddOrReplaceReadGroups \
            I={input.bam} \
            O={output.bam} \
            RGSM={wildcards.sampleName} \
            RGPL=ILLUMINA \
            RGLB=lib1 \
            RGPU=unit1 \
            1> {log.out} 2> {log.err}

        samtools index {output.bam} 1>> {log.out} 2>> {log.err}
        """
