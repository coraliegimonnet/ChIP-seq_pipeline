---
title: "README - ChIP-seq_pipeline"
author: "Coralie Gimonnet"
output:
  html_document:
    toc: true
    number_sections: true
    highlight: textmate
    theme: cerulean
    fig_caption: false
date: "June 2019"
---

This pipeline was developed for single-end fastq reads from ChIP-seq experiment.

This pipeline is based on different tools :

- [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for fastq analysis
- [Trim Galore](https://github.com/FelixKrueger/TrimGalore) for trimming (
Warning: Trim galore is used as recommended in manual for WGBS, so stringency to cut adapter was fixed at 1bp)
- [PePr](https://github.com/shawnzhangyx/PePr) for narrow peak calling
- [epic](https://github.com/biocore-ntnu/epic) for broad peak calling

**Warnings**: 
This pipeline assumes all software are installed on PATH
Paired-end data are not supported in this pipeline!


It is writen with [Nextflow framework](https://www.nextflow.io/).

Print help of ChIP-seq_pipeline.nf :

	nextflow run ChIP-seq_pipeline.nf --help

Parameters :

	--genome:	fasta file of reference genome
	--reads:	single-end reads need to be write like this 'toto.fq.gz' (surrounded by quotes)
	--outdir:	directory where all results are stored
	--configFile:	a configuration file which is necessary for peak calling

Optionnal parameters :

	--script:	path to perl script (necessary for narrow peak calling with PePr)
	--fastqc:	make a FastQC analysis
	--notrim:	skip trimming step
	--index:	path to index genome files (surrounded by quotes)
	--sharp;	for narrow peak detection (default : FALSE)
	--windowSize:	fixe a window size for PePr peak detection (default : automatically fixed by PePr)
	--chromSize:	path to chromosome size of reference genome (necessary for broad peak calling)
	--readSize:	length of reads (necessary for broad peak calling)
	--genomeSize:	size of genome in Mb



The configuration file needs to be writen like this :

	C1_H3K4me3_subsample.fastq.gz,C1_INPUT_subsample.fastq.gz,C1,control
	T1_H3K4me3_subsample.fastq.gz,T1_INPUT_subsample.fastq.gz,T1,treatment
	C2_H3K4me3_subsample.fastq.gz,C2_INPUT_subsample.fastq.gz,C2,control
	T2_H3K4me3_subsample.fastq.gz,T2_INPUT_subsample.fastq.gz,T2,treatment


Example of usage:

	nextflow run ChIP-seq_pipeline.nf --genome genome.fasta --outdir test_mai --configFile configFile.txt --reads '*fastq.gz' --script pepr_settings.pl --sharp --index 'bowtie-index/*' 


If an incident happens, you could rerun your command line with `-resume` option.