<img src="/home/cgimonnet/Documents/pipeline/boa.png" align="right" width="100" height="34"/>
<img src="/home/cgimonnet/Documents/pipeline/INRA_logo.jpg" align="right" width="300" height="54"/>

# ChIP-seq Pipeline

This pipeline was developed for single-end fastq reads from ChIP-seq experiment. 

This pipeline is based on different tools :

- [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for fastq analysis
- [Trim Galore](https://github.com/FelixKrueger/TrimGalore) for trimming (
Warning: Trim galore is used as recommended in manual for WGBS, so stringency to cut adapter was fixed at 1bp)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for alignment
- [PePr](https://github.com/shawnzhangyx/PePr) for narrow peak calling
- [epic](https://github.com/biocore-ntnu/epic) for broad peak calling
- [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) for count

Previous analysis shows PePr is better to analyse narrow peaks whereas epic is better to analyse broad peaks. 
By default, broad peak detection with epic is made.

**Warnings**: 
 This pipeline assumes all software are installed on PATH
Paired-end data are not supported in this pipeline!
Docker is not available for this pipeline for the moment.

It is writen with [Nextflow framework](https://www.nextflow.io/).

This pipeline need to be run in two step :
- 1 : ChIP-seq_pipeline.nf for quality, trimming, alignment and peak calling
- 2 : ChIP-seq_pipeline_2.nf for union between control and treatment peaks and counting

## Quality, trimming, alignment and peak calling

Peak calling is made with two different tools :

- epic : a tool adaptated for broad peaks [default mode]
- PePr ; a tool adaptated for narrow peaks

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

	nextflow run ChIP-seq_pipeline.nf --genome 'genome.fasta' --outdir 'results' --configFile configFile.txt --reads '*fastq.gz' --script 'scripts/' --sharp --index 'bowtie-index/*' 


If an incident happens, you could rerun your command line with `-resume` option.

# Union and count

Union between control and treatment peaks is made to have all specific peaks and common peaks between this two conditions. Next, reads in bam file are counting on peak coordonates (bed file).

Print help of ChIP-seq_pipeline_2.nf :

	nextflow run ChIP-seq_pipeline_2.nf --help

Parameters :

	--genome:	fasta file of reference genome
	--bam:	bam files produce by ChIP-seq_pipeline.nf
	--control:	bed file of control sample produce by ChIP-seq_pipeline.nf
	--treatment:	bed file of treatment sample produce by ChIP-seq_pipeline.nf
	--scripts:	path to scripts necessary for analysis

Optional parameters:

	--sharp:	setting to specify mark is narrow [default: FALSE]
	--project:	project name
	--outdir:	directory where results are stored
	--mark:	mark of interest [default: H3K27me3]

Example of usage:

	nextflow run ChIP-seq_pipeline_2.nf --genome 'genome.fasta' --outdir 'results' --sharp --control 'controle.bed' --treatment 'treatment.bed' --bam '*.bam' --script 'scripts/'

If an incident happens, you could rerun your command line with `-resume` option.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Authors

* **Coralie Gimonnet** 

