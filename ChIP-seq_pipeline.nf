#!/usr/bin/env nextflow

/*
#############################################################
### ChIP-seq pipeline : Trimming - Mapping - Peak calling ###
#############################################################
Script started in January 2019
Authors :
	Coralie Gimonnet <coralie.gimonnet@inra.fr>
This pipeline is adaptated for single-end ChIP-seq data.
*/

/* Default parameters */

params.genome = "genome.fasta"
params.reads = "*.fastq.gz"
params.outdir = "results"
params.configFile = "configFile.txt"
params.script = "pepr_settings.pl"

/* Optional parameters initialization */ 
params.notrim = false
params.help = false
params.fastqc = false
params.stringency = 10
params.broad = false
params.windowSize = false //PePr
params.index = false
params.sharp = false
params.chromSize = false
params.readSize = false
params.genomeSize = false

fasta = file(params.genome)
reads = file(params.reads)
configFile = file(params.configFile)
perl_script = file(params.script)

log.error """\
			ChIP-seq PIPELINE 
		================================================================================
		default parameters :
			genome : ${params.genome}	
			outdir : ${params.outdir}
			configFile : ${params.configFile}
			reads : ${params.reads} (must surrounded with quotes)
	Paired-end data are not supported in this pipeline!
		optional paramaters :
			singleEnd : specifies input are single end reads
			fastqc : make a FASTQC analysis
			notrim : skip trimming step
			index : ${params.index}
			sharp : default = FALSE
			windowSize : ${params.windowSize}
			chromSize : ${params.chromSize} (necessary for broad peak calling)
			script : ${params.script} (necessary to create settings file for PePr using)
			readSize : ${params.readSize}
			genomeSize: ${params.genomeSize}
		Warning : This pipeline assumes all software are on the PATH
		"""
		.stripIndent()

if (params.help) exit 1

/* 
* Create the reads channel
*/

Channel
    .from (reads)
	.ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
	.into { reads_files_fastqc ; read_files_trimming ; reads_mapping } 

/* 
* Create a channel for peak detection file 
* which lists all files needed for peak calling
*/

Channel
	.from (configFile.readLines())
	.ifEmpty { exit 1, "Missing config file for peak calling. Specify with --configFile"}
	.map { line ->
		def list = line.split(',')
		def chip_sample_id = list[0]
		def input_sample_id = list[1]
		def analysis_id = list[2]
		def experiment = list[3]
		[ chip_sample_id, input_sample_id, analysis_id, experiment ]
	}
	.into { file_prepa_pepr ; pepr_para ; epic_para }


/* 
* Step 1. FASTQC Analysis
*/

if (params.fastqc)
{
	process fastqc {
		publishDir "${params.outdir}/fastqc", mode: 'copy',
			saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}
		input:
			file reads from reads_files_fastqc
		output:
			file '*_fastqc.{zip,html}' into fastqc_results
		script:
		"""
			fastqc -q $reads
		"""
	}
}
else
{
	fastqc_results = Channel.from(false)
}


/*
* Step 2. Trim Galore!
*/

if (params.notrim){
	trimmed_reads = read_files_trimming
	trimgalore_results = Channel.from(false)
}
else
{
	if (params.fastqc) 
	{
		process trimming {
			publishDir "${params.outdir}/trim_galore", mode : 'copy',
				saveAs: {filename -> 
					if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
					else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
				}
			input:
				file reads from read_files_trimming
			output:
				file('*.fq.gz') into trimmed_reads
				file "*trimming_report.txt" into trimgalore_results
				file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports
			script:
	        """
	            trim_galore --fastqc --gzip $reads --stringency ${params.stringency}
		    """
		}
	}
	
	else
	{
		process trimming {
			publishDir "${params.outdir}/trim_galore",mode : 'copy',
				saveAs: {filename ->
					if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
				}
			input:
				file(reads) from read_files_trimming
			output:
				file('*.fq.gz') into trimmed_reads
				file "*trimming_report.txt" into trimgalore_results
			script:
			"""
		        trim_galore --gzip $reads --stringency ${params.stringency}
		    """

		}
	} 
}	


/*
* Step 3. Build genome index with Bowtie2
*/
if (params.index)
{
	//bowtie_index = file("${params.index}.fa")
	bowtie_index = Channel
		.fromPath("${params.index}.*.bt2")
		.ifEmpty { exit 1, "Bowtie2 index not found : ${params.index}" }
}
else
{
	process buildIndex {
		cpus 8
		tag "${fasta.baseName}"
		publishDir "${params.outdir}/bowtie-index", mode : 'copy'
		input:
			file fasta from fasta
		output:
			file "${fasta.baseName}.*" into bowtie_index
		script:
		"""
			bowtie2-build ${fasta} ${fasta.baseName} --threads ${task.cpus}
		"""
	}
}


/*
* Step 4. Mapping with Bowtie2
*/

process mapping {
	cpus 24
	tag "${fasta.baseName}"
	publishDir "${params.outdir}/mapping", mode : 'copy',
		saveAs: {filename ->
			if (filename.indexOf(".txt") > 0) "logs/$filename"
		}
	input:
		file index from bowtie_index.collect()
		file(reads) from trimmed_reads
	output:
		file "*.txt" into log_aln
		file "*.bam" into bam_aln
	script: 
	"""
		bowtie2 --end-to-end -p ${task.cpus} --very-sensitive -U $reads -x ${fasta.baseName} | samtools view -bS -@4 - | samtools sort -l9 -@8 -O bam -T tmp - >| ${reads.simpleName}.bwt2.bam
		samtools flagstat ${reads.simpleName}.bwt2.bam > ${reads.simpleName}_flagstat.txt
	"""
}


/*
* Step 5. Deduplication
*/

process deduplicate {
	tag "${bam.simpleName}"
	cpus 8 
	publishDir "${params.outdir}/dedup", mode : 'copy',
		saveAs: {filename ->
			if (filename.indexOf(".txt") > 0) "logs/$filename"
			else if (filename.indexOf(".bam") > 0) "bam/$filename"
			else if (filename.indexOf(".bai") > 0) "bam/$filename"
		}
	input :
		file bam from bam_aln
		//set chip_sample_id, input_sample_id, analysis_id, experiment from file_prepa_pepr
	output:
		file "*.dedup.bam" into bam_dedup, bam_pepr_prep
		file "*_flagstat.dedup.txt" into logs_dedup
		//file "file_pepr_prepare.txt" into file_settings
		file "*.bai" into dedup_bai
	script:
	"""
		samtools sort -@ ${task.cpus} $bam | samtools rmdup -s - ${bam.simpleName}.dedup.bam
		samtools index ${bam.simpleName}.dedup.bam
		samtools flagstat ${bam.simpleName}.dedup.bam > ${bam.simpleName}_flagstat.dedup.txt
	"""
}


/*
* Step 6. Peak detection with PePr (sharp peaks)
*/
if (params.sharp)
{
	process configFilePePr {
		cpus 8
		//publishDir "${params.outdir}", mode : 'copy'
		input:
			file configFile from configFile
			file perl_script from perl_script
		output:
			file "*.txt" into pepr_config
		script:
		"""
			perl $perl_script $configFile ${task.cpus} ${params.windowSize}
		"""
	}
	process PePrDetection {
		tag "${file_pepr_prepare[0].simpleName}"
		publishDir "${params.outdir}/peak_calling_pepr", mode : 'copy'
		input:	
			file file_pepr_prepare from pepr_config.collect()
			file (bam) from bam_dedup.collect()
			set chip_sample_id, input_sample_id, analysis_id, experiment from file_prepa_pepr
		output:
			file "*" into narrowPeakFiles, narrowPeakID
			file "*.txt" into logs_pepr
		script:
			//if (experiment == ${file_pepr_prepare.simpleName})
			//{
			def exp = experiment == '' ? '' : "-p ${file_pepr_prepare[0]}"
			"""
				PePr $exp
			""" 
			//}
	}
}
else
{
	chromSize = Channel
		.fromPath(${params.chromSize})
		.ifEmpty { exit 1, "Chromosome size file not found : ${params.chromSize}" }
	if (!params.readSize)
	{
		println "Please specify a read size parameter for peak detection"
		System.exit(1)
	}
	if (!params.genomeSize)
	{
		println "Please specify a genome size parameter for peak detection"
		System.exit(1)
	}
	process epicDetection {
		cpus 8
		tag "${bam[0].simpleName}"
		publishDir "${params.outdir}/peak_calling_epic", mode : 'copy'
		input:
			file (chrsize) from chromSize
			file (bam) from bam_dedup.collect()
			set chip_sample_id, input_sample_id, analysis_id, experiment from epic_para
		output:
			file "*.bed" into broadPeakFiles, broadPeakID
		script:
			def ctrl = input_sample_id == '' ? '' : "--control ${input_sample_id}.dedup.bam"
			"""
				epic --bam \\
				--chromSize $chrsize \\
				--number-cores ${task.cpus} \\
				$ctrl \\
				--treatment ${chip_sample_id}.dedup.bam \\
				--fragment-size ${params.readSize} \\
				--gaps-allowed 2 \\
				--false-discovery-rate 0.05
				--effective-genome-size ${params.genomeSize} \\
				--bed ${chip_sample_id}.bed
			"""
	}
}

/*
* Step 6. Peak detection (à transformer pour epic car fonctionne pour pepr mais ne gère pas les réplicats)
*/
/*
process peakCalling {
	cpus 16
	tag "${bam[0].simpleName}"
	publishDir "${params.outdir}/peak_calling", mode : 'copy',
		saveAs: {filename ->
			if (filename.indexOf(".bed") > 0) "PePr_call/$filename"
			else if (filename.indexOf(".txt") > 0) "PePr_logs/$filename"
		}
	input:
		file (bam) from bam_dedup.collect()
		file conf from pepr_config
	output:
		file "*.bed" into narrowPeakFiles, narrowPeakID
		file "*.txt" into logs_pepr
	script:
		def ctrl = input_sample_id == '' ? '' : "-i ${input_sample_id}.dedup.bam"
		if (params.windowSize)
		{
		"""
			PePr -p $config \\
			-f $bam \\
			--peaktype=SHARP \\
			-w ${params.windowSize} \\
			--num-processors ${task.cpus}
		"""
		}
		else
		{
		"""
			PePr -c ${chip_sample_id}.dedup.bam \\
			$ctrl \\
			-f bam \\
			--peaktype=SHARP \\
			--num-processors ${task.cpus}
		"""
		}
}
*/


/*
* Step 7. Software version
*/

process softwareVersion {
	publishDir "${params.outdir}/version", mode : 'copy'
	output:
		file 'software_version.txt' into software_version
	script:
	"""
		echo 'FastQC version:' > software_version.txt
		fastqc --version >> software_version.txt
		echo 'Trim Galore! version:' >> software_version.txt
		trim_galore --version >> software_version.txt
		echo 'Cutadapt version:' >> software_version.txt
		cutadapt --version >> software_version.txt
		echo 'bowtie2 version:' >> software_version.txt
		bowtie2 --version >> software_version.txt
		echo 'epic version:' >> software_version.txt
		epic --version >> software_version.txt
		echo 'PePr version:' >> software_version.txt
		PePr --version >> software_version.txt
		echo 'Samtools version:' >> software_version.txt
		samtools --version >> software_version.txt
	"""
}