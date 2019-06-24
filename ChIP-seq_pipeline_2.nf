#!/usr/bin/env nextflow

/*
##########################################################
### ChIP-seq pipeline : Identification - Union - Count ###
##########################################################
Script started in June 2019
Authors :
	Coralie Gimonnet <coralie.gimonnet@inra.fr>
This pipeline is adaptated for single-end ChIP-seq data.
This script comes after ChIP-seq_pipeline.nf.
*/

/* Default parameters */

params.genome = "genome.fasta"
params.bam = "*.bam"
params.control = "control.bed"
params.treat = "treatment.bed"
params.script = "scripts"
params.project = "project"
params.outdir = "results"
params.mark = "H3K4me3"

params.sharp = false
params.help = false


ctrl = file(params.control)
treat = file(params.treat)
genome = file(params.genome)
script = file(params.script)


log.error """\
			ChIP-seq PIPELINE 2
		===================================================================================
		default parameters :
			genome : ${params.genome}	
			outdir : ${params.outdir}
			bam : ${params.bam}
			control : ${params.control}
			treatment : ${params.treat}
			scripts : ${params.scripts} 
			mark : ${params.mark}
			project : ${params.project}
		optional paramaters :
			sharp : ${params.sharp}
		Warning : This pipeline assumes all software are on the PATH
		"""
		.stripIndent()

if (params.help) exit 1

/*
* Create Channel
*/

Channel
	.fromPath (params.bam)
	.ifEmpty { exit 1, "Cannot find any bam file for count"}
	.set { bamCount }

Channel
	.fromPath (params.control)
	.ifEmpty { exit 1, "Cannot find any bed file for control"}
	.set { bedControl }

Channel
	.fromPath (params.treat)
	.ifEmpty { exit 1, "Cannot find any bed file for treatment"}
	.set { bedTreat }

Channel
	.fromPath (params.script)
	.ifEmpty { exit 1, "Cannot find any scripts for concatenation and union"}
	.into { perl_script ; R_script }


/*
* Step 1. Add a peak ID in bed file
*/

process peakIDControl {
	tag "$bed.simpleName"
	input:
		file (bed) from bedControl
	output:
		file "*.id.bed" into bedIDControl
	shell:
	'''
		cat !{bed} | awk '{$3=$3"\\t""peak_"NR}1' OFS="\\t" > !{bed.simpleName}.id.bed
	'''
}
process peakIDTreat {
	tag "$bed.simpleName"
	input:
		file (bed) from bedTreat
	output:
		file "*.id.bed" into bedIDTreat
	shell:
	'''
		cat !{bed} | awk '{$3=$3"\\t""peak_"NR}1' OFS="\\t" > !{bed.simpleName}.id.bed
	'''
}

/*
* Step 2. Concatenation of the two bed files
*/

process concatenation {
	input:
		file (treatbed) from bedIDTreat
		file (controlbed) from bedIDControl
		file script from perl_script
	output:
		file ("*by_condition.bed") into bedConcat
	script:
	"""
		perl $script/peak_concatenation.pl $controlbed $treatbed ${params.project}_by_condition.bed
	"""
}

/*
* Step 3. Union between control and treatment
*/
process union {
	input:
		file (bed) from bedConcat
		file script from R_script.collect()
	output:
		file ('*.saf') into safFile
	script:
		if (params.sharp)
		{
		"""
			Rscript $script/peaks_union.R --input $bed --experiment ${params.project} --mark ${params.mark} --sharp TRUE
		"""
		}
		else
		{
		"""
			Rscript $script/peaks_union.R --input $bed --experiment ${params.project} --mark ${params.mark} --broad TRUE
		"""
		}
}

/*
* Step 4. Count of reads on coordinates
*/
process count {
	cpus 8
	publishDir "${params.outdir}/featureCounts", mode: 'copy'
	input :
		file ('*.bam') from bamCount.toList()
		file (saf) from safFile
	output :
		file ('*.counts') into countsFinal
		file ('*.counts.summary') into countsSummary
	script :
	"""
		featureCounts -T ${task.cpus} -a $saf -G ${params.genome} -o ${params.project}.counts -F SAF *.bam
	"""
}


/*
* Step 5. Software version
*/

process softwareVersion {
	publishDir "${params.outdir}/version", mode : 'copy'
	output:
		file 'software_version.txt' into software_version
	script:
	"""
		echo 'R version:' > software_version.txt
		R --version >> software_version.txt
		echo 'Perl version:' >> software_version.txt
		perl --version >> software_version.txt
		echo 'featureCounts version:' >> software_version.txt
		featureCounts -v >> software_version.txt
	"""
} 