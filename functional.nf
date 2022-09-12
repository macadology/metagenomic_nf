#!/usr/bin/env nextflow

//=========== DSL 2 syntax
nextflow.preview.dsl=2

//=========== Help message
params.help = false
def helpMessage() {
    log.info"""
    ############################################################################
    Jon's functional profile nextflow pipeline
    ----------------------------------------------------------------------------
    Usage:
    The typical command for running the pipeline is as follows:
        nextflow run main.nf
    Input arguments:
        --readtype              (raw) raw reads
                                (fastp) fastp reads *DEFAULT*
                                (decont) human decontaminated fastp reads (custom) reads from custom directory and filename. User must define --querydir and --queryglob.
        --querydir              Path to a folder containing all input fastq.
                                Required if readtype is (custom)
        --queryglob             Glob pattern of paired reads.
                                Required if readtype is (custom)
                                Optional for others.
                                E.g.: "*_{1,2}*{fastq,fastq.gz,fq,fq.gz}"
        --profilers             functional profilers to run. If not defined,
                                show matched reads. Options ->
        --fastpreads            Whether the default reads should be referenced
                                from the fastp folder (true) or raw folder (false). Default: true
        --overwrite             Overwrite existing directories. Otherwise,
                                skip. Default: false
    Other arguments:
    Edit params in nextflow.config or local.config where appropriate.

    Examples:
	nextflow run functional.nf #show the reads that'll be used as inputs
	nextflow run functional.nf -w <WORKDIR> # specify where to store work directory
    nextflow run functional.nf --readtype custom --querydir [Directory] --queryglob [glob pattern] #show reads queryglob and directory
    nextflow run functional.nf --readtype raw --profilers fastp,decont
    nextflow run functional.nf --readtype fastp --profilers decont
    nextflow run functional.nf --readtype raw --queryglob [glob pattern] --profilers fastp
    nextflow run functional.nf --krakenMMAP --profilers kraken2,bracken --krakenKeepOutput true
    nextflow run functional.nf --readtype decont --krakenMMAP --profilers kraken2,bracken
    nextflow run functional.nf --profilers fastp,kraken2,bracken
    nextflow run functional.nf --profilers srst2
    nextflow run functional.nf --profilers align --bwaIndex ref.fasta
    nextflow run functional.nf --profilers align --bwaIndexDir IndexDir --bwaIndex ref.fasta
    ############################################################################
    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}

//=========== Profilers
Set profilers_expected = ['bowtie2', 'diamond', 'annotate']
Set profilers = []
if(params.profilers.getClass() != Boolean){
    Set profilers_input = params.profilers.split(',')
    Set profiler_diff = profilers_input - profilers_expected
    profilers = profilers_input.intersect(profilers_expected)
    if( profiler_diff.size() != 0 ) {
    	log.warn "[Pipeline warning] Profiler $profiler_diff is not supported yet! Will only run $profilers.\n"
    }
}

println profilers

//=========== Parameters
